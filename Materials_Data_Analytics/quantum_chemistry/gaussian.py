import pandas as pd
import numpy as np
import DateTime as dt
from Materials_Data_Analytics.laws_and_constants import lorentzian
from Materials_Data_Analytics.core.coordinate_transformer import PdbParser
import plotly.express as px

pd.set_option('mode.chained_assignment', None)


class GaussianParser:
    """
    Class to parse information from a gaussian log file
    """

    def __init__(self, log_file: str | list[str]):

        self._log_file = log_file
        self._lines, self._restart, self._time_stamp = self._concatenate_log_files(log_file)
        self._keywords = self._get_keywords()

        # extract boolean attributes from keywords
        self._raman = any('raman' in s.lower() for s in self._keywords)
        self._freq = any('freq' in s.lower() for s in self._keywords)
        self._opt = any('opt' in s.lower() for s in self._keywords)
        self._stable = any('stable' in s.lower() for s in self._keywords)
        self._pop = any('pop' in s.lower() for s in self._keywords)
        self._solvent = any('scrf' in s.lower() for s in self._keywords)

        # extract boolean attributes from the log file
        self._complete = any('Normal termination of' in s for s in self._lines[-20:])
        self._esp = any('esp charges' in s.lower() for s in self._lines)

        # extract non-boolean attributes from the keywords
        self._functional = [k for k in self._keywords if "/" in k][0].split("/")[0].upper()
        self._basis = [k for k in self._keywords if "/" in k][0].split("/")[1]
        self._n_alpha = int([i for i in self._lines if "alpha electrons" in i][0].split()[0])
        self._n_beta = int([i for i in self._lines if "alpha electrons" in i][0].split()[3])
        self._n_electrons = self._n_alpha + self._n_beta

        # get the charge and multiplicity from the log file
        if any('Charge =' in c for c in self._lines):
            self._charge = int([c for c in self._lines if 'Charge =' in c][0][9:].split()[0])
            self._multiplicity = int([m for m in self._lines if 'Charge =' in m][0][27:])
        else:
            self._charge = None
            self._multiplicity = None

        # get the energy from the log file
        if any('SCF Done' in e for e in self._lines):
            self._energy = float([e for e in self._lines if 'SCF Done' in e][-1].split()[4]) * 2625.5
            self._unrestricted = True if [e for e in self._lines if 'SCF Done' in e][-1].split()[2][2] == "U" else False
            self._scf_iterations = len([e for e in self._lines if 'SCF Done' in e]) - len([e for e in self._lines if '>>>>>>>>>> Convergence criterion not met' in e])
        else:
            self._energy = None
            self._unrestricted = None
            self._scf_iterations = None

        # get the atom counts from the log file
        if any('Mulliken charges' in a for a in self._lines):
            _mull_start = self._lines.index([k for k in self._lines if 'Mulliken charges' in k][0]) + 2
            _mull_end = self._lines.index([k for k in self._lines if 'Sum of Mulliken charges' in k][0])
            self._atomcount = _mull_end - _mull_start
            self._atoms = [a.split()[1] for a in self._lines[_mull_start:_mull_end]]
            self._heavyatoms = [a.split()[1] for a in self._lines[_mull_start:_mull_end] if 'H' not in a]
            self._heavyatomcount = len(self._heavyatoms)
        else:
            self._atomcount = None
            self._atoms = None
            self._heavyatoms = None
            self._heavyatomcount = None

        # Get the stability report from the log file
        if " The wavefunction is stable under the perturbations considered.\n" in self._lines:
            self._stable = "stable"
        elif " The wavefunction has an internal instability.\n" in self._lines:
            self._stable = "internal instability"
        elif " The wavefunction has an RHF -> UHF instability.\n" in self._lines:
            self._stable = "RHF instability"
        else:
            self._stable = "untested"

    @property
    def homo(self) -> float:
        """
        Function to get the HOMO energy with reference to the vacuum level
        """
        return self.get_orbitals().query('occupied == True').query('energy == energy.max()')['energy'].iloc[0]
    
    @property
    def lumo(self) -> float:
        """
        Function to get the LUMO energy with reference to the vacuum level
        """
        return self.get_orbitals().query('occupied == False').query('energy == energy.min()')['energy'].iloc[0]

    @property
    def bandgap(self) -> float:
        """
        Function to get the band gap from the log file
        """
        return self.lumo - self.homo

    @property
    def n_alpha(self) -> int:
        return self._n_alpha
    
    @property
    def n_beta(self) -> int:
        return self._n_beta

    @property
    def time_stamp(self) -> dt.DateTime:
        return self._time_stamp.strftime("%Y-%m-%d %H:%M:%S")

    @property
    def scf_iterations(self) -> int:
        return self._scf_iterations
    
    @property
    def pop(self) -> bool:
        return self._pop
    
    @property
    def solvent(self) -> bool:
        return self._solvent

    @property
    def stable(self) -> str:
        return self._stable

    @property
    def restart(self) -> bool:
        return self._restart

    @property
    def esp(self) -> bool:
        return self._esp

    @property
    def heavyatomcount(self) -> int:
        return self._heavyatomcount

    @property
    def heavyatoms(self) -> list:
        return self._heavyatoms

    @property
    def atoms(self) -> list:
        return self._atoms
    
    @property
    def freq(self) -> bool:
        return self._freq

    @property
    def unrestricted(self) -> bool:
        return self._unrestricted

    @property
    def functional(self) -> str:
        if self._functional[0] == "U" or self._functional[0] == "R":
            return self._functional[1:]
        else:
            return self._functional

    @property
    def basis(self) -> str:
        return self._basis

    @property
    def energy(self) -> float:
        return self._energy

    @property
    def charge(self) -> int:
        return self._charge

    @property
    def raman(self) -> bool:
        return self._raman

    @property
    def opt(self) -> bool:
        return self._opt

    @property
    def multiplicity(self) -> int:
        return self._multiplicity

    @property
    def keywords(self) -> list:
        return self._keywords

    @property
    def log_file(self) -> str:
        return self._log_file

    @property
    def complete(self) -> bool:
        return self._complete

    @property
    def atomcount(self) -> int:
        return self._atomcount
    
    def _concatenate_log_files(self, log_file: list[str] | tuple[str] | pd.Series | str) -> tuple[list[str], bool, dt.DateTime]:
        """
        function to concatenate log files
        :return:
        """
        # If the file passed is just a string
        if type(log_file) == str:
            lines = [line for line in open(log_file, 'r')]
            restart = False
            time_stamp = self._get_time_stamp(log_file)

        # If a list or tuple of log files is passed
        elif (type(log_file) == list or type(log_file) == tuple or type(log_file) == pd.Series):
            log_file_dict = {}
            lines = []
            for l in log_file:
                time_stamp = self._get_time_stamp(l)
                log_file_dict[time_stamp] = l
            for key, value in sorted(log_file_dict.items()):
                lines = lines + [line for line in open(value, 'r')]
            restart = True
            time_stamp = min(log_file_dict)
        else:
            raise ValueError("The log file must be a path, a list of paths, a tuple of paths or a pd.Series of paths")
        
        return lines, restart, time_stamp
    
    @staticmethod
    def _get_time_stamp(log_file):
        """ 
        Function to read the contents of a log file, and from that get the line containing Leave Link and construct the time stamp from that line
        """
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        if any('Leave Link  ' in l for l in lines):
            time_line = [l for l in lines if 'Leave Link  ' in l][0]
            year = time_line.split()[8][:-1]
            month = time_line.split()[5]
            day = time_line.split()[6]
            time = time_line.split()[7]
            time_stamp = dt.DateTime(year + '-' + month + '-' + day + ' ' + time)
        elif any('Normal termination of' in l for l in lines):
            time_line = [l for l in lines[-20:] if 'Normal termination of' in l][0]
            year = time_line.split()[10][:-1]
            month = time_line.split()[7]
            day = time_line.split()[8]
            time = time_line.split()[9]
            time_stamp = dt.DateTime(year + '-' + month + '-' + day + ' ' + time)
        else:
            time_stamp = None

        return time_stamp

    def _get_keywords(self):
        """
        Function to extract the keywords from self._lines
        :return:
        """
        star_value = [i for i in self._lines if "****" in i][0]
        index = self._lines.index(star_value)
        temp_lines = self._lines[index+4:index+20]
        dash_value = [i for i in temp_lines if "--------" in i][0]
        index = temp_lines.index(dash_value)
        temp_lines = temp_lines[index+1:]
        index = temp_lines.index(dash_value)
        keyword_lines = temp_lines[0:index]
        keywords_str = ''
        for i in keyword_lines:
            keywords_str = keywords_str + i[1:-1]
        keywords = [i.lower() for i in keywords_str.split()][1:]

        if 'restart' in keywords or 'Restart' in keywords:
            raise ValueError("This log file is a restart file and should be parsed with the previous log file")

        return keywords

    def get_scf_convergence(self) -> pd.DataFrame:
        """
        Function to extract the SCF convergence data from the log file
        :return: pandas data frame of the SCF convergence data
        """
        if self._opt is False:
            raise ValueError("Your log file needs to be from an optimisation")

        scf_to_ignore = [i + 2 for i, line in enumerate(self._lines) if '>>>>>>>>>> Convergence criterion not met' in line]
        scf_lines = [line for i, line in enumerate(self._lines) if i not in scf_to_ignore and 'SCF Done' in line]

        data = (pd
                .DataFrame({
                    'iteration': [i for i in range(len(scf_lines))],
                    'energy': [float(s.split()[4]) * 2625.5 for s in scf_lines],
                    'cycles': [int(s.split()[7]) for s in scf_lines]})
                .assign(de = lambda x: x['energy'].diff())
                .assign(energy = lambda x: x['energy'] - x['energy'].iloc[-1])
                )
        
        data['de'].iloc[0] = data['de'].iloc[1]
        
        return data

    def get_thermo_chemistry(self) -> pd.DataFrame:
        """
        Function to extract the thermochemistry values from a log file
        :return: thermochemistry values
        """

        if self._freq is False:
            raise ValueError("Your log file needs to be from a vibrational analysis")

        zero_point_correction = float([line for line in self._lines if "Zero-point correction" in line][0].split()[2]) * 2625.5
        energy_thermal_cor = float([line for line in self._lines if "Thermal correction to Energy=" in line][0].split()[4]) * 2625.5
        enthalpy_thermal_cor = float([line for line in self._lines if "Thermal correction to Enthalpy=" in line][0].split()[4]) * 2625.5
        gibbs_thermal_cor = float([line for line in self._lines if "Thermal correction to Gibbs Free Energy=" in line][0].split()[6]) * 2625.5
        electronic_and_zp = float([line for line in self._lines if "Sum of electronic and zero-point Energies=" in line][0].split()[6]) * 2625.5
        elec_and_thermal_e = float([line for line in self._lines if "Sum of electronic and thermal Energies=" in line][0].split()[6]) * 2625.5
        elec_and_thermal_s = float([line for line in self._lines if "Sum of electronic and thermal Enthalpies=" in line][0].split()[6]) * 2625.5
        elec_and_thermal_g = float([line for line in self._lines if "Sum of electronic and thermal Free Energies=" in line][0].split()[7]) * 2625.5

        return pd.DataFrame({
            'zp_corr': [zero_point_correction],
            'e_corr': [energy_thermal_cor],
            's_corr': [enthalpy_thermal_cor],
            'g_corr': [gibbs_thermal_cor],
            'e_elec_zp': [electronic_and_zp],
            'e_elec_therm': [elec_and_thermal_e],
            's_elec_therm': [elec_and_thermal_s],
            'g_elec_therm': [elec_and_thermal_g]
        })

    def get_bonds_from_log(self):
        """
        function to extract the bond data from the gaussian log file. Use with caution - it doesn't always seem to get
        the bond data right
        :return:
        """
        if self._opt is False:
            start_line = len(self._lines) - (self._lines[::-1].index([k for k in self._lines if '!    Initial Parameters    !' in k][0]) - 4)
            end_line = len(self._lines) - self._lines[::-1].index([k for k in self._lines if 'Trust Radius=' in k][0]) + 2
        else:
            start_line = len(self._lines) - (self._lines[::-1].index([k for k in self._lines if '!   Optimized Parameters   !' in k][0]) - 4)
            end_line = len(self._lines) - (self._lines[::-1].index([k for k in self._lines if 'Largest change from initial coordinates is atom' in k][0]) + 2)

        bond_lines = [r for r in self._lines[start_line:end_line] if '! R' in r]

        data = (pd
                .DataFrame({
                    'atom_id_1': [int(r.split()[2][2:][:-1].split(",")[0]) for r in bond_lines],
                    'atom_id_2': [int(r.split()[2][2:][:-1].split(",")[1]) for r in bond_lines],
                    'length': [float(r.split()[3]) for r in bond_lines]
                })
                .assign(
                    element_1=lambda x: [self._atoms[i - 1] for i in x['atom_id_1']],
                    element_2=lambda x: [self._atoms[i - 1] for i in x['atom_id_2']]
                ))

        return data

    def get_spin_contamination(self) -> pd.DataFrame:
        """
        Function to get the spin contamination from a log file
        :return: pandas data frame of the spin contamination
        """
        contamination_lines = [s for s in self._lines if "S**2 before annihilation" in s]
        data = pd.DataFrame({
            'iteration': [i for i in range(len(contamination_lines))],
            'before_annihilation': [float(s.split()[3][:-1]) for s in contamination_lines],
            'after_annihilation': [float(s.split()[5]) for s in contamination_lines]
        })
        return data

    def get_bonds_from_coordinates(self, cutoff: float = 1.8, heavy_atoms: bool = False, scf_iteration: int = -1):
        """
        function to get bond data from the coordinates, using a cut-off distance
        :param cutoff: The cutoff for calculating the bond lengths
        :param heavy_atoms: just get the bonds involving heavy atoms
        :param pre_optimisation: get the coordinated before the optimisation has begun?
        :return:
        """
        coordinates = self.get_coordinates(heavy_atoms=heavy_atoms, scf_iteration=scf_iteration)

        cross = (coordinates
                 .merge(coordinates, how='cross', suffixes=('_1', '_2'))
                 .assign(dx=lambda x: x['x_2'] - x['x_1'])
                 .assign(dy=lambda x: x['y_2'] - x['y_1'])
                 .assign(dz=lambda x: x['z_2'] - x['z_1'])
                 .assign(length=lambda x: (x['dx'].pow(2) + x['dy'].pow(2) + x['dz'].pow(2)).pow(0.5))
                 .query('length < @cutoff')
                 .query('length > 0')
                 .filter(items=['atom_id_1', 'atom_id_2', 'length'])
                 .reset_index(drop=True)
                 .round(4)
                 )

        cross[['atom_id_1', 'atom_id_2']] = (np.sort(cross[['atom_id_1', 'atom_id_2']].to_numpy(), axis=1))

        data = (cross
                .groupby(['atom_id_1', 'atom_id_2'])
                .agg(length=('length', 'first'))
                .reset_index()
                .assign(element_1=lambda x: [self._atoms[i - 1] for i in x['atom_id_1']])
                .assign(element_2=lambda x: [self._atoms[i - 1] for i in x['atom_id_2']])
                )

        return data
    
    def get_coordinates_through_scf(self, heavy_atoms: bool = False) -> pd.DataFrame:
        """
        function to get the coordinates through the SCF iterations
        :param heavy_atoms: just get the heavy atoms?
        :return:
        """
        if self.opt is False:
            raise ValueError("This log file needs to be from an optimisation")

        data = pd.DataFrame()
        for i in range(self._scf_iterations):
            new_data = self.get_coordinates(heavy_atoms=heavy_atoms, scf_iteration=i).assign(iteration=i)
            data = pd.concat([data, new_data])

        return data

    def get_coordinates(self, heavy_atoms: bool = False, scf_iteration: int = -1) -> pd.DataFrame:
        """
        function to get the coordinates from the log file
        :param heavy_atoms: return just the heavy atoms?
        :param scf_interation: get the coordinates at this scf iteration. If 0, then before optimisation has begun
        :return:
        """
        indices = [i for i, line in enumerate(self._lines) if 'Standard orientation:' in line]
        start_line = indices[scf_iteration] + 5
        end_line = start_line + self._atomcount

        data = (pd.DataFrame({
            'atom_id': [i for i in range(1, self._atomcount + 1)],
            'element': self._atoms,
            'x': [float(a.split()[3]) for a in self._lines[start_line:end_line]],
            'y': [float(a.split()[4]) for a in self._lines[start_line:end_line]],
            'z': [float(a.split()[5]) for a in self._lines[start_line:end_line]]
        }))

        if heavy_atoms is False:
            return data
        elif heavy_atoms is True:
            return data.query("element != 'H'")

    def get_mulliken_charges(self, heavy_atoms: bool = False, with_coordinates: bool = False, **kwargs) -> pd.DataFrame:
        """
        method to return the mulliken charges from the log file
        :param heavy_atoms: whether to give the heavy atoms or all the atoms
        :param with_coordinates: whether to also output coordinates
        :return:
        """
        start_line = len(self._lines) - self._lines[::-1].index([k for k in self._lines if 'Mulliken charges' in k][0]) + 1
        end_line = start_line + self._atomcount
        if heavy_atoms is False:
            data = pd.DataFrame({
                "atom_id": [int(a.split()[0]) for a in self._lines[start_line:end_line]],
                "element": [a.split()[1] for a in self._lines[start_line:end_line]],
                "partial_charge": [float(a.split()[2]) for a in self._lines[start_line:end_line]]
            })
        else:
            start_line = start_line + self._atomcount + 3
            end_line = start_line + self._heavyatomcount
            data = pd.DataFrame({
                "atom_id": [int(a.split()[0]) for a in self._lines[start_line:end_line]],
                "element": [a.split()[1] for a in self._lines[start_line:end_line]],
                "partial_charge": [float(a.split()[2]) for a in self._lines[start_line:end_line]]
            })

        if with_coordinates is True:
            data = (data.assign(
                x=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['x'].to_list(),
                y=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['y'].to_list(),
                z=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['z'].to_list()
            ))

        return data
    
    def get_mulliken_spin_densities(self, heavy_atoms: bool = False, with_coordinates: bool = False, **kwargs) -> pd.DataFrame:
        """
        method to return the mulliken spin densities from the log file
        :param heavy_atoms: whether to give the heavy atoms or all the atoms
        :param with_coordinates: whether to also output coordinates
        :return: pandas dataframe with the spin densities
        """
        # check whether spin density information is in the log file
        if any('Mulliken charges and spin densities:' in s for s in self._lines) is False:
            spins = False
        else:
            spins = True

        start_line = len(self._lines) - self._lines[::-1].index([k for k in self._lines if 'Mulliken charges and spin densities:' in k][0]) + 1
        end_line = start_line + self._atomcount

        if heavy_atoms is False:
            data = pd.DataFrame({
                "atom_id": [int(a.split()[0]) for a in self._lines[start_line:end_line]],
                "element": [a.split()[1] for a in self._lines[start_line:end_line]]
            })

            if spins == True:
                data = data.assign(spin_density=[float(a.split()[3]) for a in self._lines[start_line:end_line]])
            else:
                data = data.assign(spin_density = 0)

        else:
            start_line = start_line + self._atomcount + 3
            end_line = start_line + self._heavyatomcount
            data = pd.DataFrame({
                "atom_id": [int(a.split()[0]) for a in self._lines[start_line:end_line]],
                "element": [a.split()[1] for a in self._lines[start_line:end_line]]
            })

            if spins == True:
                data = data.assign(spin_density=[float(a.split()[3]) for a in self._lines[start_line:end_line]])
            else:
                data = data.assign(spin_density = 0)

        if with_coordinates is True:
            data = (data.assign(
                x=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['x'].to_list(),
                y=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['y'].to_list(),
                z=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['z'].to_list()
            ))

        return data

    def get_esp_charges(self, heavy_atoms: bool = False, with_coordinates: bool = False, **kwargs) -> pd.DataFrame:
        """
        method to return the mulliken charges from the log file
        :param heavy_atoms: whether to give the heavy atoms or all the atoms
        :param with_coordinates: whether to also output the x,y,z coordinates
        :return:
        """
        if self._esp is False:
            raise ValueError("This gaussian log file doesnt have ESP data in it!")

        start_line = len(self._lines) - self._lines[::-1].index([k for k in self._lines if 'ESP charges' in k][0]) + 1
        end_line = start_line + self._atomcount

        if heavy_atoms is False:
            data = pd.DataFrame({
                "atom_id": [int(a.split()[0]) for a in self._lines[start_line:end_line]],
                "element": [a.split()[1] for a in self._lines[start_line:end_line]],
                "partial_charge": [float(a.split()[2]) for a in self._lines[start_line:end_line]]
            })
        else:
            start_line = start_line + self._atomcount + 3
            end_line = start_line + self._heavyatomcount
            data = pd.DataFrame({
                "atom_id": [int(a.split()[0]) for a in self._lines[start_line:end_line]],
                "element": [a.split()[1] for a in self._lines[start_line:end_line]],
                "partial_charge": [float(a.split()[2]) for a in self._lines[start_line:end_line] if 'H' not in a]
            })

        if with_coordinates is True:
            data = (data.assign(
                x=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['x'].to_list(),
                y=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['y'].to_list(),
                z=self.get_coordinates(heavy_atoms=heavy_atoms, **kwargs)['z'].to_list()
            ))

        return data

    def get_raman_frequencies(self, frac_filter: float = 1) -> pd.DataFrame:
        """
        method to get the raman frequencies from the log file
        :param frac_filter: take this fraction of peaks with the highest intensities
        :return:
        """
        if self._raman is False:
            raise ValueError("This gaussian log file doesnt have raman data in it!")

        if frac_filter < 0 or frac_filter > 1:
            raise ValueError("frac_filter must be between 0 and 1!")

        frequencies = [line.split("--")[1].split() for line in self._lines if "Frequencies --" in line]
        frequencies = [item for sublist in frequencies for item in sublist]
        activities = [line.split("--")[1].split() for line in self._lines if "Raman Activ --" in line]
        activities = [item for sublist in activities for item in sublist]

        data = (pd
                .DataFrame({"frequencies": frequencies, "raman_activity": activities})
                .query('raman_activity != "************"')
                .query('frequencies != "************"')
                .assign(frequencies=lambda x: x['frequencies'].astype('float'))
                .assign(raman_activity=lambda x: x['raman_activity'].astype('float'))
                )

        cutoff = data['raman_activity'].max() * (1 - frac_filter)

        data = (data
                .query('raman_activity > @cutoff')
                .assign(cutoff=cutoff)
                )

        return data

    def get_raman_spectra(self, width: float = 20, wn_min: int = 500, wn_max: int = 2500, wn_step: float = 1, **kwargs):
        """
        method to get a theoretical spectrum from the gaussian log file
        :param width: the width of the lorentzian peaks
        :param wn_min: the minimum wave number
        :param wn_max: the maximum wave number
        :param wn_step: the number of intervals in the spectrum
        :return:
        """
        peaks = self.get_raman_frequencies(**kwargs)
        wn = [w for w in np.arange(wn_min, wn_max, wn_step)]
        intensity = [0] * len(wn)

        for index, row in peaks.iterrows():
            peak = lorentzian(wn, row['frequencies'], width, row['raman_activity'])
            intensity = [sum(intensity) for intensity in zip(intensity, peak)]

        return pd.DataFrame({'wavenumber': wn, 'intensity': intensity})
    
    def get_optimisation_trajectory(self, filename: str, path: str = '.', fit_t0 = False):
        """
        Function to get the optimisation trajectory from the log file
        """
        coordinates = self.get_coordinates_through_scf()
        PdbParser.pandas_to_pdb_trajectory(coordinates, time_col='iteration', filename=filename, path=path, fit_t0=fit_t0)
        return None
    
    def get_orbitals(self):
        """
        Function to get the orbital information from the log file
        """
        if any('Population analysis using the SCF Density' in s for s in self._lines) is False:
            raise ValueError("This log file doesnt have molecular orbital data in it")
        
        start_line = [i for i, line in enumerate(self._lines) if 'Population analysis using the SCF Density' in line][-1] + 4
        end_line = [i for i, line in enumerate(self._lines) if 'Condensed to atoms (all electrons)' in line][-1]
        chunk = self._lines[start_line:end_line]

        alpha_len = len([l for l in chunk if 'Alpha  occ. eigenvalues' in l])
        beta_len = len([l for l in chunk if 'Beta  occ. eigenvalues' in l])

        # Get the alpha molecular orbitals
        if alpha_len > 1:
            alpha_occ = [l for l in chunk if 'Alpha  occ. eigenvalues' in l]
            alpha_virt = [l for l in chunk if 'Alpha virt. eigenvalues' in l]
            alpha_occ = [s.split()[4:] for s in alpha_occ]
            alpha_virt = [s.split()[4:] for s in alpha_virt]
            alpha_occ = [item for sublist in alpha_occ for item in sublist]
            alpha_virt = [item for sublist in alpha_virt for item in sublist]
            alpha_occ = [float(a)*2625.5 for a in alpha_occ]
            alpha_virt = [float(a)*2625.5 for a in alpha_virt]
            alpha_occ_df = pd.DataFrame({'energy': alpha_occ, 'occupied': True})
            alpha_virt_df = pd.DataFrame({'energy': alpha_virt, 'occupied': False})

            alpha_df = (pd
                        .concat([alpha_occ_df, alpha_virt_df])
                        .assign(orbital_number = lambda x: [i for i in range(1, len(x) + 1)])
                        )

        # Get the beta molecular orbitals
        if beta_len > 1:
            beta_occ = [l for l in chunk if 'Beta  occ. eigenvalues' in l]
            beta_virt = [l for l in chunk if 'Beta virt. eigenvalues' in l]
            beta_occ = [s.split()[4:] for s in beta_occ]
            beta_virt = [s.split()[4:] for s in beta_virt]
            beta_occ = [item for sublist in beta_occ for item in sublist]
            beta_virt = [item for sublist in beta_virt for item in sublist]
            beta_occ = [float(a)*2625.5 for a in beta_occ]
            beta_virt = [float(a)*2625.5 for a in beta_virt]
            beta_occ_df = pd.DataFrame({'energy': beta_occ, 'occupied': True})
            beta_virt_df = pd.DataFrame({'energy': beta_virt, 'occupied': False})

            beta_df = (pd
                       .concat([beta_occ_df, beta_virt_df])
                       .assign(orbital_number = lambda x: [i for i in range(1, len(x) + 1)])
                      )

        if alpha_len > 1 and beta_len > 1:
            data = pd.concat([alpha_df.assign(electron = 'alpha'), beta_df.assign(electron = 'beta')])
        elif alpha_len > 1 and beta_len <= 1:
            data = alpha_df.assign(electron = 'paired')
        else:
            raise ValueError("This log file doesnt have molecular orbital data in it!")

        return data
    
    def get_dos_plot(self, **kwargs):
        """ Function to return a density of states plot """
        data = self.get_orbitals()
        bins = len(data) // 2
        figure = px.histogram(data, x='energy', marginal='rug', color='occupied', nbins=bins, 
                              pattern_shape='electron', labels={'energy': 'E [kJ/mol]'}, **kwargs)
        return figure
