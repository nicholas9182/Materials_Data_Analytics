import pandas as pd
import numpy as np
from Materials_Data_Analytics.laws_and_constants import lorentzian

pd.set_option('mode.chained_assignment', None)


class GaussianParser:
    """
    Class to parse information from a gaussian log file
    """

    def __init__(self, log_file: str | list[str]):

        self._log_file = log_file

        if type(log_file) == str or (type(log_file) == list and len(log_file) == 1) or (type(log_file) == tuple and len(log_file) == 1):
            if type(log_file) == str:
                self._lines = [line for line in open(log_file, 'r')]
            elif len(log_file) == 1:
                self._lines = [line for line in open(log_file[0], 'r')]
            self._restart = False
        elif type(log_file) == list or type(log_file) == tuple:
            lines = []
            for file in log_file:
                lines = lines + [line for line in open(file, 'r')]
            self._lines = lines
            self._restart = True
        else:
            raise ValueError("The log file must be a path, a list of paths or a tuple of paths")
        
        self._keywords = self._get_keywords()
        self._raman = True if len([i for i in self.keywords if 'raman' in i]) > 0 else False
        self._opt = True if 'opt' in self._keywords else False
        self._complete = True if len([r for r in self._lines if 'Normal termination of ' in r]) > 0 else False
        self._esp = True if len([r for r in self._lines if 'ESP charges:' in r]) > 0 else False
        self._functional = [k for k in self._keywords if "/" in k][0].split("/")[0].upper()
        self._basis = [k for k in self._keywords if "/" in k][0].split("/")[1]

        if len([c for c in self._lines if 'Charge =' in c]) > 0:
            self._charge = int([c for c in self._lines if 'Charge =' in c][0][9:].split()[0])
            self._multiplicity = int([m for m in self._lines if 'Charge =' in m][0][27:])
        else:
            self._charge = None
            self._multiplicity = None

        if self._complete is True:
            self._energy = float([e for e in self._lines if 'SCF Done' in e][-1].split()[4]) * 2625.5
            self._unrestricted = True if [e for e in self._lines if 'SCF Done' in e][-1].split()[2][2] == "U" else False
            self._mull_start = self._lines.index([k for k in self._lines if 'Mulliken charges' in k][0]) + 2
            self._mull_end = self._lines.index([k for k in self._lines if 'Sum of Mulliken charges' in k][0])
            self._atomcount = self._mull_end - self._mull_start
            self._atoms = [a.split()[1] for a in self._lines[self._mull_start:self._mull_end]]
            self._heavyatoms = [a.split()[1] for a in self._lines[self._mull_start:self._mull_end] if 'H' not in a]
            self._heavyatomcount = len(self._heavyatoms)
        else:
            self._energy = None
            self._unrestricted = None
            self._mull_start = None
            self._mull_end = None
            self._atomcount = None
            self._atoms = None
            self._heavyatoms = None
            self._heavyatomcount = None

        if " The wavefunction is stable under the perturbations considered.\n" in self._lines:
            self._stable = "stable"
        elif " The wavefunction has an internal instability.\n" in self._lines:
            self._stable = "internal instability"
        elif " The wavefunction has an RHF -> UHF instability.\n" in self._lines:
            self._stable = "RHF instability"
        else:
            self._stable = "untested"

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


    @property
    def stable(self):
        return self._stable

    @property
    def restart(self):
        return self._restart

    @property
    def esp(self):
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

    def get_thermo_chemistry(self) -> pd.DataFrame:
        """
        Function to extract the thermochemistry values from a log file
        :return: thermochemistry values
        """
        if self._raman is False:
            raise ValueError("Your log file needs to be from a vibrational analysis")

        zero_point_correction = float([line for line in self._lines
                                       if "Zero-point correction" in line][0].split()[2]) * 2625.5
        energy_thermal_cor = float([line for line in self._lines
                                    if "Thermal correction to Energy=" in line][0].split()[4]) * 2625.5
        enthalpy_thermal_cor = float([line for line in self._lines
                                      if "Thermal correction to Enthalpy=" in line][0].split()[4]) * 2625.5
        gibbs_thermal_cor = float([line for line in self._lines
                                   if "Thermal correction to Gibbs Free Energy=" in line][0].split()[6]) * 2625.5
        electronic_and_zp = float([line for line in self._lines
                                   if "Sum of electronic and zero-point Energies=" in line][0].split()[6]) * 2625.5
        elec_and_thermal_e = float([line for line in self._lines
                                    if "Sum of electronic and thermal Energies=" in line][0].split()[6]) * 2625.5
        elec_and_thermal_s = float([line for line in self._lines
                                    if "Sum of electronic and thermal Enthalpies=" in line][0].split()[6]) * 2625.5
        elec_and_thermal_g = float([line for line in self._lines
                                    if "Sum of electronic and thermal Free Energies=" in line][0].split()[7]) * 2625.5

        return pd.DataFrame({
            'zp_corr': [zero_point_correction],
            'e_corr': [energy_thermal_cor],
            's_corr': [enthalpy_thermal_cor],
            'g_corr': [gibbs_thermal_cor],
            'e_elec_zp': [electronic_and_zp],
            'e_elec_therm': [elec_and_thermal_e],
            's_elec_therm': [elec_and_thermal_s],
            'g_elec_therm': [elec_and_thermal_g]
        }).round(5)

    def get_bonds_from_log(self):
        """
        function to extract the bond data from the gaussian log file. Use with caution - it doesn't always seem to get
        the bond data right
        :return:
        """
        if self._opt is False:
            start_line = len(self._lines) - \
                         (self._lines[::-1].index(
                             [k for k in self._lines if '!    Initial Parameters    !' in k][0]) - 4)
            end_line = len(self._lines) - self._lines[::-1].index(
                [k for k in self._lines if 'Trust Radius=' in k][0]) + 2
        else:
            start_line = len(self._lines) - \
                         (self._lines[::-1].index(
                             [k for k in self._lines if '!   Optimized Parameters   !' in k][0]) - 4)

            end_line = len(self._lines) - \
                       (self._lines[::-1].index(
                           [k for k in self._lines if 'Largest change from initial coordinates is atom' in k][0]) + 2)

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

    def get_bonds_from_coordinates(self, cutoff: float = 1.8, heavy_atoms: bool = False, pre_optimisation: bool = False):
        """
        function to get bond data from the coordinates, using a cut-off distance
        :param cutoff: The cutoff for calculating the bond lengths
        :param heavy_atoms: just get the bonds involving heavy atoms
        :param pre_optimisation: get the coordinated before the optimisation has begun?
        :return:
        """
        coordinates = self.get_coordinates(heavy_atoms=heavy_atoms, pre_optimisation=pre_optimisation)

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

    def get_coordinates(self, heavy_atoms: bool = False, pre_optimisation: bool = False) -> pd.DataFrame:
        """
        function to get the coordinates from the log file
        :param heavy_atoms: return just the heavy atoms?
        :param pre_optimisation: get the coordinated before the optimisation has begun?
        :return:
        """
        if pre_optimisation is False:
            start_line = len(self._lines) - (self
                                             ._lines[::-1]
                                             .index([k for k in self._lines if 'Standard orientation:' in k][0]) - 4
                                             )
        else:
            start_line = self._lines.index([k for k in self._lines if 'Standard orientation:' in k][0]) + 5

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
        start_line = len(self._lines) - self._lines[::-1].index(
            [k for k in self._lines if 'Mulliken charges' in k][0]) + 1
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

    def get_raman_frequencies(self, frac_filter: float = 0.99) -> pd.DataFrame:
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
