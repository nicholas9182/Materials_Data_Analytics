import pandas as pd
from analytics.laws_and_constants import lorentzian
pd.set_option('mode.chained_assignment', None)


class GaussianParser:
    """
    Class to parse information from a gaussian log file
    """
    def __init__(self, log_file: str):

        self._log_file = log_file
        self._lines = [line for line in open(log_file, 'r')]
        self._keywords = [k for k in self._lines if '#p ' in k][0].split()[1:]
        self._charge = int([c for c in self._lines if 'Charge =' in c][0].split()[2])
        self._multiplicity = int([m for m in self._lines if 'Charge =' in m][0].split()[5])
        self._raman = True if len([r for r in self._keywords if 'raman' in r]) > 0 else False
        self._opt = True if len([r for r in self._keywords if 'opt' in r]) > 0 else False
        self._complete = True if len([r for r in self._lines if 'Normal termination of ' in r]) > 0 else False
        start_line = self._lines.index([k for k in self._lines if 'Mulliken charges' in k][0])
        end_line = self._lines.index([k for k in self._lines if 'Sum of Mulliken charges' in k][0])
        self._atomcount = end_line - start_line - 2
        self._energy = float([e for e in self._lines if 'SCF Done' in e][-1].split()[4])*2625.5

    @property
    def energy(self):
        return self._energy

    @property
    def charge(self):
        return self._charge

    @property
    def raman(self):
        return self._raman

    @property
    def opt(self):
        return self._opt

    @property
    def multiplicity(self):
        return self._multiplicity

    @property
    def keywords(self):
        return self._keywords

    @property
    def log_file(self):
        return self._log_file

    @property
    def complete(self):
        return self._complete

    @property
    def atomcount(self):
        return self._atomcount

    def get_raman_frequencies(self) -> pd.DataFrame:
        """
        method to get the raman frequencies from the log file
        :return:
        """
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

        return data

    def get_raman_spectra(self, width: float = 20, wn_min: int = 500, wn_max: int = 2500, wn_step: int = 1):
        """
        method to get a theoretical spectram from the gaussian log file
        :param width: the width of the lorentzian peaks
        :param wn_min: the minimum wave number
        :param wn_max: the maximum wave number
        :param wn_step: the number of intervals in the spectrum
        :return:
        """
        peaks = self.get_raman_frequencies()
        wn = [w for w in range(wn_min, wn_max, wn_step)]
        intensity = [0]*len(wn)

        for index, row in peaks.iterrows():
            peak = lorentzian(wn, row['frequencies'], width, row['raman_activity'])
            intensity = [sum(intensity) for intensity in zip(intensity, peak)]

        return pd.DataFrame({'wavenumber': wn, 'intensity': intensity})
