from Materials_Data_Analytics.experiment_modelling.core import ElectrochemicalMeasurement
from Materials_Data_Analytics.laws_and_constants import NA, KB, F, E
from Materials_Data_Analytics.materials.electrolytes import Electrolyte
import pandas as pd
import plotly.express as px
import io

class ImpedanceSpectra(ElectrochemicalMeasurement):
    """
    A class to represent impedance spectroscopy data.

    Main contributors:
    L. Julian Mele

    Contributors:
    Nicholas Siemons
    """
    
    def __init__(self, metadata: dict = None) -> None:

        super().__init__(metadata=metadata)   # metadata is a dictionary when you make the object (see experiment() ) you type metadata={'independent variable': XXX}
        self._data = pd.DataFrame(dtype=float)
        # check columns names
        self._e_vs_ref = None
        self._va_mv = None

    def __str__(self):
        # Insert here all the metadata from the mps file if available
        return f"ImpedanceSpectra(electrolyte={self.electrolyte}, metadata={self.metadata}, frequency={self._frequency}, impedance_data={self._impedance_data})"

    @property
    def frequency(self):
        return self._frequency
    
    @property
    def impedance_data(self):
        return self._impedance_data
    
    @property
    def admittance_data(self):
        return self._admittance_data
    
    @property
    def phase_data(self):
        return self._phase_data    

    @property
    def e_vs_ref(self):
        return self._e_vs_ref

    @property
    def va_mv(self):
        return self._va_mv

    @classmethod
    def from_biologic(cls, path: str = None, data: pd.DataFrame = None, **kwargs):
        """
        Create an ImpedanceSpectra object from a Biologic EC-Lab ASCII file.
        
        The file is expected to have a header section (with a line such as
        "Nb header lines : 66") followed by a data block whose first non-empty
        line gives the column names. The rest of the file is numerical data.
        """
        # If a DataFrame is provided directly, use it
        if path is None and data is not None:
            metadata = kwargs.pop('metadata', {})
            obj = cls(metadata=metadata, **kwargs)
            obj._data = data
            return obj
        elif path is None:
            raise ValueError("Either a file path or a DataFrame must be provided.")

        # Read the entire file
        with open(path, 'r', encoding='cp1252') as f:
            lines = f.readlines()

        # Look for "Nb header lines" in the file to determine the header length.
        header_line_count = None
        for line in lines:
            if "Nb header lines" in line:
                try:
                    header_line_count = int(line.split(":", 1)[1].strip())
                except Exception:
                    header_line_count = None
                break
        # If not found, default to 66
        if header_line_count is None:
            header_line_count = 66

        # Parse header lines and build metadata dictionary.
        header_lines = lines[:header_line_count]
        metadata = {}
        for line in header_lines:
            if ":" in line:
                key, value = line.split(":", 1)
                metadata[key.strip()] = value.strip()

        # Get data lines (skip any empty lines after the header)
        data_lines = [line for line in lines[header_line_count:] if line.strip()]
        if not data_lines:
            raise ValueError("No data found after header in file.")

        # The first nonempty line contains the column headers.
        col_header_line = header_lines[header_line_count-1]
        # Replace any occurrence of two or more spaces with a tab.
        while "  " in col_header_line:
            col_header_line = col_header_line.replace("  ", "\t")
        col_names = col_header_line.strip().split("\t")

        # The remaining lines are the numerical data.
        data_string = "\n".join(data_lines[1:])
        df = pd.read_csv(io.StringIO(data_string), delim_whitespace=True, header=None, names=col_names)

        # Standardize column names by renaming them.
        rename_map = {
            'freq/Hz': 'frequency',
            'Re(Z)/Ohm': 'real_impedance',
            '-Im(Z)/Ohm': 'minus_imag_impedance',
            '|Z|/Ohm': 'abs_impedance',
            'Phase(Z)/deg': 'phase',
            'cycle number': 'cycle_number',
            'THD Ewe/%': 'THD_Ewe',
            'NSD Ewe/%': 'NSD_Ewe',
            'NSR Ewe/%': 'NSR_Ewe',
            'THD I/%': 'THD_I',
            'NSD I/%': 'NSD_I',
            'NSR I/%': 'NSR_I',
            'time/s': 'time',
            '<Ewe>/V': 'Ewe',
            '<I>/mA': 'current',
            '|Ewe|/V': 'abs_Ewe',
            '|I|/A': 'abs_current',
            'Ns': 'N_s'
        }
        df.rename(columns=rename_map, inplace=True)

        # Create an instance and store the DataFrame in the _data attribute.
        obj = cls(metadata=metadata, **kwargs)
        obj._data = df
        return obj


    def get_nyquist_plot(self, **kwargs) -> None:
        """Plot the Nyquist plot."""
        fig = px.scatter(self.impedance_data, x='z_real', y='z_imag', **kwargs)
        fig.show()
        return fig
        
    def plot_bode(self) -> None:
        """Plot the Bode plot."""
        fig, ax1 = px.subplots()
        figure1 = px.line(self.data, x='magnitude', y='|Z|/Ohm', title='Bode Plot', labels={'|Z|/Ohm': 'Magnitude (Ohm)', 'frequency': 'Frequency (Hz)'})

        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Magnitude (Ohm)', color='tab:blue')
        ax1.plot(self.frequency, self.impedance_data['|Z|/Ohm'], 'b-')
        ax1.tick_params(axis='y', labelcolor='tab:blue')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Phase (degrees)', color='tab:red')
        ax2.plot(self.frequency, self.impedance_data['Phase(Z)/deg'], 'r-')
        ax2.tick_params(axis='y', labelcolor='tab:red')

        px.title('Bode Plot')
        px.show()

    def fit_model(self, model) -> None:
        """Fit an equivalent circuit model to the data."""
        # Implementation for fitting the model
        pass


class PEIS(ImpedanceSpectra):
    def __init__(self, electrolyte: Electrolyte = None, metadata: dict = None) -> None:
        super().__init__(electrolyte=electrolyte, metadata=metadata)
        # Additional attributes specific to PEIS

    # Additional methods specific to PEIS

class SPEIS(PEIS):
    def __init__(self, electrolyte: Electrolyte = None, metadata: dict = None) -> None:
        super().__init__(electrolyte=electrolyte, metadata=metadata)
        # Additional attributes specific to SPEIS

    # Additional methods specific to SPEIS