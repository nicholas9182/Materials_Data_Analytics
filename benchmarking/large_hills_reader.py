import pandas as pd
import numpy as np
import plumed as pl
from timer import Timer
from Materials_Data_Analytics.metadynamics.free_energy import FreeEnergySpace

if __name__ == "__main__":

    t = Timer()
    file = "/Users/nicholassiemons/Dropbox/OBT/0085/010/HILLS"

    t.start()
    r1 = pl.read_as_pandas(file)
    t.stop()

    t.start()
    col_names = open(file).readline().strip().split(" ")[2:]
    r2 = pd.read_table(file, delim_whitespace=True, comment="#", header=None, names=col_names, dtype=np.float64)
    t.stop()

