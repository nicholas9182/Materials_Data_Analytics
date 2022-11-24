import plumed as pl
import pandas as pd
import numpy as np
from timer import Timer
from analytics.metadynamics.free_energy import FreeEnergyLandscape


def benchmark_plotting_hills_v1() -> pd.DataFrame:
    """
    old way of formatting the hills that is very slow
    :return:
    """
    data = pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0091/002/HILLS")
    landscape = FreeEnergyLandscape(data)
    hills = landscape.hills
    cvs = landscape.cvs
    long_hills = (hills
                  .melt(value_vars=cvs+['height'], id_vars=['time', 'walker'])
                  .assign(time=lambda x: x['time'].round(1))
                  .groupby(['time', 'walker', 'variable'], group_keys=False)
                  .apply(lambda x: x.assign(value=lambda y: y['value'].mean()))
                  .drop_duplicates()
                  )

    return long_hills


def benchmark_plotting_hills_v2() -> pd.DataFrame:
    """
    old way of formatting the hills that is very slow
    :return:
    """
    data = pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0091/002/HILLS")
    landscape = FreeEnergyLandscape(data)
    hills = landscape.hills
    cvs = landscape.cvs
    long_hills = (hills
                  .melt(value_vars=cvs+['height'], id_vars=['time', 'walker'])
                  .assign(time=lambda x: x['time'].round(1))
                  .groupby(['time', 'walker', 'variable'], group_keys=False)
                  .apply(lambda x: x.assign(value=lambda y: y['value'].mean()))
                  .drop_duplicates()
                  )

    return long_hills


if __name__ == "__main__":

    t = Timer()

    t.start()
    r1 = benchmark_plotting_hills_v1()
    t.stop()

    t.start()
    r2 = benchmark_plotting_hills_v2()
    t.stop()

    print(r1.equals(r2))
