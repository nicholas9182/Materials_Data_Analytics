import plumed as pl
import pandas as pd
from timer import Timer
from analytics.metadynamics.free_energy import FreeEnergySpace


def benchmark_plotting_hills_v1(time_res: int) -> pd.DataFrame:
    """
    old way of formatting the _hills that is very slow using apply and assign
    :return:
    """
    data = pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0091/002/HILLS")
    landscape = FreeEnergySpace(data)
    hills = landscape._hills
    cvs = landscape.cvs
    long_hills = (hills
                  .melt(value_vars=cvs+['height'], id_vars=['time', 'walker'])
                  .assign(time=lambda x: x['time'].round(time_res))
                  .groupby(['time', 'walker', 'variable'], group_keys=False)
                  .apply(lambda x: x.assign(value=lambda y: y['value'].mean()))
                  .drop_duplicates()
                  )

    return long_hills


def benchmark_plotting_hills_v2(time_res: int) -> pd.DataFrame:
    """
    option where we just take a random row in each group by dropping duplicates across walker, time and variable. Not ideal but seems much faster
    :return:
    """
    data = pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0091/002/HILLS")
    landscape = FreeEnergySpace(data)
    hills = landscape._hills
    cvs = landscape.cvs
    long_hills = (hills
                  .melt(value_vars=cvs+['height'], id_vars=['time', 'walker'])
                  .assign(time=lambda x: x['time'].round(time_res))
                  .drop_duplicates(subset=['time', 'walker', 'variable'])
                  )

    return long_hills


def benchmark_plotting_hills_v3(time_res: int) -> pd.DataFrame:
    """
    try using the groupby and .mean()
    :return:
    """
    data = pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0091/002/HILLS")
    landscape = FreeEnergySpace(data)
    hills = landscape._hills
    cvs = landscape.cvs
    long_hills = (hills
                  .melt(value_vars=cvs+['height'], id_vars=['time', 'walker'])
                  .assign(time=lambda x: x['time'].round(time_res))
                  .groupby(['time', 'walker', 'variable'], group_keys=False)
                  .mean()
                  .reset_index()
                  )

    return long_hills


if __name__ == "__main__":

    t = Timer()

    t.start()
    r1 = benchmark_plotting_hills_v1(1)
    t.stop()

    t.start()
    r2 = benchmark_plotting_hills_v2(1)
    t.stop()

    t.start()
    r3 = benchmark_plotting_hills_v3(1)
    t.stop()
