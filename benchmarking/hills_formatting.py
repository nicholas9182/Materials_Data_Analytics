import plumed as pl
import pandas as pd
from timer import Timer


def benchmark_formatting_hills_v1() -> pd.DataFrame:
    """
    old way of formatting the hills that is very slow
    :return:
    """
    data = pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0091/002/HILLS")
    hills = (pd.DataFrame(data)
             .loc[:, ~data.columns.str.startswith('sigma')]
             .drop(columns=['biasf'])
             .assign(time=lambda x: x['time']/1000)
             .groupby('time', group_keys=False).apply(lambda x: x.assign(walker=lambda y: range(0, y.shape[0])))
             )

    hills = hills[hills['time'] < max(hills['time'])]

    return hills


def benchmark_formatting_hills_v2() -> pd.DataFrame:
    """
    new way which calculates some other properties. Seems much faster.
    :return:
    """
    data = pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0091/002/HILLS")

    hills = (pd.DataFrame(data)
             .loc[:, ~data.columns.str.startswith('sigma')]
             .drop(columns=['biasf'])
             .assign(time=lambda x: x['time']/1000)
             )

    hills = hills[hills['time'] < max(hills['time'])]

    walker_num = hills[hills['time'] == min(hills['time'])].shape[0]
    t_steps = hills[['time']].drop_duplicates().shape[0]
    hills['walker'] = [i for i in range(0, walker_num)]*t_steps

    return hills


if __name__ == "__main__":

    t = Timer()

    t.start()
    r1 = benchmark_formatting_hills_v1()
    t.stop()

    t.start()
    r2 = benchmark_formatting_hills_v2()
    t.stop()

    print(r1.equals(r2))
