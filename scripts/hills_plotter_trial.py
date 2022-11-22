import pandas as pd
import plumed as pl
import matplotlib.pyplot as plt


def plot_hills(df):
    plt.plot(df['time'], df['height'])
    plt.show()


hills = (pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0085/004/HILLS_test")
         .groupby('time', group_keys=True)
         .apply(lambda x: x.assign(walker=lambda y: range(0, y.shape[0])))
         .query('walker == 0')
         .pipe(plot_hills)
         )
