import plumed as pl
from analytics.metadynamics.free_energy import FreeEnergySpace

hills2 = pl.read_as_pandas("../test_trajectories/ndi_na_binding/HILLS")
landscape = FreeEnergySpace(hills2)

av_hills = (landscape
            ._hills
            .assign(time=lambda x: x['time'].round(1))
            .groupby(['time'])
            .mean()
            .reset_index()
            )

print(av_hills[['time', 'height']])

max_hills = (landscape
             ._hills
             .groupby(['time'])
             .max()
             .reset_index()
             )

print(max_hills[['time', 'height']])
