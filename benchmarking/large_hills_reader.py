import plumed as pl
from timer import Timer
from analytics.metadynamics.free_energy import FreeEnergySpace

if __name__ == "__main__":

    t = Timer()

    t.start()
    r1 = pl.read_as_pandas("/Users/nicholassiemons/Dropbox/OBT/0085/010/HILLS")
    t.stop()

    t.start()
    landscape = FreeEnergySpace(r1)
    t.stop()

    t.start()
    figures = landscape.get_hills_figures(time_resolution=1, height_power=0.5)
    t.stop()
