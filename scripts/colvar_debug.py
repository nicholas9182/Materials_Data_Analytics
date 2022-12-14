import os
from analytics.metadynamics.free_energy import FreeEnergySpace, FreeEnergyLine, MetaTrajectory, FreeEnergySurface

space = FreeEnergySpace("/Users/nicholassiemons/Dropbox/OBT/0091/001/HILLS")

files_d1 = ["/Users/nicholassiemons/Dropbox/OBT/0091/001/FES_D1/" + f for f in os.listdir("/Users/nicholassiemons/Dropbox/OBT/0091/001/FES_D1/")]
files_cm1 = ["/Users/nicholassiemons/Dropbox/OBT/0091/001/FES_CM1/" + f for f in os.listdir("/Users/nicholassiemons/Dropbox/OBT/0091/001/FES_CM1/")]
space.add_line(FreeEnergyLine.from_plumed(files_d1))
space.add_line(FreeEnergyLine.from_plumed(files_cm1))

path = "/Users/nicholassiemons/Dropbox/OBT/0091/001/"
files = [path + f for f in os.listdir(path) if 'COLVAR_REWEIGHT' in f and 'bck' not in f]
[space.add_metad_trajectory(MetaTrajectory(f)) for f in files]

data = space.get_line('CM4', bins=2)
data = data.data
print(data)
