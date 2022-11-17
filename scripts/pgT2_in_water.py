import MDAnalysis as mda

# Initiate Universe
u = mda.Universe("/Users/nicholassiemons/Dropbox/OBT/0083/000/EM.tpr")

polymer = u.atoms.select_atoms('not resname SOL or resname NA or resname CL')
print(polymer.atoms)
