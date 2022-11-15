import MDAnalysis as mda

# make the universe
u = mda.Universe("../test_trajectories/alk_mon_in_water/Prod.gro")

# get the number of atoms and the number of residues in the system
print(f"Number of atoms is {len(u.atoms)}")
print(f"Number of residues is {len(u.residues)}")

# get the residue names and residue id's in the universe
res_names = u.atoms.residues.resnames
res_ids = u.atoms.residues.resids
print(f"Residue names are {res_names}")
print(f"Residue ids are {res_ids}")

# selection language trials
print(len(u.select_atoms("name OW")))
print(len(u.select_atoms("name CS1 CS2 CS3 CS4 CD1 CD2 CD3 CD4")))
print(len(u.select_atoms("name H*")))
