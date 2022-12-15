import MDAnalysis as mda
from MDRestraintsGenerator import search, restraints

assert mda.__version__ == '1.1.1'

u = mda.Universe('../npt_prod1/npt_prod1.tpr', 'npt_prod1_center.xtc')

ligand_atoms = search.find_ligand_atoms(u, l_selection="resname LIG",
                                        p_align="protein and name CA C N")

# find protein atoms
atom_set = []
for l_atoms in ligand_atoms:
    psearch = search.FindHostAtoms(u, l_atoms[0],
                                   p_selection="protein and name CA")
    psearch.run(verbose=True)
    atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])

# Create the boresch finder analysis object
boresch = restraints.FindBoreschRestraint(u, atom_set)

boresch.run(verbose=True)

boresch.restraint.plot()
boresch.restraint.write()

dG = boresch.restraint.standard_state()

with open('dG_off.dat', 'w') as writer:
    writer.write(f'{dG}')

