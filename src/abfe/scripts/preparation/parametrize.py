import argparse
import os
import shutil
import subprocess

from rdkit import Chem
import parmed as pmd
import BioSimSpace as BSS

from openff.interchange import Interchange
from openff.toolkit.topology import Molecule
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff import ForceField


def parameterize(input_molecule, output_dir, hmr=False, input_molecule_name="MOL", ff="openff"):
    print("------------------------------------------------------")
    print("processing : %s" % input_molecule)
    print("setting HMR to : %s" % hmr)
    print("Molecule name will be set to : %s" % input_molecule_name)
    print("------------------------------------------------------")

    # Make a directory for molecule
    dir_name = output_dir
    if (not os.path.exists(output_dir)):
        subprocess.getoutput(f"mkdir {dir_name}")
    # if ".sdf" in input_molecule:
    #    dir_name = input_molecule.replace(".sdf","")
    # elif ".mol2" in input_molecule:
    #    dir_name = input_molecule.replace(".mol2","")
    # dir_name = os.path.split(dir_name)[-1]

    # subprocess.getoutput(f"rm -fr {dir_name}")
    if(ff=="openff"):
        sdf_file_path = (input_molecule)
        molecule: Molecule = Molecule.from_file(sdf_file_path,"sdf")
        topology: Topology = molecule.to_topology()

        sage = ForceField("openff-2.0.0.offxml")

        interchange = Interchange.from_smirnoff(force_field=sage, topology=topology)
        interchange.positions = molecule.conformers[0]

        # openmm_system = interchange.to_openmm()
        pmd_pdb_file = os.path.join(dir_name, "out.pdb")
        pmd_top_file = os.path.join(dir_name, "out.prmtop")

        interchange.to_gro(pmd_pdb_file)
        interchange.to_gro(os.path.join(dir_name, "out.gro"))
        interchange.to_top(os.path.join(dir_name, "out.top"))
        interchange.to_prmtop(pmd_top_file)
	
	# Define the path to the gro file
        gro_coord_file = os.path.join(dir_name, "out.gro")
    elif(ff=="gaff"):
        # input
        system = BSS.IO.readMolecules(input_molecule)
        mol = system.getMolecules()[0]

        # Load a molecule from file.
        process = BSS.Parameters.gaff(mol)
        molecule = process.getMolecule()

        pmd_top_file, pmd_rst_file, pmd_pdb_file, gro_coord_file, gro_top_file = BSS.IO.saveMolecules(dir_name + "/out", molecule,
                                                                                                      ["prm7", "rst7", "pdb", "Gro87", "GroTop"])
    else:
        raise ValueError("I don't know this FF!")

    pmd_top = pmd.load_file(pmd_top_file)
    pmd_pdb = pmd.load_file(gro_coord_file)

    for atom in pmd_pdb.atoms:
        atom.residue.name = input_molecule_name
    for atom in pmd_top.atoms:
        atom.residue.name = input_molecule_name
    if hmr:
        pmd.tools.HMassRepartition(pmd_top, 3).execute()

    subprocess.getoutput("mkdir %s" % os.path.join(dir_name, "for_gromacs"))
    subprocess.getoutput("mkdir %s" % os.path.join(dir_name, "for_amber"))
    pmd_pdb.box = [1, 1, 1, 90, 90, 90] #dummy box

    pmd_top.save(os.path.join(dir_name, "for_gromacs", "MOL.top"), overwrite=True)
    pmd_pdb.save(os.path.join(dir_name, "for_gromacs", "MOL.pdb"), overwrite=True)
    pmd_pdb.save(os.path.join(dir_name, "for_gromacs", "MOL.gro"), overwrite=True)

    pmd_top.save(os.path.join(dir_name, "for_amber", "MOL.prmtop"), overwrite=True)
    pmd_pdb.save(os.path.join(dir_name, "for_amber", "MOL.inpcrd"), overwrite=True)
    pmd_pdb.save(os.path.join(dir_name, "for_amber", "MOL.pdb"), overwrite=True)


def gen_ffparams(input_molecule: str, output_dir: str, input_molecule_name: str = "LIG", hmr: bool = False, ff="openff"):
    if ".sdf" in input_molecule or ".mol2" in input_molecule:
        print(input_molecule, input_molecule_name, hmr)

        if os.path.isfile(input_molecule):
            parameterize(input_molecule=input_molecule, output_dir=output_dir, hmr=hmr, input_molecule_name=input_molecule_name, ff=ff)
        else:
            print(f"Input file not valid : {input_molecule}")
            raise ValueError(f"Input file not valid : {input_molecule}")
            parameterize(input_molecule, hmr=hmr, input_molecule_name=input_molecule_name)
    else:
        raise IOError("input mol must be of type .sdf or .mol2")


#######################################################################################################################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='OpenFF topology generator for small molecules',
    )

    parser.add_argument('-i', type=str, default='arg_default', nargs='?', help='input_molecule_file', required=True)
    parser.add_argument('-o', type=str, default='./', nargs='?', help='output dir', required=True)

    parser.add_argument('-mol_name', type=str, default='LIG', help='Molecule Name', required=False)

    parser.add_argument('-hmr', type=bool, default=False, help='Makes hydrogens heavy', required=False)

    args = parser.parse_args()

    gen_ffparams(input_molecule=args.i, output_dir=args.o, input_molecule_name=args.mol_name, hmr=args.hmr)
