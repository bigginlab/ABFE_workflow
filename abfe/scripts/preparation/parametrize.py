#!/home/magarkar/anaconda3/envs/openff/bin/python

import argparse,os,subprocess
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import get_data_file_path
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.interchange import Interchange
from openff.toolkit.topology import Topology
import parmed as pmd

def parameterize(input_molecule, output_dir, hmr=False,input_molecule_name="MOL"):
    print("------------------------------------------------------")
    print("processing : %s"%input_molecule)
    print("setting HMR to : %s"%hmr)
    print("Molecule name will be set to : %s"% input_molecule_name)
    print("------------------------------------------------------")
    
    # Make a directory for molecule
    dir_name =output_dir
    #if ".sdf" in input_molecule: 
    #    dir_name = input_molecule.replace(".sdf","")
    #elif ".mol2" in input_molecule:
    #    dir_name = input_molecule.replace(".mol2","")
    #dir_name = os.path.split(dir_name)[-1]
    
    #subprocess.getoutput(f"rm -fr {dir_name}")
    if(not os.path.exists(output_dir)):
        subprocess.getoutput(f"mkdir {dir_name}")

    sdf_file_path = (input_molecule)
    molecule: Molecule = Molecule.from_file(sdf_file_path)
    topology: Topology = molecule.to_topology()

    sage = ForceField("openff-2.0.0.offxml")

    interchange = Interchange.from_smirnoff(force_field=sage, topology=topology)
    interchange.positions = molecule.conformers[0]

    #openmm_system = interchange.to_openmm()

    interchange.to_gro(os.path.join(dir_name,"out.pdb"))
    interchange.to_gro(os.path.join(dir_name,"out.gro"))
    interchange.to_top(os.path.join(dir_name,"out.top"))
    interchange.to_prmtop(os.path.join(dir_name,"out.prmtop"))

    pmd_top = pmd.load_file(os.path.join(dir_name,"out.prmtop"))
    pmd_pdb = pmd.load_file(os.path.join(dir_name,"out.pdb"))

    for atom in pmd_pdb.atoms:
        atom.residue.name=input_molecule_name
    for atom in pmd_top.atoms:
        atom.residue.name=input_molecule_name

    if hmr:
        pmd.tools.HMassRepartition(pmd_top,3).execute()

    subprocess.getoutput("mkdir %s"%os.path.join(dir_name,"for_gromacs"))
    subprocess.getoutput("mkdir %s"%os.path.join(dir_name,"for_amber"))

    pmd_top.save(os.path.join(dir_name,"for_gromacs","MOL.top"),overwrite=True)
    pmd_pdb.save(os.path.join(dir_name,"for_gromacs","MOL.pdb"),overwrite=True)
    pmd_pdb.save(os.path.join(dir_name,"for_gromacs","MOL.gro"),overwrite=True)

    pmd_top.save(os.path.join(dir_name,"for_amber","MOL.prmtop"),overwrite=True)
    pmd_pdb.save(os.path.join(dir_name,"for_amber","MOL.inpcrd"),overwrite=True)
    pmd_pdb.save(os.path.join(dir_name,"for_amber","MOL.pdb"),overwrite=True)

def gen_openff(input_molecule:str, output_dir:str, input_molecule_name:str="LIG", hmr:bool=False):
    
    if ".sdf" in input_molecule or ".mol2" in input_molecule:
        print(input_molecule,input_molecule_name, hmr)
        
        if os.path.isfile(input_molecule):
            parameterize(input_molecule=input_molecule, output_dir=output_dir,hmr=hmr,input_molecule_name=input_molecule_name) 
        else:
            print(f"Input file not valid : {input_molecule}")
            raise ValueError(f"Input file not valid : {input_molecule}")
            parameterize(input_molecule,hmr=hmr,input_molecule_name=input_molecule_name) 
    else:
        raise IOError("input mol must be of type .sdf or .mol2")

        

    
    
#######################################################################################################################


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='OpenFF topology generator for small molecules',
    )

    parser.add_argument('-i',type=str,default='arg_default',nargs='?',help='input_molecule_file',required=True)
    parser.add_argument('-o',type=str,default='./',nargs='?',help='output dir',required=True)

    parser.add_argument('-mol_name',type=str,default='LIG',help='Molecule Name', required=False)
    
    parser.add_argument('-hmr',type=bool,default=False, help='Makes hydrogens heavy', required=False)

    args = parser.parse_args()

    gen_openff(input_molecule=args.i, output_dir=args.o, input_molecule_name=args.mol_name, hmr=args.hmr)

