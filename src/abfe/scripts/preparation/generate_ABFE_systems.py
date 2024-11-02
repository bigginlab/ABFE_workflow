import argparse
import ast
import glob
import os
import shutil
import subprocess
import sys

import BioSimSpace as bss

from abfe.scripts.preparation.parametrize import gen_ffparams
from abfe.scripts.preparation.topology_fixer import fix_topology


# from conda.base.context import locate_prefix_by_name


#########################################################################################33
def solvate(md_system, outdir: str):
    print("Solvating the system in: ", outdir)

    padding = 12 * bss.Units.Length.angstrom

    #box_min, box_max = md_system.getAxisAlignedBoundingBox()
    #base_length = max(2 * molecule._getAABox().halfExtents())
    #box_size = [y - x for x, y in zip(box_min, box_max)]
    #box_length = (max(box_size) + padding)
    #box, angles = bss.Box.truncatedOctahedron(box_length.value() * bss.Units.Length.angstrom)

    solvated = bss.Solvent.tip3p(md_system, shell=padding) #box=box, angles=angles)
    out_solv = outdir + "/solvated"
    bss.IO.saveMolecules(out_solv, solvated, ["GroTop", "Gro87"])


def prepare_md_system(protein="", ligand="", cofactor=""):
    if protein and ligand and cofactor:
        md_system = protein + ligand + cofactor
    elif protein and ligand and not cofactor:
        md_system = protein + ligand
    elif protein and not ligand and not cofactor:
        md_system = protein
    elif not protein and ligand and not cofactor:
        md_system = ligand
    else:
        print("Can not proceed with the provided inputs")
        sys.exit()

    return md_system


def process_ligand(ligand, out_dir: str, name="LIG", hmr: bool = True, ff="openff"):
    # Add a line to convert molecule to sdf by babel for safe processing
    tmp_mol = None
    if ".sdf" in ligand:
        ligand_filename = ligand.replace(".sdf", "")

        # Ensure sdf standard - for biosimspace
        #mol = next(Chem.SDMolSupplier(ligand, removeHs=False))
        #wri = Chem.SDWriter(ligand)
        #wri.write(mol)
        #time.sleep(5) #file latency

    elif ".mol2" in ligand:
        ligand_filename = ligand.replace(".mol2", "")
    else:
        ligand_filename = ligand

    if not os.path.isfile(os.path.join(ligand_filename, "for_gromacs", "MOL.top")):
        gen_ffparams(input_molecule=ligand, output_dir=out_dir, input_molecule_name=name, hmr=hmr, ff=ff)

    if(tmp_mol is not None):
        os.remove(tmp_mol)

    top_file = os.path.join(out_dir, "for_gromacs", "MOL.top")
    gro_file = os.path.join(out_dir, "for_gromacs", "MOL.gro")

    sys_ligand = bss.IO.readMolecules([top_file, gro_file])

    return sys_ligand, out_dir + "/for_gromacs/MOL"


def process_protein(protein_pdb_path, out_dir: str):
    if (not os.path.exists(out_dir)):
        os.mkdir(out_dir)

    env_prefix = os.environ["CONDA_PREFIX"]

    pdb_fixed = out_dir + '/' + os.path.basename(protein_pdb_path.replace(".pdb", "_fix.pdb"))
    subprocess.getoutput(f"{env_prefix}/bin/pdbfixer {protein_pdb_path} --output={pdb_fixed} --add-atoms=all --replace-nonstandard")

    subprocess.getoutput(
        f"gmx pdb2gmx -f {pdb_fixed} -merge all -ff amber99sb-ildn -water tip3p -o {out_dir}/conf.gro -p {out_dir}/topol.top -i {out_dir}/posre.itp -ignh")

    protein = bss.IO.readMolecules([out_dir + "/conf.gro", out_dir + "/topol.top"])
    protein.repartitionHydrogenMass(factor=3, water="no")

    bss.IO.saveMolecules(out_dir + "/heavy", protein, ["GroTop"])
    sys_protein = bss.IO.readMolecules([out_dir + "/conf.gro", out_dir + "/heavy.top"])

    return sys_protein


def prepare_for_abfe(out_dir, ligand_dir, sys_dir):
    complex_out = out_dir + "/complex"
    ligand_out = out_dir + "/ligand"
    if (not os.path.exists(complex_out)): os.makedirs(complex_out)
    if (not os.path.exists(ligand_out)): os.makedirs(ligand_out)

    shutil.copyfile(src=sys_dir + "/solvated.gro", dst=complex_out + "/complex.gro")
    shutil.copyfile(src=sys_dir + "/solvated_fix.top", dst=complex_out + "/complex.top")

    for itp_file in glob.glob(sys_dir + "/*.itp"):
        shutil.copy(src=itp_file, dst=complex_out)

    shutil.copyfile(src=ligand_dir + "/solvated.gro", dst=ligand_out + "/ligand.gro")
    shutil.copyfile(src=ligand_dir + "/solvated_fix.top", dst=ligand_out + "/ligand.top")

    for itp_file in glob.glob(ligand_dir + "/*.itp"):
        shutil.copy(src=itp_file, dst=ligand_out)


def clean_up(ligand="", cofactor=""):
    subprocess.getoutput("rm *#* *.itp *.top *.gro")
    if ligand:
        ligand_name = os.path.split(ligand)[-1]
        ligand_name = ligand_name.replace(".sdf", "")
        ligand_name = ligand_name.replace(".mol2", "")
        subprocess.getoutput(f"rm -fr {ligand_name}")
    if cofactor:
        cofactor_name = os.path.split(cofactor)[-1]
        cofactor_name = cofactor_name.replace(".sdf", "")
        cofactor_name = cofactor_name.replace(".mol2", "")
        subprocess.getoutput(f"rm -fr {cofactor_name}")


def prepare_input_files(protein_pdb: str, ligand_sdf: str, cofactor_sdf: str, out_dir: str, ff:str ="openff"):
    work_dir = out_dir + "/tmp"

    if (not os.path.exists(work_dir)):
        os.mkdir(work_dir)

    if ligand_sdf:
        print("Processing Ligand")
        sys_ligand, lig_file = process_ligand(ligand_sdf, out_dir=work_dir + "/ligand", name="LIG", ff=ff)
    else:
        sys_ligand = ""
    if cofactor_sdf:
        print("Processing Cofactor")
        sys_cofactor, _ = process_ligand(cofactor_sdf, out_dir=work_dir + "/cof", name="COF", ff=ff)
    else:
        sys_cofactor = ""
    if protein_pdb:
        print("Processing Protein")
        sys_protein = process_protein(protein_pdb, out_dir=work_dir + "/protein", )
    else:
        sys_protein = ""

    # Construct MD system:
    md_system = prepare_md_system(protein=sys_protein, ligand=sys_ligand, cofactor=sys_cofactor)
    system_dir = work_dir + "/system"

    solvate(md_system, outdir=system_dir)
    solvate(sys_ligand, outdir=work_dir + "/ligand")

    fix_topology(input_topology_path=system_dir + "/solvated.top", out_topology_path=system_dir + "/solvated_fix.top")
    fix_topology(input_topology_path=work_dir + "/ligand/solvated.top", out_topology_path=work_dir + "/ligand/solvated_fix.top")

    # Construct ABFE system:
    prepare_for_abfe(out_dir=out_dir, ligand_dir=work_dir + "/ligand", sys_dir=system_dir)
    # clean_up(ligand=ligand_sdf, cofactor=cofactor_sdf)


#############################################################################################

if __name__ == "__main__":
    # ARGPARSE
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', "--protein_pdb_path", help='Input protein file', required=True)
    parser.add_argument('-l', "--ligand_sdf_dir", help='Input ligand file', required=True)
    parser.add_argument('-c', "--cofactor_sdf_path", help='Input cofactor file', required=False, default=None)
    parser.add_argument('-o', "--output_dir_path", help='Output directory', required=False, default=".")
    parser.add_argument('-ff', "--smallmol_ff", help='Force Field small molecule parametrization', required=False, default="openff")

    args = parser.parse_args()

    print("Input: ", args.protein_pdb_path, args.ligand_sdf_dir, args.cofactor_sdf_path, args.output_dir_path)

    prepare_input_files(protein_pdb=args.protein_pdb_path, ligand_sdf=args.ligand_sdf_dir, cofactor_sdf=ast.literal_eval(args.cofactor_sdf_path),
                        out_dir=args.output_dir_path, ff=args.smallmol_ff)
