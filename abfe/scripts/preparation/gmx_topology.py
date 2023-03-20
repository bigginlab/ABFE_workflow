#!/user/bin python
import os
from typing import Union, Iterable
PathLike = Union[os.PathLike, str, bytes]

def get_molecule_names(input_topology:PathLike) -> list:
    """IT get the molecules names specified inside input_topology

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology

    Returns
    -------
    list
        A list of the molecules presented in the topology
    """
    molecules_flag = False
    molecules = []
    with open(input_topology, "r") as topology_file:
        for line in topology_file.readlines():
            if '[ molecules ]' in line:
                molecules_flag = True
            if molecules_flag:
                if not "[" in line and not line.startswith(";"):
                    if len(line.split()) > 1:
                        molecules.append(line.split()[0])
    return molecules

def add_posres_section(input_topology:PathLike, molecules:Iterable[str], out_file:PathLike="topol2.top"):
    """This will add to the original topology file the corresponded POSRES section to the 
    provided molecules:
    Examples of added lines:

    #ifdef POSRES
    #include "posres_{molecule}.itp"
    #endif


    #ifdef POSRES_{molecule}
    #include "posres_{molecule}.itp"
    #endif

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology
    molecules : Iterable[str]
        The list of name of the molecules for which the topology section will be add
    out_file : PathLike, optional
        The path to output the modified topology file, by default "topol2.top"
    """
    with open(input_topology, "r") as f:
        top_lines = f.readlines()

    look_out_flag = False
    out_line = []
    for line in top_lines:
        for molecule in molecules:
            if f"{molecule}" in line and " 3\n" in line:
                look_out_flag = True
                mol_name = line.split()[0]
            if look_out_flag and ('[ moleculetype ]' in line or '[ system ]' in line):
                out_line.append("\n#ifdef POSRES\n")
                out_line.append(f'#include "posres_{mol_name}.itp"\n')
                out_line.append("#endif\n\n")
                # if ("LIG" in mol_name or "MOL" in mol_name):
                #     out_line.append("\n#ifdef POSRES_LIG\n")
                #     out_line.append(f'#include "posres_{mol_name}.itp"\n')
                #     out_line.append("#endif\n\n")
                look_out_flag = False
        out_line.append(line)

    with open(out_file, "w") as w:
        for line in out_line:
            w.write(line)

def make_posres_files(input_topology:PathLike, molecules:Iterable[str], out_dir:PathLike):
    """Make a position restraint file out of input_topology for all the molecules specified
    on molecules. The force constant used will be: fx = 2500; fy = 2500; fz = 2500
    for all atoms (included hydrogens, TODO: Check if this is what is expected, or only on heavy atoms)

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology
    molecules : Iterable[str]
        The list of name of the molecules for which the posres file will be created
    out_dir : PathLike
        The path where the posres files will be written
    """
    for molecule in molecules:
        atom_flag = False
        bonds_flag = False

        with open(input_topology, "r") as f:
            top_lines = f.readlines()

        posres_filename = f"posres_{molecule}.itp"
        with open(os.path.join(out_dir,posres_filename), "w") as posres_file:
            posres_file.write("[ position_restraints ]\n")

            for i in range(len(top_lines)):
                if f"{molecule}  " in (top_lines[i]) and " 3\n" in (top_lines[i]):
                    atom_flag = False
                    bonds_flag = False

                    for j in range(i + 1, len(top_lines)):
                        if '[ atoms ]' in top_lines[j]:
                            atom_flag = True
                        if '[ bonds ]' in top_lines[j]:
                            bonds_flag = True
                            break
                        if atom_flag and not bonds_flag:
                            if not "[" in top_lines[j] and not top_lines[j].startswith("\n") and not top_lines[j].startswith(";"):
                                if float(top_lines[j].split()[7]) > 3:
                                    posres_str = f"{top_lines[j].split()[0]} 1 2500 2500 2500\n"
                                    posres_file.write(posres_str)

def fix_topology(input_topology: PathLike, out_dir: PathLike, exclusion_list:list = ["SOL", "NA", "CL", "MG", "ZN"]):
    """It will go through input_topology, create the posres files for the identified molecules
    not in exclusion list, and finally add the corresponded include statements in the topology file.
    on out_topology_path you will have
    It will call, sequentially, to:

    #. :meth:`abfe.scripts.gmx_topology.get_molecule_names`
    #. :meth:`abfe.scripts.gmx_topology.make_posres_files`
    #. :meth:`abfe.scripts.gmx_topology.add_posres_section`

    In out_dir you will have:

    #. Fixed topology file name_of_input_topology_fix.top
    #. Posres itp files.

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology
    out_dir : PathLike
        Where the fixed and generated itp posres files will be written
    exclusion_list : list, optional
        Molecules names to do not take into account during the posres file generation, by default ["SOL", "NA", "CL", "MG", "ZN"]
    """
    name, ext = os.path.splitext(os.path.basename(input_topology))
    out_topology = os.path.join(out_dir, f"{name}_fix{ext}")
    print(get_molecule_names(input_topology))
    molecules = list(set(get_molecule_names(input_topology)) - set(exclusion_list))
    print(molecules)
    make_posres_files(input_topology, out_dir = out_dir, molecules = molecules)
    add_posres_section(input_topology, molecules, out_file=out_topology)


if __name__ == "__main__":...

    # # input_topology_path = sys.argv[1]
    # # out_topology_path = sys.argv[2]
    # exclusion_list = ["SOL", "NA", "CL", "MG", "ZN"]
    # fix_topology(
    #     input_topology = '/home/users/alejandro/GIT/ABFE_workflow/examples/prepearing_system/builder/system/solvated.top',
    #     out_dir='test',
    #     exclusion_list=exclusion_list,

    # )

    # # fix_topology(input_topology_path=input_topology_path, out_topology_path=out_topology_path, exclusion_list=exclusion_list)