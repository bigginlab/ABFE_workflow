#!/user/bin python
import os
import sys


def add_it_in_topology_file(input_topology, molecules, out_file="topol2.top"):
    out_line = []

    topology_file = open(input_topology, "r")
    top_lines = topology_file.readlines()

    look_out_flag = 0

    for line in top_lines:
        for molecule in molecules:
            if f"{molecule}" in line and " 3\n" in line:
                look_out_flag = 1
                mol_name = line.split()[0]
            if look_out_flag == 1 and ('[ moleculetype ]' in line or '[ system ]' in line):
                out_line.append("\n#ifdef POSRES\n")
                out_line.append('#include "posres_%s.itp"\n' % mol_name)
                out_line.append("#endif\n\n")
                if ("LIG" in mol_name or "MOL" in mol_name):
                    out_line.append("\n#ifdef POSRES_LIG\n")
                    out_line.append('#include "posres_%s.itp"\n' % mol_name)
                    out_line.append("#endif\n\n")
                look_out_flag = 0

        out_line.append(line)

    w = open(out_file, "w")
    for line in out_line:
        w.write(line)
    w.close()


def write_posres_files(input_topology, out_dir: str, molecules=[]):
    for molecule in molecules:
        print(molecule)
        molecule_flag = 0
        atom_flag = 0
        bonds_flag = 0

        topology_file = open(input_topology, "r")
        top_lines = topology_file.readlines()

        posres_filename = 'posres_' + molecule + ".itp"
        posres_file = open(out_dir + "/" + posres_filename, "w")
        posres_file.write("[ position_restraints ]\n")

        for i in range(0, len(top_lines)):
            if f"{molecule}  " in (top_lines[i]) and " 3\n" in (top_lines[i]):
                atom_flag = 0
                bonds_flag = 0

                for j in range(i + 1, len(top_lines)):

                    if '[ atoms ]' in top_lines[j]:
                        atom_flag = 1
                    if '[ bonds ]' in top_lines[j]:
                        bonds_flag = 1
                        break
                    if atom_flag == 1 and bonds_flag == 0:
                        if not "[" in top_lines[j] and not top_lines[j].startswith("\n") and not top_lines[j].startswith(";"):
                            if float(top_lines[j].split()[7]) > 3:
                                posres_str = top_lines[j].split()[0] + " 1 2500 2500 2500\n"
                                posres_file.write(posres_str)

        posres_file.close()


def what_are_the_molecules(input_topology, exclusion_list=[]):
    molecules_flag = 0
    molecules = []
    topology_file = open(input_topology, "r")

    for line in topology_file.readlines():
        if '[ molecules ]' in line:
            # print(line)
            molecules_flag = 1
        if molecules_flag == 1:
            if not "[" in line and not line.startswith(";"):
                if len(line.split()) > 1:
                    if line.split()[0] not in exclusion_list:
                        molecules.append(line.split()[0])

    return (molecules)


def fix_topology(input_topology_path: str, out_topology_path: str, exclusion_list=["SOL", "NA", "CL", "MG", "ZN"]):
    molecules = what_are_the_molecules(input_topology_path, exclusion_list)

    print(molecules)
    out_dir = os.path.dirname(out_topology_path)
    write_posres_files(input_topology_path, out_dir=out_dir, molecules=molecules)
    add_it_in_topology_file(input_topology_path, molecules, out_file=out_topology_path)


if __name__ == "__main__":
    input_topology_path = sys.argv[1]
    out_topology_path = sys.argv[2]
    exclusion_list = ["SOL", "NA", "CL", "MG", "ZN"]

    fix_topology(input_topology_path=input_topology_path, out_topology_path=out_topology_path, exclusion_list=exclusion_list)

    # Calls
