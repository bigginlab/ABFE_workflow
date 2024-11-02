import argparse

import pandas as pd

"""
calculate dF for full system
"""

def main():
    # results list
    results = []

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_lig_path', default='./',
                        help='results from ligand.tsv')
    parser.add_argument('--in_comp_path', default='',
                        help='results from complex.tsv')
    parser.add_argument('--out_csv_path', default='./',
                        help='output path for writing files')
    args = parser.parse_args()

    out_csv_path = args.out_csv_path
    in_complex_path = args.in_comp_path
    in_ligand_path = args.in_lig_path

    # read in csv:
    complex_df = pd.read_csv(in_complex_path, index_col=0, )
    ligand_df = pd.read_csv(in_ligand_path, index_col=0)

    # Create new dict containing the results and calculate complete process:
    complex_dict = complex_df.to_dict()
    ligand_dict = ligand_df.to_dict()

    dG_totLig_MBAR = -float(ligand_dict["vdw"]["MBAR"]) - float(ligand_dict["coul"]["MBAR"]) + float(complex_dict["boresch"]["MBAR"])
    dG_totComp_MBAR = -float(complex_dict["restraints"]["MBAR"]) - float(complex_dict["vdw"]["MBAR"]) - float(complex_dict["coul"]["MBAR"])

    dG_totLig_TI = -float(ligand_dict["vdw"]["TI"]) - float(ligand_dict["coul"]["TI"]) + float(complex_dict["boresch"]["TI"])
    dG_totComp_TI = -float(complex_dict["restraints"]["TI"]) - float(complex_dict["vdw"]["TI"]) - float(complex_dict["coul"]["TI"])

    complex_dict.update({"total": {"sys": "complex", "windows": "-", "MBAR": dG_totComp_MBAR, "TI": dG_totComp_TI, }})
    ligand_dict.update({"total": {"sys": "ligand", "windows": "-", "MBAR": dG_totLig_MBAR, "TI": dG_totLig_TI, }})

    collapse = []
    for k, v in ligand_dict.items():
        v.update({"step": k})
        collapse.append(v)

    for k, v in complex_dict.items():
        v.update({"step": k})
        collapse.append(v)

    collapse.append({"step": "ABFE", "sys": "ABFE", "windows": "-", "MBAR": dG_totComp_MBAR - dG_totLig_MBAR, "TI": dG_totComp_TI - dG_totLig_TI, })

    # convert to pandas df.
    fin_df = pd.DataFrame(collapse).sort_values("sys")
    fin_df = fin_df.reset_index(drop=True)
    fin_df["MBAR"] = fin_df["MBAR"].astype(float)
    fin_df["TI"] = fin_df["TI"].astype(float)

    fin_df = fin_df.round(2)
    fin_df.to_csv(out_csv_path)


if __name__ == "__main__":
    main()