"""
calculate dF for complex system
"""

import json
import argparse
import numpy as np
from abfe.scripts.alchemlyb_analysis import analyze_ligand


if __name__ == "__main__":
    # results list
    results = []

    parser = argparse.ArgumentParser()
    parser.add_argument('--xvgpath', default='./',
                        help='input xvg path')
    parser.add_argument('--confpath', default='../../../snake_conf.json',
                        help='input xvg path')
    parser.add_argument('--outpath', default='./',
                        help='output path for writing files')
    parser.add_argument("--boresch_data",
                        help='boresch restraint correction')
    args = parser.parse_args()


    conf = json.load(open(args.confpath,"r"))
    system_steps_windows = {
        "restraints-xvg": conf['n_rest_windows_complex'],
        'vdw-xvg':conf['n_vdw_windows_complex'],
        'coul-xvg': conf['n_coul_windows_complex']
    }

    res_df = analyze_ligand(prefix=args.xvgpath, 
                            system_steps_windows=system_steps_windows,
                            system_name="complex",
                            lower=1000,
                            min_samples=200)
    
    #include boresch correction
    col = []
    res_V = float(np.loadtxt(args.boresch_data))
    for row in res_df.index:
        if(row == "sys"):
            col.append("ligand")
        elif(row == "windows"):
            col.append("-")
        else:
            col.append(res_V)
    res_df['boresch'] = col
    
    res_df.to_csv(args.outpath+"/dg_results.csv")
