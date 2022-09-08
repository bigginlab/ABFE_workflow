"""
calculate dF for complex system
"""

import json
import argparse

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
    args = parser.parse_args()


    conf = json.load(open(args.confpath,"r"))
    system_steps_windows = {
        "restraints-xvg": conf['n_rest_windows_complex'],
        'vdw-xvg':conf['n_vdw_windows_ligand'],
        'coul-xvg': conf['n_coul_windows_ligand']
    }

    res_df = analyze_ligand(prefix=args.xvgpath, 
                            system_steps_windows=system_steps_windows,
                            system_name="complex",
                            lower=1000,
                            min_samples=200)
    
    res_df.to_csv(args.outpath+"/dg_results.csv")
