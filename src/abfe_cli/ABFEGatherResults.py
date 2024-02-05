#!/usr/bin/env python3

import os
import argparse
import logging


loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.NOTSET)

from abfe.scripts import final_receptor_results


def main():
    # ARGPARSE
    parser = argparse.ArgumentParser(description="this tool can be used to gather final or temporary results of the ABFE workflow.")
    parser.add_argument('-i', "--in_dir", help='Input receptor folder, containing ligand approaches', required=True, type=str)
    parser.add_argument('-o', "--out_dir", help='Output approach folder for free energy .tsvs', required=True, type=str)

    args = parser.parse_args()

    if os.path.isdir(args.in_dir):
        in_dir = str(args.in_dir)
    else:
        raise IOError("could not find in_dir: " + str(args.in_dir))
    if os.path.isdir(args.out_dir):
        out_dir = str(args.out_dir)
    else:
        raise IOError("could not find out_dir: "+str(args.out_dir))

    print("Trying to gather ready results from ", in_dir)
    print("Searching: ")
    out_df_final_results, out_df_single_detailed_results =    final_receptor_results.get_final_results(out_dir=out_dir, in_root_dir=in_dir)
    print()
    print("Get the final results from here: ", out_df_final_results)

if __name__ == "__main__":
    main()
