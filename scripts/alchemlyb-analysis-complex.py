import math
import warnings
import glob
import csv
import argparse

import pandas as pd
import numpy as np

from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
from alchemlyb.preprocessing import statistical_inefficiency, slicing
from alchemlyb.estimators import TI, MBAR


def run_alchemlyb(xvgs: list, overlap_path: str = None, lower: int = None,
                  upper: int = None, min_samples: int = 500,
                  temperature: float = 298.15):
    """
    Function to get MBAR and TI estimates using alchemlyb from an input set of
    xvgs

    Parameters
    ----------
    xvgs : list of str
        list of filenames for input xvg files
    overlap_path : str
        path to write overlap matrix (if None, no matrix will be written)
        [None]
    lower : int
        starting time to sample dhdl xvgs from [None]
    upper : int
        inclusive end time to sample dhdl xvgs from [None]
    min_samples : int
        minimum number of samples to analyze, if statistical inefficiency
        returns fewer samples than this, then min_samples will be picked
        instead [500]
    temperature : float
        simulation temperature [298.15]

    Returns
    -------
    deltaG : dict
        two entry dictionary containing the MBAR and TI free energy and
        associated variance error estimate in kcal/mol
    """

    # first get the cutoff values by running statistical inefficiency on
    # all the data points and checking if it returns fewer or more data points
    # than min_samples. For consistency, we run the statistical inefficiency
    # on the dHdl values of the lambda that is changing (first entry)
    sub_steps = []

    for xvg in xvgs:
        df = slicing(extract_dHdl(xvg, T=temperature), lower=lower,
                     upper=upper)
        df_ineff = statistical_inefficiency(df, series=df.iloc[:, 0])

        ineff_step = math.ceil(len(df) / len(df_ineff))
        step_cutoff = int(len(df) / min_samples)

        if ineff_step <= step_cutoff:
            sub_steps.append(ineff_step)
        else:
            sub_steps.append(int(step_cutoff))

    print(f"number of samples per window: {[int(len(df)/i) for i in sub_steps]}")

    dhdls = pd.concat([slicing(extract_dHdl(xvg, T=temperature), lower=lower,
                               upper=upper, step=step)
                       for xvg, step in zip(xvgs, sub_steps)])
    u_nks = pd.concat([slicing(extract_u_nk(xvg, T=temperature), lower=lower,
                               upper=upper, step=step)
                       for xvg, step in zip(xvgs, sub_steps)])

    ti = TI().fit(dhdls)
    mbar = MBAR(maximum_iterations=1000000).fit(u_nks)

    deltaG = {'MBAR': (mbar.delta_f_.iloc[0, -1] * 0.593,
                       mbar.d_delta_f_.iloc[0, -1] * 0.593),
              'TI': (ti.delta_f_.iloc[0, -1] * 0.593,
                     ti.d_delta_f_.iloc[0, -1] * 0.593)}

    return deltaG


def analyze_complex(prefix: str,
                    xvg_prefix: str = 'dhdl',
                    lower: int = None, upper: int = None,
                    min_samples: int = None,
                    temperature: float = 298.15):
    """
    Function to run an FEP analysis for an FEP cycle

    Parameters
    ----------
    prefix : str
        Path to complex.
    xvg_prefix : str
        The prefix name for the sequentially numbered xvg files.
    lower : int
        Time at which to start sampling xvg files [None]
    upper : int
        Time at which to stop sampling xvg files [None]
    min_sample : int
        Minimum number of samples to use in each xvg file (post-subsampling)
        [None]
    temperature : float
        Simulation temperature [298.15]
    """

    complex_steps = ['restraints-xvg', 'vdw-xvg', 'coul-xvg']
    complex_windows = [12, 21, 11]
    complex_results = []

    for step, windows in zip(complex_steps, complex_windows):
        xvgs = [f'{prefix}/{step}/{xvg_prefix}.{i}.xvg'
                for i in range(windows)]
        dG = run_alchemlyb(xvgs, lower=lower, upper=upper,
                           min_samples=min_samples, temperature=temperature)
        ddG_estimator = abs(dG['MBAR'][0] - dG['TI'][0])
        if ddG_estimator > 0.5:
            wmsg = (f'ddG_estimator > 0.5 kcal/mol: {ddG_estimator} '
                    f'{prefix}/{step}')
            warnings.warn(wmsg)

        complex_results.append(dG['MBAR'][0])

    print("test")
    print(complex_results)
    return complex_results


if __name__ == "__main__":
    # results list
    results = []

    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--xvgpath', default='./',
                            help='input xvg path')
        parser.add_argument('--outpath', default='./',
                            help='output path for writing files')
        args = parser.parse_args()
        return args

    args = parse_args()

    results.append(analyze_complex(prefix=args.xvgpath, lower=1000,
                                   min_samples=200))

    # Print everything in a nice csv formatted way
    with open(f'{args.outpath}/dg_results.csv', 'w', newline='\n') as csvfile:
        csvwriter = csv.writer(csvfile)
        for entry in allruns_results:
            csvwriter.writerow(entry)
