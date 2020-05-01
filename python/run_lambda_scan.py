#!/usr/bin/env python
"""
Script that launches a batch job for a lambda_2 vs lambda_1 scan
"""
import os
import subprocess
import numpy as np


BATCH_SCRIPT = os.path.join(os.environ['WORK'], 'chic_pol_global_fit',
                            'scripts', 'batch_run_fit.sh')
DATA_DIR = os.path.join(os.environ['WORK'], 'chic_pol_global_fit', 'data')

def get_batch_base_command(outdir, use_costh, use_psi, unphys_corrs):
    """Get the base of all batch commands"""
    return [
        "sbatch", BATCH_SCRIPT,
        os.path.join(outdir, "lambda_scan.root"), # --outfile is added by batch script
        "--datadir", DATA_DIR,
        "--nographs", "true", "--randomScan", "false",
        "--useCosthRatios", str(use_costh).lower(),
        "--usePsiPolarizations", str(use_psi).lower(),
        "--usePhysicalCorrections", str(unphys_corrs).lower(),
    ]


def get_scan_commands(min_l1, max_l1, min_l2, max_l2, delta_l, n_jobs):
    """Get the parts of the command that specify which lambda ranges to scan"""
    # Scan the whole lambda_1 range in thin slices of lambda_2
    n_points_1 = int(round((max_l1 - min_l1) / delta_l))
    if min_l1 + delta_l * n_points_1 != max_l1:
        print('Range for lambda 1 is [{}, {}], grid spacing is {}. Cannot '
              'determine an integer number of scan points. Actual grid '
              'spacing will be slightly different'.format(min_l1, max_l1, delta_l))

    lambda_1_com = [
        "--flow1", str(min_l1),
        "--fhigh1", str(max_l1),
        "--nscan1", str(n_points_1 + 1) # include the end-points
    ]

    n_points_2 = int(round((max_l2 - min_l2) / delta_l))
    if min_l2 + delta_l * n_points_2 != max_l2:
        print('Range for lambda 2 is [{}, {}], grid spacing is {}. Cannot '
              'determine an integer number of scan points. Actual grid '
              'spacing will be slightly different'.format(min_l2, max_l2, delta_l))

    sub_coms = []
    for scan_vals in np.array_split(np.linspace(min_l2, max_l2, n_points_2 + 1), n_jobs):
        job_min = np.min(scan_vals)
        job_max = np.max(scan_vals)
        n_steps = len(scan_vals)

        sub_coms.append(
            lambda_1_com +
            ["--flow2", str(job_min), "--fhigh2", str(job_max), "--nscan2", str(n_steps)]
        )

    return sub_coms


def main(args):
    """Main"""
    base_command = get_batch_base_command(args.outdir,
                                          not args.no_costh_ratios,
                                          not args.no_psi_polarization,
                                          args.use_physical_corr)

    scan_commands = get_scan_commands(args.lambda_low_1, args.lambda_high_1,
                                      args.lambda_low_2, args.lambda_high_2,
                                      args.delta_lambda, args.n_jobs)

    for scan_job in scan_commands:
        full_com = base_command + scan_job
        if args.dryrun:
            print(' '.join(full_com))
        else:
            subprocess.call(full_com)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Scrip to launch a lambda scan '
                                     'batch job')
    parser.add_argument('outdir', help='Directory into which the resulting files'
                        ' from the batch job will be put')
    parser.add_argument('--no-costh-ratios', help='Do not the chic costh ratio '
                        'constraints', default=False, action='store_true')
    parser.add_argument('--no-psi-polarization', help='Do not use the psi '
                        'polarization constraints', default=False,
                        action='store_true')
    parser.add_argument('--use-physical-corr', help='Use corrections '
                        'calculated from lambda values outside the physically '
                        'allowed region', default=False, action='store_true')
    parser.add_argument('-l1', '--lambda-low-1', help='Minimal value of lambda 1'
                        ' to be used in the scanning', default=-1.0,
                        type=float)
    parser.add_argument('-h1', '--lambda-high-1', help='Maximal value of lambda '
                        '1 to be used in the scanning', default=2.0,
                        type=float)
    parser.add_argument('-l2', '--lambda-low-2', help='Minimal value of lambda 2'
                        ' to be used in the scanning', default=-1.5,
                        type=float)
    parser.add_argument('-h2', '--lambda-high-2', help='Maximal value of lambda '
                        '2 to be used in the scanning', default=2.0,
                        type=float)
    parser.add_argument('-dl', '--delta-lambda', help='Distance between two scan'
                        ' points', type=float, default=0.01)
    parser.add_argument('-n', '--n-jobs', help='Number of batch jobs to create',
                        type=int, default=50)
    parser.add_argument('--dryrun', help='Only print the commands that would be '
                        'launched instead of actually launching them',
                        action='store_true', default=False)


    clargs = parser.parse_args()
    main(clargs)
