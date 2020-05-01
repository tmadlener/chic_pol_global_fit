#!/bin/bash
#SBATCH -J CorrelationScan
#SBATCH -D /afs/hephy.at/work/t/tmadlener/chic_pol_global_fit/
#SBATCH -o /afs/hephy.at/work/t/tmadlener/chic_pol_global_fit/batch_logfiles/run_correlation_scan_%A_%a.out

source ${CHIB_CHIC_POLFW_DIR}/scripts/bash_functions.sh

set -x

## input arguments
scanfile=${1}
fitresult=${2}
nscans=${3}

set +x

exe="/afs/hephy.at/work/t/tmadlener/chic_pol_global_fit/bin/run_pTM_scan"

outdir=$(dirname ${scanfile})

run_sandboxed ${outdir} ${exe} --outfile $(basename ${scanfile}) --npoints 1 --ptmmin 5 --noparams true --fitresult ${fitresult} --nscans ${nscans}
