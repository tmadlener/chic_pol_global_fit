#!/bin/bash
#SBATCH -J GlobalFit
#SBATCH -D /afs/hephy.at/work/t/tmadlener/chic_pol_global_fit/
#SBATCH -o /afs/hephy.at/work/t/tmadlener/chic_pol_global_fit/batch_logfiles/run_global_fit_%A_%a.out

FITTER="/afs/hephy.at/work/t/tmadlener/chic_pol_global_fit/bin/run_fit"

source ${CHIB_CHIC_POLFW_DIR}/scripts/bash_functions.sh

set -x
## the scan file is necessary to determine the output directory for
## run_sandboxed
scanfile=${1}
shift

outdir=$(dirname ${scanfile})

run_sandboxed ${outdir} ${FITTER} --outfile $(basename ${scanfile}) $@
set +x
