#!/usr/bin/env bash
set -e

# Some location things
START_DIR=$(pwd)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
BASE_DIR=$(realpath ${SCRIPT_DIR}/../)
PYTHON_DIR=${BASE_DIR}/python

# Some adjustable settings
PARAMS_SETTINGS=${BASE_DIR}/misc/settings_1S0_only_NLO.txt
OUTBASE="SDCs_Carlos_interpol/fit_1S0_only"
SDC_ORDER="NLO"

# Semi fixed settings, only change these if you know what you want
DATADIR=${BASE_DIR}/data
SDCDIR=$(realpath ${BASE_DIR}/../SDCs_Chung/)

# Fixed settings
FIT_EXE=${BASE_DIR}/bin/run_fit_nrqcd
SCAN_EXE=${BASE_DIR}/bin/run_pTM_scan_nrqcd
FIT_RESULT_FILE="fit_results_nrqcd_global_fit.root"
SCAN_RESULT_FILE="param_scan_ptm5.root"

# Computed filenames, directories and settings
result_dir=${BASE_DIR}/${OUTBASE}_${SDC_ORDER}
fit_result=${result_dir}/${FIT_RESULT_FILE}
scan_result=${result_dir}/${SCAN_RESULT_FILE}

# This is where the fun starts
set -x
mkdir -p ${result_dir}

# run the fit
${FIT_EXE} --datadir ${DATADIR} --sdcdir ${SDCDIR}\
    --order ${SDC_ORDER} \
    --outdir ${result_dir} \
    --paramsSettings ${PARAMS_SETTINGS} \

# run the scan at pT/M = 5
${SCAN_EXE} --ptmmin 5 --npoints 1 --nscans 100000 \
    --fitresult ${fit_result} --outfile ${scan_result} --sdcdir ${SDCDIR}

# Make a fit report
python ${PYTHON_DIR}/make_fit_result_report.py ${result_dir}
cd ${result_dir}
pdflatex fit_results_report.tex
rm fit_results_report.{aux,log}
cd ${START_DIR}

set +x
