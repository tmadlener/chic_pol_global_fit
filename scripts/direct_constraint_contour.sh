#!/bin/bash
# Script to produce the lambda_2 vs lambda_1 contour from the costh ratio measurments

# Make the script callable from everywhere by getting the directory of the
# script and then defining the main directory of the framework based on that
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
MAINDIR=${SCRIPTDIR}/../


OUTDIR=${MAINDIR}/results/
mkdir -p ${OUTDIR}

SCANFILE=${OUTDIR}/costh_chi2_scan.root
CONTOURFILE=${OUTDIR}/contours_direct_constraints.root


# Run the costh fit. Nothing to do here, everything is fixed
${MAINDIR}/bin/run_costh_fit --outfile ${SCANFILE}

python ${MAINDIR}/python/contour.py --variable-x "lambda_1" --variable-y "lambda_2" \
       --outfile ${CONTOURFILE} --conf-levels 0.683,0.955,0.997 \
       ${SCANFILE}
