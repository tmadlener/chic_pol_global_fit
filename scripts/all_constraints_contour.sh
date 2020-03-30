#!/bin/bash
# Script to produce the lambda_2 vs lambda_1 contour from the cross section and
# psi polarization measurements adding also the costh ratios measurements

# Make the script callable from everywhere by getting the directory of the
# script and then defining the main directory of the framework based on that
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
MAINDIR=${SCRIPTDIR}/../


OUTDIR=${MAINDIR}/results/
mkdir -p ${OUTDIR}/contours/

SCANFILE=${OUTDIR}/scan_regular_lambda_all_constraints.root
CONTOURFILE=${OUTDIR}/contours/contours_all_constraints.root

${MAINDIR}/bin/run_fit --outfile ${SCANFILE} \
              --flow1 -1.5 --fhigh1 1.5 --nscan1 16 \
              --flow2 -1.65 --fhigh2 1.65 --nscan2 16 \


python ${MAINDIR}/python/contour.py --variable-x "f_long_c1" --transform-x "frac_to_lam" \
       --variable-y "f_long_c2" --transform-y "frac_to_lam" --conf-levels 0.683,0.955,0.973 \
       --outfile ${CONTOURFILE} ${SCANFILE}
