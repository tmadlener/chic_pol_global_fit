#!/bin/bash
# Script to produce the lambda_2 vs lambda_1 contour from the cross section and
# psi polarization measurements adding also the costh ratios measurements

# Make the script callable from everywhere by getting the directory of the
# script and then defining the main directory of the framework based on that
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
MAINDIR=${SCRIPTDIR}/../

low=$1
high=$2
steps=$3

OUTDIR=${MAINDIR}/results/
mkdir -p ${OUTDIR}/contours/

SCANFILE=${OUTDIR}/scan_regular_lambda_all_constraints_lambda2_${low}_${high}_${steps}.root
CONTOURFILE=${OUTDIR}/contours/contours_all_constraints.root

${MAINDIR}/bin/run_fit --outfile ${SCANFILE} \
              --flow1 -1.5 --fhigh1 1.5 --nscan1 301 \
              --flow2 ${low} --fhigh2 ${high} --nscan2 ${steps} \


# python ${MAINDIR}/python/contour.py --variable-x "f_long_c1" --transform-x "frac_to_lam" \
#        --variable-y "f_long_c2" --transform-y "frac_to_lam" --conf-levels 0.683,0.955,0.973 \
#        --outfile ${CONTOURFILE} ${SCANFILE}
