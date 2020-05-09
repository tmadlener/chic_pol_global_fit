#!/bin/bash
# Script to run the fit with the three different options and producing the files
# containing the bands vs. pT/M
set -e

RUN_FIT=1
RUN_SCAN=1

OUT_BASE="w_ATLAS_dimuon_no_low_pT"

declare -A FIT_OPTIONS
FIT_OPTIONS["equal_total_cs"]=1
FIT_OPTIONS["equal_betas"]=2
FIT_OPTIONS["free_betas"]=0

FIT_EXE=bin/run_fit
FIT_EXE_ARGS="--noscan true"

SCAN_EXE=bin/run_pTM_scan
SCAN_EXE_ARGS="--ptmmin 2 --ptmmax 10 --npoints 50 --nscans 50000 --noparams true"

PLOT_FR_SCRIPT=python/make_fit_result_plots.py
BAND_SCRIPT=python/ptm_dep_band_graphs.py
PLOT_1D_SCRIPT=python/make_1d_band_plots.py
PLOT_2D_SCRIPT=python/make_2d_corr_plots.py

## Variables for bands
VARIABLES=(lth_jpsi lth_chic1 lth_chic2 lth_psip lth_jpsi_chic dlth_jpsi_psip\
                    r_chic1_jpsi r_chic2_jpsi r_chic_jpsi r_psip_jpsi\
                    dlth_chic2_chic1 chic2_chic1_cs_br)

BAND_VARS=""
for var in ${VARIABLES[@]}; do
    BAND_VARS=${BAND_VARS},${var}
done
# remove leading comma
BAND_VARS=$(echo ${BAND_VARS} | sed 's/\,//')


for fit_opt in ${!FIT_OPTIONS[@]}; do
    echo "============================================================"
    echo "Running: "${fit_opt}
    echo "============================================================"

    if [ ${RUN_FIT} -eq 1 -o ${RUN_SCAN} -eq 1 ]; then
        # make the executable with the correct settings
        make clean
        fit_opt_val=${FIT_OPTIONS["${fit_opt}"]}
        make all -j4 ADD_FLAGS=-DFIT_OPTION=${fit_opt_val}
    fi

    fitfile=results_${OUT_BASE}/fit_results_${fit_opt}.root
    graphfile=results_${OUT_BASE}/graphs_and_models_${fit_opt}.root
    plotdir=plots_${OUT_BASE}_${fit_opt}

    mkdir -p $(dirname ${fitfile})

    if [ ${RUN_FIT} -eq 1 ]; then
        ${FIT_EXE} ${FIT_EXE_ARGS} --outfile ${fitfile} --graphfile ${graphfile}
    fi
    python ${PLOT_FR_SCRIPT} --outdir ${plotdir}/fit_plots/ ${graphfile}

    scanfile=results_${OUT_BASE}/ptm_scan_ptm_2_10_50_50k_${fit_opt}.root
    bandfile=results_${OUT_BASE}/bands_1d_ptm_2_10_50_50k_${fit_opt}.root
    scanfile_ptm5=results_${OUT_BASE}/scan_ptm_5_2500k_${fit_opt}.root
    bandfile_pl=$(echo ${bandfile} | sed 's/.root/_phys_lambdas.root/')

    if [ ${RUN_SCAN} -eq 1 ]; then
        ${SCAN_EXE} ${SCAN_EXE_ARGS} --fitresult ${fitfile} --outfile ${scanfile}
        python ${BAND_SCRIPT} --variables ${BAND_VARS} --outfile ${bandfile} ${scanfile}
        python ${BAND_SCRIPT} --variables ${BAND_VARS} --physical-lambdas\
               --outfile ${bandfile_pl}  ${scanfile}

        ${SCAN_EXE} --ptmmin 5 --npoints 1 --nscans 2500000 --noparams true \
                    --fitresult ${fitfile} --outfile ${scanfile_ptm5}
    fi

    python ${PLOT_1D_SCRIPT} --outdir ${plotdir}/bands_v_ptm/ ${bandfile}
    python ${PLOT_1D_SCRIPT} --outdir ${plotdir}/phys_lambdas/bands_v_ptm/ ${bandfile_pl}

    # python ${PLOT_2D_SCRIPT} --outdir ${plotdir}/corr_2d/ \
    #        --input results/correlation_graphs_2d_ptm_5_${fit_opt}.root
    # python ${PLOT_2D_SCRIPT} --outdir ${plotdir}/phys_lambdas/corr_2d/ \
    #        --input results/correlation_graphs_2d_ptm_5_${fit_opt}_phys_lambdas.root

done

mkdir -p results_${OUT_BASE}/tables
python python/make_fit_par_table.py results_${OUT_BASE} --outfile results_${OUT_BASE}/tables/fit_results.tex

# just in case set back to clean state and force recompilation with the desired
# fit option
make clean
