# Global fit of j/psi, psi(2S) and chic cross sections and polarizations

This repository holds the code and the input data to run a global fit on j/psi,
psi(2S) and chic cross section and polarization data.

## Getting all inputs for the paper
The script [`scripts/run_all_options.sh`](scripts/run_all_options.sh) can be
used to run all three fit variant and produce majority of all plots and tables
that can be found in the paper. The only plots that are missing are the 2d
correlation plots, since the intermediate scan files they produce take up a lot
of disk space and since they can be easily parallelized these have been run
using a batch system. Additionally all the plots will only contain one variant,
but all the inputs to produce combination plots are present.

### Customization options for `run_all_options.sh`
- `RUN_FIT`: Either `0` or `1`. If `1` the global fit will be run for each
  variant. If `0` it assumes that this has already happened and the
  corresponding output file is present in the output directory.
- `RUN_SCAN`: Either `0` or `1`. If `1` the post fit scan of the parameter
  distribution to obtain the 1d bands as a function of `pT/M` and the input
  files to derive the values at `pT/M=5` that will be put into the summary
  table.
- `OUT_BASE`: The postfix of the output directory. It will be placed in the
  current working directory under `result_${OUT_BASE}`.

Note that running this script will automatically recompile all the c++
executables for each of the three fit scenarios. Additionally it will do a `make
clean` at the end to avoid having remnants of a previous run still linger
around.

## Building the c++ executables
Simply running `make` should build all the executables necessary to run the fit.
Most imporantly it builds `run_fit` in the `bin` folder. Possible targets for
`make` are:

- `all`: build all executables
- `debug`: build all executables with `-O0` and debug symbols
- `bin/[executable]`: To specifically build one of the executables defined by
  `.cc` files in the `src` folder
- `clean`: Get a clean build environment (equivalent to removing the `bin`
  folder and all its contents)
  
### Building for the different fit variants
It is possible to switch the parametrization of the cross section models at
compile time using the `FIT_OPTION` pre-processor variable. From a clean build
environment it is possible to specify this for make using the `ADD_FLAGS`
argument like this: `make [target] ADD_FLAGS=-DFIT_OPTION=X`, where `X` has to
be replaced with one of the following options:

- `0`: The `6 beta` version, where the psi and chic states each get two
  independent beta parameters (one for the total, one for the longitudional
  cross section)
- `1`: The total cross sections of all states are fitted to the same value and
  only the longitudinal beta are allowed to be different. This is the `4 beta`
  version
- `2`: All betas are replaced by one beta (i.e. `1 beta`). This implies that the
  total cross section of all states has the same pT/M behavior, and that all
  polarizations and feed-down fractions are pT/M independent.

Using a different value will make the build fail.

## Executables
The following c++ executables are present in the repository. Each `.cc` file in
the `src` folder will get an executable in the `bin` folder of the same name
(minus the `.cc`).

### [`run_fit`](src/run_fit.cc)
Runs the global fit and haves the following options. When scanning the inputs
are the lambdas of the chic1 and the chic2, even though internally the fit is
done in terms of the longitudinal fractions of the two states. For the grid scan
the corresponding values are determined such that the grid is regular in the
lambdas instead of the fractions.

- `--outfile` (`std::string`: `"results/scan_file.root"`): The file into which
  the fit results and the `TTree` containing the scanning values are stored.
- `--graphfile` (`std::string`: `"results/fit_result_graphs.root"`): The file
  into which the corrected data graphs and the best fit curves are stored.
- `--nographs` (`bool`: `false`): Do not produce the file containing the data
  graphs and best fit models.
- `--datadir` (`std::string`: `"./data/"`): The path to the input data files
  (e.g. when not running in  the root directory of the repository).
- `--useCosthRatios` (`bool`: `true`): Use the chic2 / chic1 ratios as a
  function of costh in the global fit.
- `--usePsiPolarizations` (`bool`: `true`): Use the j/psi and psi(2S)
  polarization data in the global fit.
- `--fixLambdasChicBestFit` (`bool`: `false`): Run two fits. The first to get
  the best fit values for the longitudinal fractions of the chic1 and chic2 and
  a second one with these values fixed. Only the results of the second fit are
  stored in the output file.
- `--useUnphysicalCorrections` (`bool`: `false`): Use corrections according to
  the lambdas of the chic1 and chic2 even if these values are outside of the
  physically allowed range, where the corrections are not necessarily
  applicable. In the default case such corrections are not used, but instead the
  corresponding extreme corrections are used.
- `--nscan1` (`unsigned`: `51`): Number of scan points for lambda(chic1)
- `--nscan2` (`unsigned`: `51`): Number of scan points for lambda(chic2)
- `--flow1` (`double`: `-1`): Minimum value of lambda(chic1)
- `--fhigh1` (`double`: `1`): Maximum value of lambda(chic1)
- `--flow2` (`double`: `-1`): Minimum value of lambda(chic2)
- `--fhigh2` (`double`: `1`): Maximum value of lambda(chic2)
- `--randomScan` (`bool`: `false`): Instead of doing a scan on a grid of the
  chic longitudinal fractions, do a random scan, where the parameter values are
  randomly drawn from the multivariate normal distribution described by the best
  fit parameter values and the covariance matrix at the minimum.
- `--npoints` (`unsigned`: `1000000`): The number of random scanning points.
- `--noscan` (`bool`: `false`): Do not do a chi2 scan.

### [`run_pTM_scan`](src/run_pTM_scan)
Does random parameter extractions and evaluates the fit models using these
parameter values. This in turn can be used to determine uncertainty bands for
different variables as a function of pT/M or to extract multidimensional
parameter correlations. The main assumption in this case is that the global
likelihood has parabolic behavior around its minimum as as the random parameter
values are extracted from a multivariate normal distribution centered at the
best fit values and using the covariance matrix around the minimum.

- `--fitresult` (`std::string`): File containing the best fit parameter values
  and the covariance matrix at the minimum. (Output from `run_fit`)
- `--outfile` (`std::string`: `"ptm_dep_scan_results.root"`): Output file into
  which the scanning points are stored
- `--ptmmin` (`double`: `0`): Minimum pT/M value (inclusive)
- `--ptmmax` (`double`: `10.0`): Maximum pT/M value (inclusive)
- `--npoints` (`size_t`: `50`): Number of pT/M scanning points
- `--nscans` (`size_t`: `20000`): Number of scans at each pT/M scan points
- `--noparams` (`bool`: `false`): Also store the parameter values for each scan
  point instead of just storing the computed values of the fit model.

## Python scripts
Some python scripts are located in [`python`](python). Most of these scripts
rely on the utilities that are defined in
[chib_chic_polFW](https://github.com/tmadlener/chib_chic_polFW/tree/master/python)
and the fact that `setup.sh` has been called for that repository in order for
the `PYTHONPATH` to also point to the corresponding modules.

The most important scripts are:
- [`make_fit_result_plots.py`](python/make_fit_result_plots.py) for plotting the
  fit results in comparison with the data graphs.
- [`ptm_dep_band_graphs.py`](python/ptm_dep_band_graphs.py) for turning the output
  files of `run_pTM_scan` into `TGraphAsymmErrors` as a function of `pT/M` for
  easier plotting afterwards.
- [`make_1d_band_plots.py`](python/make_1d_band_plots.py) for plotting the 1d bands
  produced by `ptm_dep_band_graphs.py`
- [`correlation_graphs.py`](python/correlation_graphs.py) for producing the 2d
  correlation graphs from the results of a high statistics run of `run_pTM_scan`
  (Of course it will also work for a regular run, but it will simply lump
  together all `pT/M` bins and pretend that only one value is present)
- [`make_2d_corr_plots.py`](python/make_2d_corr_plots.py) for plotting the 2d correlation plots.
- [`make_fit_par_table.py`](python/make_fit_par_table.py) for producing the summary
  tables. Note that this works best (and is also done automatically) when used
  in conjunction with `run_all_options.sh`

Other scripts are mainly there for historical reasons from when everything was
still under development. They have been left in place in case they could become
interesting again in the future. They are all related to the plan of obtaining
lambda(chic2) vs lambda(chic1) contours directly from scanning the chi2 space.
All of these work, but are much more restricted in generality compared to the
others. In the end we decided to not use this method to make the paper more
consistent in the treatment of derived quantities.

These scripts are:
- [`contour.py`](python/contour.py) for obtaining contours at arbitrary
  confidence levels from a chi2 scan.
- [`plot_contours.py`](python/plot_contours.py) for plotting the contours
- [`plot_scan_debug.py`](python/plot_scan_debug.py) for producing some debugging
  plots for a chi2 scan
- [`run_lambda_scan.py`](python/run_lambda_scan.py) for launching chi2 scan jobs
  onto the `slurm` batch system

### Other scripts
Some bash scripts are located in [`scripts`](scripts). The most important one
has already been described [above](#getting-all-inputs-for-the-paper). The
others are briefly described here. They are not as polished as the python
scripts and probably need some sort of adaption before they can be used.
 
- [`batch_correlation_scan.sh`](scripts/batch_correlation_scan.sh) Batch script
  for `slurm` that runs the scan for obtaining the 2d correlation plots. Note
  that the path to the built `run_pTM_scan` is hardcoded and needs to be adapted
  in any case.
- [`batch_run_fit.sh`](scripts/batch_run_fit.sh) Base batch script file for
  running the fit on the `slurm` batch system. It only takes the output file as
  an argument and simply passes all other arguments on directly to `run_fit`.
  This is not intended to be used as standalone script but rather in conjunction
  with `run_lambda_scan.py` which manages the
  submission of batch jobs for a specified grid on lambda(chic2) and
  lambda(chic1) on a higher level and is much easier to use. Note that also here
  the path to the `run_fit` executable is hardcoded and will need adaption
  before it can be used.
