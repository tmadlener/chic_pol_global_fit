# Global fit of j/psi, psi(2S) and chic cross sections and polarizations

This repository holds the code and the input data to run a global fit on j/psi,
psi(2S) and chic cross section and polarization data.

## Building
Simply running `make` should build all the executables necessary to run the fit.
Most imporantly it builds `run_fit` in the `bin` folder. Possible targets for
`make` are:

- `all`: build all executables
- `debug`: build all executables with `-O0` and debug symbols
- `bin/[executable]`: To specifically build one of the executables defined by
  `.cc` files in the `src` folder
- `clean`: Get a clean build environment (equivalent to removing the `bin`
  folder and all its contents)
  
### Building for pT-independent polarization
It is possible to switch the parametrization of the cross section models at
compile time using the `FIT_OPTION` pre-processor variable. From a clean build
environment it is possible to specify this for make using the `ADD_FLAGS`
argument like this: `make [target] ADD_FLAGS=-DFIT_OPTION=X`, where `X` has to
be replaced with one of the following options:

- `0`: The default version, where the psi and chic states each get two
  independent beta parameters (one for the total, one for the longitudional
  cross section)
- `1`: The total cross sections of all states are fitted to the same value and
  only the longitudinal beta are allowed to be different. 
- `2`: All betas are replaced by one beta. This implies that the total cross
  section of all states has the same pT/M behavior, and that all polarizations
  and feed-down fractions are pT/M independent.

Using a different value will make the build fail.

## Executables
The following c++ executables are present in the repository. Each `.cc` file in
the `src` folder will get an executable in the `bin` folder of the same name
(minus the `.cc`).

### `run_fit`
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

### `run_pTM_scan`
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
