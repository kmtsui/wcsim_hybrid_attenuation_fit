# wcsim_hybrid_attenuation_fit

Simple analysis to extract detector parameters from wcsim_hybrid MC.

## Setup

Setup your ROOT and $WCSIMDIR before compilation. It is assumed that $WCSIMDIR contains the WCSIM src/ and include/ directory, and libWCSimRoot.so. It is required to have Minuit/Minuit2 in ROOT for minimization.

From within the cloned repository

```
$ mkdir build; cd build; cmake ../
$ make install
```

After the build you can setup your environment

```
$ source Linux/setup.sh
```
  
## analysis_absorption

The analysis is done in two steps. First use `analysis_absorption` to perform data reduction.

    $ analysis_absorption -f wcsim_output.root 

The program assumes a single source of photons in simulation and store the basic PMT hits and PMT geometry information in `TTree` format. Modify the code if you want to store extra information.

## optical_fit

The fitter `optical_fit` imports the data samples produced by `analysis_absorption` and fits for a set of detector parameters defined by users.

    $ optical_fit -o fitoutput.root -c config.toml

The sample and fit configurations are defined in the toml file `config.toml`. See `app/config/config.toml` for detailed explanation.

`AnaSample` class handles the input of analysis samples and loads the PMT hits and geometry. Modify `AnaEvent` and `AnaTree` classes to load extra information if necessary. `AnaSample` bins the PMT hits according to the binning defined in `config.toml`.

`AnaFitParameters` class defines the fit parameters which parameterize the number of PMT hits observed in the each sample bin.

`Fitter` class does the chi2 minimization and saves the output in `fitoutput.root`.
