# wcsim_hybrid_attenuation_fit

Analysis framework to extract detector parameters from WCSIM hybridPMT MC.

## Installation

Requirements:
- c++11 compiler
- cmake 2.8+
- ROOT 5.34.34+ or 6+
- (Optional) WCSIM hybridPMT branch

Setup your ROOT before compilation. It is required to have Minuit2 in ROOT for minimization. OpenMP installation is recommended to take advantage of parallelism when doing the fit.

From within the cloned repository

```
$ mkdir build; cd build; cmake ../
$ make install
```

To build with WCSIM support, setup `$WCSIMDIR` which contains the WCSIM `src/` and `include/` directory, and `libWCSimRoot.so`.

```
$ cmake ../ -DUSE_WCSIM=1
```

After the build you can setup your environment

```
$ source Linux/setup.sh
```
  
## WCSIM_TreeConvert

The analysis is done in two steps. First use `WCSIM_TreeConvert` to perform data reduction on WCSIM output.

    $ WCSIM_TreeConvert -f wcsim_output.root 

The program assumes a diffuser simulation and store the basic PMT hits and PMT geometry (relative to the source) information in `TTree` format. Modify the code if you want to store extra information.

In WCSIM hybridPMT branch, there are the mPMTs in addition to the ordinary PMTs. The two types of PMT hits and geometries are stored in separate trees.

## optical_fit

The fitter `optical_fit` imports the data samples produced by `WCSIM_TreeConvert` and fits for a set of detector parameters defined by users.

    $ optical_fit -o fitoutput.root -c config.toml

The sample and fit configurations are defined in the toml file. See `$OPTICALFIT/var/OPTICALFIT/config/config.toml` for detailed explanation. You can copy the entire `config` folder to somewhere else and run the code there.

`AnaSample` class handles the input of analysis samples to load the PMT hits and geometry. Modify `AnaEvent` and `AnaTree` classes to load extra information if necessary. `AnaSample` bins the PMT hits (p.e.) according to the binning defined in `config.toml`.

`AnaFitParameters` class defines the fit parameters which parameterize the number of hits expected in each PMT. The expected numbers are compared with the observation in `AnaSample` to compute a chi2.

`Fitter` class does the chi2 minimization with respect to the fit parameters and saves the output in `fitoutput.root`.

The fitter is adapted from T2K xsllhFitter at https://gitlab.com/cuddandr/xsLLhFitter
