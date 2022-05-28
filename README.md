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
mkdir build; cd build; cmake ../
make install
```

To build with WCSIM support, setup environment variable `$WCSIMDIR` which contains the WCSIM `src/` and `include/` directory, and `libWCSimRoot.so`.

```
cmake ../ -DUSE_WCSIM=1
make install
```

After the build you can setup your environment

```
source Linux/setup.sh
```
  
## WCSIM_TreeConvert

The analysis is done in two steps. First use `WCSIM_TreeConvert` to perform data reduction on WCSIM output.
```
WCSIM_TreeConvert -f wcsim_output.root 
```
Available arguments are
- `-f` : input file name
- `-o` : output file name
- `-l` : laser wavelength. For calculating speed of light in water
- `-w` : apply diffuser angular profile reweight
- `-z` : apply top-bottom-asymmetry reweight with input slope of attenuation length
- `-p` : set water parameters ABWFF,RAWFF for -z reweight
- `-b` : only read B&L PMTs data
- `-d` : read raw Cherenkov hits and perform ad-hoc digitization
- `-t` : use separated triggers for B&L and mPMTs
- `-v` : turn on detailed verbose
- `-s` : specify start event
- `-e` : specify end event
- `-r` : random seed value

The program assumes a light injector simulation with fixed source position and store the basic PMT hits and PMT geometry (relative to the source) information in `TTree` format. Modify the code if you want to store extra information.

In WCSIM hybridPMT branch, there are the mPMTs in addition to the ordinary PMTs. The two types of PMT hits and geometries are stored in separate trees.

In the output file, there are two types of trees: `pmt_typeX` is the PMT geometry tree, `hitRate_pmtTypeX` is the PMT hit tree, where `X=0,1` corresponds to the B&L PMT and mPMT respectively.

The `pmt_typeX` tree contains the branches:
- `R/D` : distance to source 
- `costh/D` : cosine of photon incident angle relative to PMT
- `cosths/D` : PMT costheta angle relative to source
- `phis/D` : PMT phi angle relative to source
- `omega/D`: solid angle subtended by PMT
- `PMT_id/I` : unique PMT id
- `mPMT_id/I` : ID of small PMT inside a mPMT module
- `costhm/D` : photon incident costheta angle relative to central small PMT for mPMT
- `phim/D` : photon incident phi angle relative to central small PMT for mPMT
- `dz/D` : z-position relative to source

The `hitRate_pmtTypeX` tree contains the branches:
- `nPE/D` : number of PE
- `timetof/D` : time-to-flight subtracted hit time
- `PMT_id/I` : PMT id that the hit belongs to

The `-d` option stores addition truth information of photon hit time/charge, and photon reflection/scattering history.

## optical_fit

The fitter `optical_fit` imports the data samples produced by `WCSIM_TreeConvert` and fits for a set of detector parameters defined by users.
```
optical_fit -o fitoutput.root -c $OPTICALFIT/var/OPTICALFIT/config/config.toml
```
Available arguments are
- `-o` : output file name
- `-c` : input config file name
- `-s` : random seed value
- `-n` : number of threads 
- `-t` : number of toy fits. Repeat the fits with tweaks in certain parameters

The sample and fit configurations are defined in the toml file. See `$OPTICALFIT/var/OPTICALFIT/config/config.toml` for detailed explanation. You can copy the entire `config` folder to somewhere else and run the code there. The config and binning files are assumed to live in the same folder.

`AnaSample` class handles the input of analysis samples to load the PMT hits and geometry. Modify `AnaEvent` and `AnaTree` classes to load extra information if necessary. `AnaSample` bins the PMT according to the binning file defined in `config.toml`, later which a binned Possion likelihood function is used to calculate a chi2.

`AnaFitParameters` class defines the fit parameters which parameterize the number of hits expected in each PMT. The expected numbers are compared with the observation in `AnaSample` to compute the chi2.

`Fitter` class does the chi2 minimization with respect to the fit parameters and saves the output in `fitoutput.root`. The output includes the event histogram `evhist_sam*_[data/pred]`, the chi2 evolution `chi2_*_periter`, the post-fit parameter values `res_vector`, covariance matrix `res_cov_matrix` and correlation matrix `res_cor_matrix`. The output also contains a `PMTTree` which includes all the PMTs involved in the fit, their observed and predicted PE after the fit.

The fitter is adapted from T2K xsllhFitter at https://gitlab.com/cuddandr/xsLLhFitter

### MCMC posterior
An optional MCMC algorithm is available to sample the likelihood surface around the best-fit point. It uses the Metropolis-Hastings Algorithm for acceptance-rejection, and by default an adaptive step proposal defined in `src/TSimpleMCMC.hh` (directly copied from https://github.com/ClarkMcGrew/root-simple-mcmc).

### GPU support
GPU threading is now available for CUDA (NVIDIA) architecture. To enable GPU support, source the cuda environment, then configure cmake with `-DUSE_CUDA=1`. It requires cmake 3.8+, and I only tested it on cuda V11.2.152.

Currently only the event reweighing in each MINUIT iteration is done in GPU. The memory copy between host (CPU) and device (GPU) is mostly done inside src/CacheManager, which is modified from the T2K GUNDAM fitter at https://github.com/nadrino/gundam. Other relevant codes are enclosed within `#ifdef USING_CUDA #endif`.

In the initialization stage, each event (PMT) weight, fit parameter and the necessary event information, and histogram binning is copied from host to device. The `CacheManager` class handles most of the logistics; the `Weight*` classes do the actual event reweighing; and the `CacheIndexedSums` class sums all the event weights into the correct histogram bins.

In each MINUIT iteration, the fit parameter values are copied from host to device; the event weights are calulated and summed into bins on device; and the results are copied back from device to host for chi2 calculation. If you define a new type of fit parameter for event reweighing or a new kind of histogram for chi2 calculation, those new components have to be defined inside CacheManager for proper GPU computation.

## Container
Container image is available for docker and singularity. See `container/README.md` for instructions.
