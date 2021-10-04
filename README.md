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

To build with WCSIM support, setup `$WCSIMDIR` which contains the WCSIM `src/` and `include/` directory, and `libWCSimRoot.so`.

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
See available arguments with `WCSIM_TreeConvert -h`.

The program assumes a diffuser simulation and store the basic PMT hits and PMT geometry (relative to the source) information in `TTree` format. Modify the code if you want to store extra information.

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

The `hitRate_pmtTypeX` tree contains the branches:
- `nPE/D` : number of PE
- `timetof/D` : time-to-flighted subtracted hit time
- `PMT_id/I` : PMT id that the hit belongs to

## optical_fit

The fitter `optical_fit` imports the data samples produced by `WCSIM_TreeConvert` and fits for a set of detector parameters defined by users.
```
optical_fit -o fitoutput.root -c $OPTICALFIT/var/OPTICALFIT/config/config.toml
```
See available arguments with `optical_fit -h`.

The sample and fit configurations are defined in the toml file. See `$OPTICALFIT/var/OPTICALFIT/config/config.toml` for detailed explanation. You can copy the entire `config` folder to somewhere else and run the code there. The config and binning files are assumed to live in the same folder.

`AnaSample` class handles the input of analysis samples to load the PMT hits and geometry. Modify `AnaEvent` and `AnaTree` classes to load extra information if necessary. `AnaSample` bins the PMT according to the binning file defined in `config.toml`, later which a binned Possion likelihood function is used to calculate a chi2.

`AnaFitParameters` class defines the fit parameters which parameterize the number of hits expected in each PMT. The expected numbers are compared with the observation in `AnaSample` to compute the chi2.

`Fitter` class does the chi2 minimization with respect to the fit parameters and saves the output in `fitoutput.root`. The output includes the event histogram `evhist_sam*_[data/pred]`, the chi2 evolution `chi2_*_periter`, the post-fit parameter values `res_vector`, covariance matrix `res_cov_matrix` and correlation matrix `res_cor_matrix`. The output also contains a `PMTTree` which includes all the PMTs involved in the fit, their observed and predicted PE after the fit.

The fitter is adapted from T2K xsllhFitter at https://gitlab.com/cuddandr/xsLLhFitter

## Container
Container image is available for docker and singularity. See `container/README.md` for instructions.