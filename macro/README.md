## FitAngularResponsePol.c
Load post-fit angular response parameters from optical_fit output, then fit the histogram as piecewise continuous polynomials in a similar fashion as fiTQun

## makeLEDProfileInput.c
Produce source angular profile and covariance matrix as input to optical_fit

## make_scatter_map.c
From output from WCSIM_TreeConvert, calculate the number of indirect photons for each PMT in the signal and control regions, then save the ratio as TH1D which is used as input to `scatter_map` option of optical_fit

## truth_alpha.c 
Calculate the true light attenuation length as a function of photon wavelength and WCSIM water parameters

## PlotGui.c
A not so useful ROOT GUI for plotting post-fit parameters. Either supply the fitouput root file as function argument, or open the file with GUI `Open File` button. 

## build_template_and_spline.c
Build the indirect photon template for the each PMT, and genearate the spline of Latt and Lsct parameters
