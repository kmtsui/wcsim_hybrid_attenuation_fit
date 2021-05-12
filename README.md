# wcsim_hybrid_attenuation_fit

Simple analysis to extract water attenuation length from wcsim_hybrid diffuser MC

Setup your ROOT and $WCSIMDIR before compilation

From within the cloned repository

```
$ mkdir build; cd build; cmake ../
$ make install
```

After the build you can setup your environment

```
$ source Linux/setup.sh
```
  
The analysis is done in two steps. First use analysis_absorption to perform data reduction.

    $ analysis_absorption -f wcsim_output.root 

Then use the root macro fit_water_attenuation.c to do the fit

    $ root fit_water_attenuation.c
