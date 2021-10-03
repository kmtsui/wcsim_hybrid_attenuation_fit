## singularity usage
In case you don't have the proper dependencies or don't want to install them at all, you can now use the container image hosted on https://cloud.sylabs.io/. The image is built with the recipe `optical_fit.def`.
```
singularity pull library://kmtsui/optical_fit/optical_fit_wcsim:latest
```
After the pull, you can enter the container and use the app as usual.
```
singularity shell optical_fit_wcsim.sif
WCSIM_TreeConvert -h
optical_fit -h
```
The `optical_fit` and `WCSim` packages are `git clone` to `/wcsim_hybrid_attenuation_fit` and `/WCSim`, and they are built inside the respective `build` directory.

If you want to modify and recompile any package, you can build the container with `--sandbox` option and `shell --writable` it.
```
singularity build --sandbox optical_fit library://kmtsui/optical_fit/optical_fit_wcsim:latest
singularity shell --writable optical_fit
```
If there is any mounting problem, try the `--no-home` option. You will probably lose access to the local file system, so you might need to mount the file system yourself, or exit and re-enter the container with the usual `shell` command afterwards.
```
singularity shell --writable --no-home optical_fit
```

## docker usage
Docker image is also available on docker hub. The difference is that the  `optical_fit` and `WCSim` packages are not `git clone` but installed from tar. Every time the `experimental` branch is updated, an image built is triggered and handled by https://github.com/kmtsui/optical_fit_docker/.
```
docker run -it --rm --name optical_fit_container kmtsui/opticalfit:latest bash
```