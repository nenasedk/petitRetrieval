# petitRadTrans Atmospheric Retrieval

This package is based around the examplar emission retrieval from the petitRadTrans paper [Molliere 2019, https://arxiv.org/abs/1904.11504]. 

It has been modified to be more modular and to use a nested sampling approach for sampling the parameter space, 
as opposed to the MCMC approach of the original work. Credit is also due to the HELIOS-t atmospheric retrieval 
code, from which inspiration was drawn.

## Requirements
* [petitRadTrans](https://petitradtrans.readthedocs.io/en/latest/content/notebooks/ret_emission_master.html)
* [pyMultiNest](https://github.com/JohannesBuchner/PyMultiNest), including building the [Fortran Multinest libraries](https://github.com/JohannesBuchner/MultiNest). Recommend building pymultinest from source.
* MPI and mpi4py for parallelization

## Usage
Configuration of the retrieval is done through the config.py file. 
This file contains paths to data inputs and outputs,parameter selection and hyperparameters to run the nested sampling code. 
The input data files should be formatted in a space separated .dat file, with columns
```
wavelength [micron], F_pl/F_s, F_pl/F_s [error]
```
if using a contrast measurement or
```
wavelength [micron], F_pl [uJy], F_pl_err [uJy]
```
if using a photometrically calibrated flux measurement.
All prior setup is done in a dictionary in the config file, where priors and boundaries on planet properties can be set.
The user must ensure the parameters in the ATMOSPHERE list match those included in the prior dictionary. 
While the number of Live Points is a user set parameter, it is recommended to use on the order of 50 times the number of free parameters.

Additional consideration should be given to priors.py, which contains the prior distributions for the parameters used, 
as well as the log-likelihood function.

This code should be run in parallel over a large number of cores. Using MPI, run:
```
mpiexec -n $NUMTHREADS python retrieve_emission.py
```

Plotting is best done using the multinest_marginals_corner.py tool.
To use this, navigate to the directory of the output files from the retrieval,
and run
```
/PATH/TO/PETITRETRIEVAL/multinest_marginals_corner.py PREFEX
```
Where prefix is the prefix prepended to all of the pymultinest output files.
This will be automated in a future release.

Author: Evert Nasedkin
Contact: nasedkinevert@gmail.com
