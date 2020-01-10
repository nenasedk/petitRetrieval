# petitRadTrans Atmospheric Retrieval

This package is based around the examplar emission retrieval from the petitRadTrans paper [Molliere 2019, https://arxiv.org/abs/1904.11504]. 

https://petitradtrans.readthedocs.io/en/latest/content/notebooks/ret_emission_master.html

It has been modified to be more modular and to use a nested sampling approach for sampling the parameter space, 
as opposed to the MCMC approach of the original work. Credit is also due to the HELIOS-t atmospheric retrieval 
code, from which inspiration was drawn.

## Requirements
* petitRadTrans
* pyMultiNest, including building the fortran Multinest libraries
* MPI and mpi4py for parallelization

## Useage
Configuration of the retrieval is done through the config.py file. This file contains paths to data inputs and outputs,
parameter selection and hyperparameters to run the nested sampling code. The input data files should be formatted in a 
space separated .dat file, with columns
```
wavelength [micron], F_pl/F_s, F_pl/F_s [error]
```
An option for calibrated flux files rather than contrast are under development.
Additional consideration should be given to priors.py, which contains the prior distributions for the parameters used, 
as well as the log-likelihood function.

This code should be run in parallel over a large number of cores. Using MPI, run:
```
mpiexec -n $NUMTHREADS python retrieve_emission.py
```

Author: Evert Nasedkin

contact: nasedkinevert@gmail.com
