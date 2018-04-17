# GPU_functions

## Description

These are some cuda codes based on the cuda N-body examples distributed long time ago. The code is a straighforward brute force and it is only fast because it runs on the GPU.

Included are IDL and c codes that perform the same task for comparison purposes.

The IDL functions show how to call the c/cuda code from IDL.

### PotGPU

Computes gravitational potential for a particle distribution. Useful to perform remove unbound particles in halos/subhalos.

### GPU_distances

This is basically the same code in PotGPU but changed to compute distance transforms. The code was used for this paper: 

http://adsabs.harvard.edu/abs/2015MNRAS.454..463A

