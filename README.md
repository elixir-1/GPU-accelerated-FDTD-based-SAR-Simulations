# GPU-accelerated-FDTD-based-SAR-Simulations

Compilation Guide
1. For serial implementation, use intel compilers and compile it using 
  icc -mcmodel="large" serial_fdtd.c
2. For parallel implementation, compile using
  nvcc -arch=compute_35 -code=sm_35 parallelfdtd.cu
