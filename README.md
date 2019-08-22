# GPU-accelerated-FDTD-based-SAR-Simulations

Parallel Specific Absorption Rate using 3D Finite-Difference Time-Domain algorithm on CUDA framework.

Introduction
3D FDTD is computationally expensive due to its volumetric nature. Using CUDA, a GPU based implementation accelerates the algorithm. This faster implementation aids faster calculation of SAR for given frequency.

Requirements
CUDA
Thrust library
Intel compilers

src Modules
1. serial - serial C implementation of FDTD based SAR
2. parallel - FDTD based SAR accelerated using CUDA-C

Compilation Guide
1. For serial implementation, use intel compilers and compile it using 
  icc -mcmodel="large" serial_fdtd.c
2. For parallel implementation, compile using
  nvcc -arch=compute_35 -code=sm_35 parallelfdtd.cu
