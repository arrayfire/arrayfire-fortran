arrayfire_fortran
=================

This repository contains the files required to use ArrayFire from Fortran^.

Prerequisites
---------------

- The latest version of ArrayFire. You can [download here](http://www.accelereyes.com/download_arrayfire)
    - All the pre-requisites for ArrayFire still apply.

- `gfortran`

- `make`

- `Linux` Right now other distributions are not supported

Contents
---------------

- src/: Contains the source files for the ArrayFire Fortran wrapper
    - fortran_wrapper.cpp: The C++ part of the wrapper
    - arrayfire.f95: The fortran part of the wrapper

- lib/, lib64/: The location where the wrapper library is stored

- examples: contains a few examples demonstrating the usage


Usage
----------------

After you the necessary pre-requisites, edit the following paramets

- Open common.mk and change `AF_PATH` to the right location


### Linux

- To build the Fortran Wrapper for ArrayFire run
    - `make cuda all`   to use CUDA   (generates `libafcu_fortran.so`)
    - `make opencl all` to use OpenCL (generates `libafcl_fortran.so`)

- To build the examples do one of the following from the examples directory
    - `make cuda  ` to use CUDA   (generates `examplename_cuda`)
    - `make opencl` to use OpenCL (generates `examplename_ocl` )

Documentation
---------------

- The documentation can be found [over here](http://www.accelereyes.com/arrayfire/fortran/)

License
---------------

- Please check the LICENSE file in the root directory
