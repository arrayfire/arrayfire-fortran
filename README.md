arrayfire-fortran
=================

This project provides Fotran bindings for ArrayFire.

Prerequisites
---------------

- The latest version of ArrayFire. You can get ArrayFire in the following ways.
    - [Binary Installer](http://www.arrayfire.com/download)
    - [Install from source](http://github.com/arrayfire/arrayfire)

- `gfortran`

- `make`

Contents
---------------

- `src/`: Contains the source files for the ArrayFire Fortran wrapper
    - `fortran_wrapper.cpp` The C++ part of the wrapper
    - `arrayfire.f95` The fortran part of the wrapper

- `lib/` The location where the wrapper library is stored

- `examples`: contains a few examples demonstrating the usage


Usage
----------------

After you the necessary pre-requisites, edit the following paramets in `common.mk`

- Change `AF_PATH` to the right location
- Change `AF_LIB_NAME` to point to the right backend.


### Linux

- To build the Fortran Wrapper for ArrayFire run
    - `make all`(generates `libaf_fortran.so`)

- To build the examples do one of the following from the examples directory
    - `make all`   (generates `examplename`)

Documentation
---------------

- Work under progress

License
---------------

This project is licensed under BSD 3 clause license.

Please check the LICENSE file in the root directory for more information.
