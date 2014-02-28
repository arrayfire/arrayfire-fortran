AF_FORT_PATH=$(shell pwd)

-include common.mk

AF_FORT_LIB=$(AF_FORT_LIB_PATH)/lib$(AF_FORT).so

all: $(AF_FORT_LIB)

$(AF_FORT_LIB): $(AF_FORT_PATH)/src/fortran_wrapper.cpp
	gfortran -shared -fPIC $< $(AF_CFLAGS) -L$(AF_LIB_PATH) -l$(AF) -o $@

clean:
	rm -f $(AF_FORT_LIB)


.PHONY: opencl cuda
