AF_FORT_PATH=$(shell pwd)

-include $(AF_FORT_PATH)/common.mk

all: $(AF_FORT_LIB) $(AF_FORT_FILE)

$(AF_FORT_FILE): $(AF_FORT_PATH)/src/arrayfire.f90
	@echo Copying $(shell (basename $@))
	@cp $< $@

$(AF_FORT_LIB): $(AF_FORT_PATH)/src/fortran_wrapper.cpp
	@echo Building $(shell (basename $@))
	@gfortran -shared -fPIC $< $(AF_CFLAGS) -L$(AF_LIB_PATH) -l$(AF_LIB_NAME) -o $@

clean:
	rm -f $(AF_FORT_LIB)
	rm -f $(AF_FORT_FILE)
