AF_FORT_PATH=$(shell pwd)

-include $(AF_FORT_PATH)/common.mk

all: $(AF_FORT_LIB)

$(AF_FORT_LIB): $(AF_FORT_PATH)/src/fortran_wrapper.cpp
	gfortran -shared -fPIC $< $(AF_CFLAGS) -L$(AF_LIB_PATH) -l$(AF_LIB_NAME) -o $@

clean:
	rm -f $(AF_FORT_LIB)
	make -C $(AF_FORT_PATH)/examples clean
