AF_FORT_PATH?=$(shell (dirname `pwd`))

-include $(AF_FORT_PATH)/common.mk

AF_FORT_MOD = $(AF_FORT_PATH)/src/arrayfire.f95
LDFLAGS += -Wl,--no-as-needed -L$(AF_FORT_LIB_PATH) -l$(AF_FORT) -L$(AF_LIB_PATH) -l$(AF)
LDFLAGS += -Wl,-rpath,$(AF_FORT_LIB_PATH),-rpath,$(abspath $(AF_FORT_LIB_PATH))
LDFLAGS += -Wl,-rpath,$(AF_LIB_PATH),-rpath,$(abspath $(AF_LIB_PATH))

SRC:=$(wildcard *.f95)
BIN:=$(patsubst %.f95, %_$(EXT), $(SRC))

cuda: $(BIN)

opencl: $(BIN)

%_ocl: %.f95 $(AF_FORT_MOD)
	gfortran -L$(AF_PATH)/$(LIB) $(CFLAGS) $(LDFLAGS) $(AF_FORT_MOD) -o $@ $<

%_cuda: %.f95 $(AF_FORT_MOD)
	gfortran -L$(AF_PATH)/$(LIB) $(CFLAGS) $(LDFLAGS) $(AF_FORT_MOD) -o $@ $<

clean:
	rm -f *_ocl *_cuda *.mod
