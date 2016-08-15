AF_FORT_PATH?=$(shell (dirname `pwd`))

-include $(AF_FORT_PATH)/common.mk

AF_FORT_MOD = $(AF_FORT_LIB_PATH)/arrayfire.f90
LDFLAGS += -Wl,--no-as-needed -L$(AF_FORT_LIB_PATH) -l$(AF_FORT) -L$(AF_LIB_PATH) -l$(AF_LIB_NAME)
LDFLAGS += -Wl,-rpath,$(AF_FORT_LIB_PATH),-rpath,$(abspath $(AF_FORT_LIB_PATH))
LDFLAGS += -Wl,-rpath,$(AF_LIB_PATH),-rpath,$(abspath $(AF_LIB_PATH))

SRC:=$(wildcard *.f90)
BIN:=$(patsubst %.f90, %, $(SRC))

all: $(BIN)

%: %.f90 $(AF_FORT_MOD)
	gfortran -L$(AF_PATH)/$(LIB) $(CFLAGS) $(LDFLAGS) $(AF_FORT_MOD) -o $@ $<

clean:
	rm -f $(BIN)
