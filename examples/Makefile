AF_FORT_EXAMPLE_DIR=$(shell pwd)
AF_FORT_PATH?=$(shell (dirname $(AF_FORT_EXAMPLE_DIR)))
AF_FORT_BIN_DIR?=$(AF_FORT_PATH)/bin

-include $(AF_FORT_PATH)/common.mk

AF_FORT_MOD = $(AF_FORT_LIB_PATH)/arrayfire.f90
LDFLAGS += -Wl,--no-as-needed -L$(AF_FORT_LIB_PATH) -l$(AF_FORT) -L$(AF_LIB_PATH) -l$(AF_LIB_NAME)
LDFLAGS += -Wl,-rpath,$(AF_FORT_LIB_PATH),-rpath,$(abspath $(AF_FORT_LIB_PATH))
LDFLAGS += -Wl,-rpath,$(AF_LIB_PATH),-rpath,$(abspath $(AF_LIB_PATH))

SRC:=$(wildcard $(AF_FORT_EXAMPLE_DIR)/*.f90)
BIN:=$(patsubst $(AF_FORT_EXAMPLE_DIR)/%.f90, $(AF_FORT_BIN_DIR)/%, $(SRC))

all: $(BIN)

$(AF_FORT_BIN_DIR)/%: $(AF_FORT_EXAMPLE_DIR)/%.f90 $(AF_FORT_MOD)
	@echo Building $(shell (basename $@))
	@gfortran -L$(AF_PATH)/$(LIB) $(CFLAGS) $(LDFLAGS) $(AF_FORT_MOD) -o $@ $<

clean:
	rm -f $(BIN)
