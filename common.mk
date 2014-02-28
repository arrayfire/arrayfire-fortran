AF_PATH ?= /opt/arrayfire
CUDA_PATH ?= /usr/local/cuda/
OCL_PATH ?= $(AF_PATH)

ifneq ($(shell uname), Linux)
$(error Only Linux supported for fortran)
endif

ifeq ($(shell uname -m), x86_64)
  LIB:=lib64
else
  LIB:=lib
endif

AF_CFLAGS  = -I$(AF_PATH)/include
ifeq ($(findstring opencl, $(MAKECMDGOALS)), opencl)
	AF_CFLAGS += -DAFCL -I$(OCL_PATH)/include
	AF_FORT=afcl_fortran
	AF=afcl
	EXT=ocl
else
	AF_CFLAGS += -I$(CUDA_PATH)/include
	AF_FORT=afcu_fortran
	AF=afcu
	EXT=cuda
endif

AF_FORT_LIB_PATH = $(AF_FORT_PATH)/$(LIB)/
AF_LIB_PATH      = $(AF_PATH)/$(LIB)/
