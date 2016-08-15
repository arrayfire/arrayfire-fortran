AF_PATH?=/opt/arrayfire
AF_LIB_NAME?=af
AF_FORT=af_fortran

LIB:=lib
AF_CFLAGS  = -I$(AF_PATH)/include

AF_FORT_LIB_PATH = $(AF_FORT_PATH)/$(LIB)/
AF_LIB_PATH      = $(AF_PATH)/$(LIB)/

AF_FORT_LIB=$(AF_FORT_LIB_PATH)/lib$(AF_FORT).so
AF_FORT_FILE=$(AF_FORT_LIB_PATH)/arrayfire.f90
