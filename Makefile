.PHONY: all help info clean build_quadpack build_cephes build_amos build_numpy_mt19937 vp_benchmark_default vp_benchmark_rngstreams

CC                 := gcc
CFLAGS             := -O3 -m64 -mavx2 -Wall -std=c99 -fopenmp
LDFLAGS            := -lgfortran -lm

# Change the paths according to libary installation in your system
UNURAN             := -I./unuran-1.8.1/build/include ./unuran-1.8.1/build/lib/libunuran.a
RNGSTREAMS         := -I./rngstreams-1.0.1/build/include ./rngstreams-1.0.1/build/lib/librngstreams.a


OPENSPECFUN        := ./openspecfun/libopenspecfun.a
CEPHES             := ./cephes/libprob.a

QUADPACK           := $(wildcard quadpack/*.f)
QUADPACK_TARGETS   := $(QUADPACK:%.f=%)
QUADPACK_OBJS      := $(QUADPACK:%.f=%.o)

MT19937_NUMPY_SCRS := $(wildcard mt19937_numpy/*.c)
MT19937_NUMPY_TARS := $(MT19937_NUMPY_SCRS:%.c=%)
MT19937_NUMPY_OBJS := $(MT19937_NUMPY_SCRS:%.c=%.o)
MT19937_NUMPY_INC  := -I./mt19937_numpy

SHELL              := /bin/bash

PYTHON             := python3


all: vp_benchmark_default

help: info

info:
	@echo -e "Quadpack Targets\n----------------"
	@$(foreach var,$(QUADPACK_TARGETS),echo $(var);)
	@echo
	@echo -e "NumPy MT19937 Targets\n---------------------"
	@$(foreach var,$(MT19937_NUMPY_TARS),echo $(var);)

clean:
	@rm -rf *.o *.x __pycache__ data.py quadpack/*.o cephes/*.o cephes/libprob.a *.csv *.csv *.log mt19937_numpy/*.o
	@$(MAKE) -C openspecfun clean --no-print-directory

build_quadpack:
	@$(foreach var,$(QUADPACK_TARGETS),gfortran -c $(var).f -o $(var).o ;)

build_cephes:
	@$(MAKE) -C cephes --no-print-directory

build_amos:
	@$(MAKE) -C openspecfun --no-print-directory

build_numpy_mt19937:
	@$(foreach var,$(MT19937_NUMPY_TARS),$(CC) -c $(var).c -o $(var).o $(CFLAGS) $(LDFLAGS);)

vp_benchmark_rngstreams: build_quadpack build_cephes build_amos
	@echo "Building and Running vp benchmark with RngStreams URNG"
	@echo
	@$(CC) $(QUADPACK_OBJS) -o vp_benchmark_rngstreams.x vp_benchmark.c -DUSE_RNGSTREAMS $(CEPHES) $(OPENSPECFUN) $(UNURAN) $(RNGSTREAMS) $(CFLAGS) $(LDFLAGS)
	@./vp_benchmark_rngstreams.x

vp_benchmark_default: build_quadpack build_cephes build_amos build_numpy_mt19937
	@echo "Building and Running vp benchmark with the MT19937 URNG"
	@echo
	@$(CC) $(QUADPACK_OBJS) $(MT19937_NUMPY_OBJS) -o vp_benchmark_default.x vp_benchmark.c $(MT19937_NUMPY_INC) $(CEPHES) $(OPENSPECFUN) $(UNURAN) $(CFLAGS) $(LDFLAGS)
	@./vp_benchmark_default.x
