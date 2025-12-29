# Note: setup.f90 should be compiled and executed separately
# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# gfortran flags
COMPILERFLAGS += -ffree-line-length-0 -O3 -fimplicit-none  -Wall -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface -fbounds-check -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -fdefault-real-8 -fdefault-integer-8 -g

# ifort flags
#COMPILERFLAGS += -r8 -g -traceback -check all -fp-stack-check -fast -r8

# ifx flags
#COMPILERFLAGS += -r8 -xHost -O3 -check all,nouninit -traceback -g

# Project name
NAME := cartesian_flux_transport

# Configuration settings
FC := gfortran
#FC := ifort
#FC := ifx
AR := ar rcs
LD := $(FC)
RM := rm -f

# List of all source files
SRCS := src/cartesian-flux-transport/io.f90 \
	src/cartesian-flux-transport/kernel.f90

TEST_SRCS := src/cartesian-flux-transport/solver.f90 \
	     src/cartesian-flux-transport/tests/tests_find_dt.f90
	     
# Add source and tests directories to search paths
vpath % .: src/cartesian-flux-transport
vpath % .: src/cartesian-flux-transport/tests

# Define a map from each file name to its object file
obj = $(src).o
$(foreach src, $(SRCS) $(TEST_SRCS), $(eval $(src) := $(obj)))

# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
TEST_OBJS := $(addsuffix .o, $(TEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
TEST_EXE := $(patsubst %.f90, %.exe, $(TEST_SRCS))

# Declare all public targets
.PHONY: all clean
all: $(LIB) $(TEST_EXE)

# Create the static library from the object files
$(LIB): $(OBJS)
	$(AR) $@ $^

# Link the test executables
$(TEST_EXE): %.exe: %.f90.o $(LIB)
	$(LD) -o $@ $^

# Create object files from Fortran source
$(OBJS) $(TEST_OBJS): %.o: %
	$(FC) $(COMPILERFLAGS) -c -o $@ $<

# Define all module interdependencies
io.mod := src/cartesian-flux-transport/io.f90
kernel.mod := src/cartesian-flux-transport/kernel.f90

$(src/cartesian-flux-transport/solver.f90): $(io.mod)
$(src/cartesian-flux-transport/solver.f90): $(kernel.mod)
$(src/cartesian-flux-transport/tests/tests_find_dt.f90): $(kernel.mod)

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) $(LIB) $(wildcard *.mod)
