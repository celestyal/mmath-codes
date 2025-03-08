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
SRCS := io.f90 kernel.f90

TEST_SRCS := solver.f90

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
io.mod := io.f90.o
kernel.mod := kernel.f90.o

solver.f90.o: $(io.mod)
solver.f90.o: $(kernel.mod)

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) $(LIB) $(wildcard *.mod)
