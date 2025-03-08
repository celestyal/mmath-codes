# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# gfortran flags
COMPILERFLAGS += -ffree-line-length-0 -O3 -fimplicit-none  -Wall -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface -fbounds-check -Wunused-parameter -fwhole-file -fcheck=all -pedantic -fbacktrace -fdefault-real-8 -fdefault-integer-8 -g

# ifort flags
#COMPILERFLAGS += -r8 -g -traceback -check all -fp-stack-check -fast -r8

# ifx flags
#COMPILERFLAGS += -r8 -xHost -O3 -check all,nouninit -traceback -g

# Project name
NAME := flux_2d

# Configuration settings
FC := gfortran
#FC := ifort
#FC := ifx
AR := ar rcs
LD := $(FC)
RM := rm -f

# List of all source files
SRCS := flux2d_io.f90 flux2d_solv.f90

TEST_SRCS := flux2d.f90 flux2d_gen_fields.f90

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
flux2d_io.mod := flux2d_io.f90.o
flux2d_solv.mod := flux2d_solv.f90.o

flux2d.f90.o: $(flux2d_io.mod)
flux2d.f90.o: $(flux2d_solv.mod)

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) $(LIB) $(wildcard *.mod)
