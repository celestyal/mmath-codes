# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# gfortran flags
COMPILERFLAGS += -ffree-line-length-0 -O3 -fimplicit-none  -Wall -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface -fbounds-check  -Wunused-parameter -fwhole-file -fcheck=all -pedantic  -fbacktrace -fdefault-real-8

# ifort flags
#COMPILERFLAGS += -r8 -g -traceback -check all -fp-stack-check #-fast -r8

# ifx flags
#COMPILERFLAGS +=  -r8 -xHost -O3

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
SRCS := mag_grid.f90 \
				file_io.f90 \
				vector_potential.f90 \
				numerical_solver.f90 \
				types.f90 \
				velocity.f90
TEST_SRCS := flux_2d.f90

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
mag_grid.mod := mag_grid.f90.o
file_io.mod := file_io.f90.o
vector_potential.mod := vector_potential.f90.o
numerical_solver.mod := numerical_solver.f90.o
types.mod := types.f90.o
velocity.mod := velocity.f90.o

file_io.f90.o: $(types.mod)
numerical_solver.f90.o: $(file_io.mod)

flux_2d.f90.o: $(mag_grid.mod)
flux_2d.f90.o: $(file_io.mod)
flux_2d.f90.o: $(vector_potential.mod)
flux_2d.f90.o: $(numerical_solver.mod)
flux_2d.f90.o: $(types.mod)
flux_2d.f90.o: $(velocity.mod)


# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) $(LIB) $(wildcard *.mod)
