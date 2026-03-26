# Author: Arthur Murphy
# Contributors: Itxaso Muñoz-Aldalur, Oliwier Misztal

FC     = gfortran
FFLAGS = -O2 -Wall -Wextra

# Make sure that MPI is installed
# 1. Search for mpi.mod dynamically in common linux paths, taking only the first match
# 2. Check if the file was found; if empty, throw an error (exiting make)
# 3. Extract just the directory path from the file path & add it with -I to compiler flags

MPI_MOD_FILE := $(shell find /usr /opt -name mpi.mod -type f 2>/dev/null | head -n 1)

ifeq ($(strip $(MPI_MOD_FILE)),)
$(error "MPI is required but mpi.mod was not found! Please install OpenMPI (e.g. 'sudo apt install libopenmpi-dev')")
endif

MPI_INC_DIR := $(dir $(MPI_MOD_FILE))
FFLAGS += -I$(MPI_INC_DIR)


OBJ_DIR = ../bin
EXE     = $(OBJ_DIR)/main_serial.x

# List object files taking into account dependencies
LIB_OBJ = $(OBJ_DIR)/parameters.o \
          $(OBJ_DIR)/io.o \
          $(OBJ_DIR)/initial-conf.o \
          $(OBJ_DIR)/energy.o \
          $(OBJ_DIR)/energy_all_atoms.o \
          $(OBJ_DIR)/monte_carlo.o \
          $(OBJ_DIR)/observables.o

MAIN_OBJ = $(OBJ_DIR)/main_serial.o

all: $(EXE)

$(EXE): $(LIB_OBJ) $(MAIN_OBJ)
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -o $@ $(LIB_OBJ) $(MAIN_OBJ)

# Pattern rule for modules inside lib/
$(OBJ_DIR)/%.o: lib/%.f90
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@

# Explicit dependencies
$(OBJ_DIR)/energy.o: lib/energy.f90 \
                     $(OBJ_DIR)/parameters.o

$(OBJ_DIR)/energy_all_atoms.o: lib/energy_all_atoms.f90 \
                               $(OBJ_DIR)/parameters.o

$(OBJ_DIR)/io.o: lib/io.f90 \
                 $(OBJ_DIR)/parameters.o

$(OBJ_DIR)/initial-conf.o: lib/initial-conf.f90 \
                           $(OBJ_DIR)/parameters.o

$(OBJ_DIR)/monte_carlo.o: lib/monte_carlo.f90 \
                          $(OBJ_DIR)/parameters.o \
                          $(OBJ_DIR)/initial-conf.o \
                          $(OBJ_DIR)/energy.o \
                          $(OBJ_DIR)/energy_all_atoms.o

$(OBJ_DIR)/observables.o: lib/observables.f90 \
                          $(OBJ_DIR)/energy.o \
                          $(OBJ_DIR)/energy_all_atoms.o

# Main program depends on everything
$(OBJ_DIR)/main_serial.o: main_serial.f90 $(LIB_OBJ)
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FFLAGS) -J$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@

.PHONY: all clean run figures pipeline

run: all
	@mkdir -p ../results
	$(EXE)

figures:
	python lib/plot_results.py

clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(EXE)

pipeline:
	$(MAKE) run
	$(MAKE) figures

