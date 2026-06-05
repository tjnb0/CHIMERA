# Fortran compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall -ffree-line-length-none -march=native -funroll-loops -fopenmp \
         -I/usr/lib64/gfortran/modules
MODDIR = mod
OBJDIR = obj
SRCDIR = src
EXE = mhd_sim

# HDF5 linking flags
HDF5_LIBS = -L/usr/lib64 -lhdf5_fortran -lhdf5

# Detect the OS
ifeq ($(OS),Windows_NT)
    # Windows
    CREATE_DIR = @if not exist $(1) mkdir $(1)
    RM = del /Q
    EXE_SUFFIX = .exe
else
    # Linux
    CREATE_DIR = @mkdir -p $(1)
    RM = rm -f
    EXE_SUFFIX = 
endif

# Define object files
OBJS = \
	$(OBJDIR)/mhd_config.o \
	$(OBJDIR)/mhd_write_h5.o \
	$(OBJDIR)/mhd_field_ops.o \
	$(OBJDIR)/mhd_change_states.o \
	$(OBJDIR)/mhd_derivatives.o \
	$(OBJDIR)/mhd_ghost_bcs.o \
	$(OBJDIR)/mhd_flux.o \
	$(OBJDIR)/mhd_init.o \
	$(OBJDIR)/main.o 

# Default target
all: $(EXE)

# Link final executable
$(EXE)$(EXE_SUFFIX): $(OBJS)
	$(CREATE_DIR) $(MODDIR)
	$(CREATE_DIR) $(OBJDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ $^ $(HDF5_LIBS)

# Dependencies
$(OBJDIR)/main.o: $(SRCDIR)/main.f90 $(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_init.o \
                  $(OBJDIR)/mhd_field_ops.o $(OBJDIR)/mhd_change_states.o \
                  $(OBJDIR)/mhd_derivatives.o $(OBJDIR)/mhd_flux.o $(OBJDIR)/mhd_write_h5.o \
				  $(OBJDIR)/mhd_ghost_bcs.o 
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_init.o: $(SRCDIR)/mhd_init.f90 $(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_field_ops.o
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_flux.o: $(SRCDIR)/mhd_flux.f90 $(OBJDIR)/mhd_field_ops.o $(OBJDIR)/mhd_config.o 
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_derivatives.o: $(SRCDIR)/mhd_derivatives.f90 
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_config.o: $(SRCDIR)/mhd_config.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_write_h5.o: $(SRCDIR)/mhd_write_h5.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_field_ops.o: $(SRCDIR)/mhd_field_ops.f90 $(OBJDIR)/mhd_config.o
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_ghost_bcs.o: $(SRCDIR)/mhd_ghost_bcs.f90 $(OBJDIR)/mhd_config.o
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_change_states.o: $(SRCDIR)/mhd_change_states.f90 $(OBJDIR)/mhd_config.o
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@


# Clean rule
clean:
	$(RM) $(EXE)$(EXE_SUFFIX)
	@if [ -d $(OBJDIR) ]; then $(RM) $(OBJDIR)/*.o; fi
	@if [ -d $(MODDIR) ]; then $(RM) $(MODDIR)/*.mod; fi

.PHONY: all clean