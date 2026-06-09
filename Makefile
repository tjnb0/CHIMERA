# Fortran compiler and flags
FC     = gfortran
FFLAGS = -O3 -Wall -ffree-line-length-none -march=native -funroll-loops -fopenmp \
         -I/usr/lib64/gfortran/modules

MODDIR = mod
OBJDIR = obj
SRCDIR = src
OUTDIR = outputs
EXE    = chimera

# HDF5 linking flags
HDF5_LIBS = -L/usr/lib64 -lhdf5_fortran -lhdf5

# OS detection
ifeq ($(OS),Windows_NT)
    RM         = del /Q
    EXE_SUFFIX = .exe
    MKDIR      = if not exist "$(1)" mkdir "$(1)"
else
    RM         = rm -f
    EXE_SUFFIX =
    MKDIR      = mkdir -p $(1)
endif

# Object files (in dependency order)
OBJS = \
	$(OBJDIR)/mhd_config.o        \
	$(OBJDIR)/mhd_write_h5.o      \
	$(OBJDIR)/mhd_field_ops.o     \
	$(OBJDIR)/mhd_change_states.o \
	$(OBJDIR)/mhd_derivatives.o   \
	$(OBJDIR)/mhd_bc.o            \
	$(OBJDIR)/mhd_flux.o          \
	$(OBJDIR)/mhd_init.o          \
	$(OBJDIR)/main.o

# Test executables (no HDF5 required)
TESTDIR   = tests/fortran
TEST_EXES = \
	$(TESTDIR)/test_change_states \
	$(TESTDIR)/test_field_ops     \
	$(TESTDIR)/test_bcs

# -- Targets -----------------------------------------------------------------

.PHONY: all clean tests

all: $(EXE)$(EXE_SUFFIX)

# Directory creation - order-only prerequisites.
# Make checks these exist before any rule that lists them after |
# but does not use their timestamps to decide whether to rebuild.
$(OBJDIR)/ $(MODDIR)/ $(OUTDIR)/:
	@$(call MKDIR,$@)

# Link
$(EXE)$(EXE_SUFFIX): $(OBJS)
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ $^ $(HDF5_LIBS)

# -- Compile rules (| means order-only: directory must exist, not timestamped) --

$(OBJDIR)/mhd_config.o: $(SRCDIR)/mhd_config.f90 | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_write_h5.o: $(SRCDIR)/mhd_write_h5.f90 | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_field_ops.o: $(SRCDIR)/mhd_field_ops.f90 $(OBJDIR)/mhd_config.o | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_change_states.o: $(SRCDIR)/mhd_change_states.f90 $(OBJDIR)/mhd_config.o | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_derivatives.o: $(SRCDIR)/mhd_derivatives.f90 $(OBJDIR)/mhd_config.o | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_bc.o: $(SRCDIR)/mhd_bc.f90 $(OBJDIR)/mhd_config.o | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_flux.o: $(SRCDIR)/mhd_flux.f90 $(OBJDIR)/mhd_field_ops.o $(OBJDIR)/mhd_config.o | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/mhd_init.o: $(SRCDIR)/mhd_init.f90 $(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_field_ops.o | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(OBJDIR)/main.o: $(SRCDIR)/main.f90 $(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_init.o \
                  $(OBJDIR)/mhd_field_ops.o $(OBJDIR)/mhd_change_states.o \
                  $(OBJDIR)/mhd_derivatives.o $(OBJDIR)/mhd_flux.o \
                  $(OBJDIR)/mhd_write_h5.o $(OBJDIR)/mhd_bc.o | $(OBJDIR)/ $(MODDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

# -- Test rules --------------------------------------------------------------

tests: $(TEST_EXES)

$(TESTDIR)/test_change_states: $(TESTDIR)/test_change_states.f90 \
		$(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_change_states.o | $(TESTDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ $< \
		$(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_change_states.o

$(TESTDIR)/test_field_ops: $(TESTDIR)/test_field_ops.f90 \
		$(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_field_ops.o \
		$(OBJDIR)/mhd_derivatives.o | $(TESTDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ $< \
		$(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_field_ops.o \
		$(OBJDIR)/mhd_derivatives.o

$(TESTDIR)/test_bcs: $(TESTDIR)/test_bcs.f90 \
		$(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_bc.o | $(TESTDIR)/
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ $< \
		$(OBJDIR)/mhd_config.o $(OBJDIR)/mhd_bc.o

$(TESTDIR)/:
	@mkdir -p $@

# -- Clean -------------------------------------------------------------------

# Leaves outputs/ intact
clean:
	$(RM) $(EXE)$(EXE_SUFFIX)
	$(RM) $(TEST_EXES)
	@if [ -d $(OBJDIR) ]; then $(RM) $(OBJDIR)/*.o;   fi
	@if [ -d $(MODDIR) ]; then $(RM) $(MODDIR)/*.mod; fi
