# Makefile for rm_synthesis (Fortran RM-synthesis package)
# Quick build without CMake - targets: make, make clean, make install

.PHONY: all clean clean-all install uninstall help build_dir

# Compiler and flags
FC := gfortran
FFLAGS := -std=gnu -fallow-argument-mismatch -ffree-line-length-none
OPTFLAGS := -O3 -march=native
DEBUGFLAGS := -g -fbacktrace -fbounds-check

# Build mode: release or debug
MODE ?= release
ifeq ($(MODE),debug)
  FFLAGS += $(DEBUGFLAGS)
else
  FFLAGS += $(OPTFLAGS)
endif

# Optional OpenMP support (set OMP=1 to enable)
OMP ?= 0
ifeq ($(OMP),1)
	FFLAGS += -fopenmp
endif

# Mode/OMP specific artifact tag so build outputs do not conflict
MODE_TAG := $(MODE)_omp$(OMP)

# Directories
SRCDIR := src
BUILDDIR := build/$(MODE_TAG)
BINDIR := bin
MODDIR := $(BUILDDIR)/modules

# CFITSIO library
CFITSIO_LIB ?= -lcfitsio
LIBS := $(CFITSIO_LIB)

# Source files
MODSRC := $(SRCDIR)/rm_synthesis_mod.f90
MAINSRC := $(SRCDIR)/rm_synthesis.f
INCSRC := $(SRCDIR)/myfits_info.f $(SRCDIR)/printerror.f

SOURCES := $(MODSRC) $(MAINSRC) $(INCSRC)
OBJFILES := $(BUILDDIR)/rm_synthesis_mod.o $(BUILDDIR)/rm_synthesis.o

# Target executable (mode-specific plus default convenience path)
EXECUTABLE_MODE := $(BINDIR)/rm_synthesis_$(MODE_TAG)
EXECUTABLE := $(BINDIR)/rm_synthesis

# Default target
all: $(EXECUTABLE)

$(BUILDDIR):
	@mkdir -p $(BUILDDIR) $(MODDIR) $(BINDIR)

# Module compilation
$(BUILDDIR)/rm_synthesis_mod.o: $(MODSRC) | $(BUILDDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

# Main program compilation
$(BUILDDIR)/rm_synthesis.o: $(MAINSRC) $(BUILDDIR)/rm_synthesis_mod.o | $(BUILDDIR)
	$(FC) $(FFLAGS) -I$(MODDIR) -J$(MODDIR) -c $< -o $@

# Linking
$(EXECUTABLE_MODE): $(OBJFILES) | $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)
	@cp -f $@ $(EXECUTABLE)
	@echo "✓ Executable created: $@"
	@echo "✓ Updated default executable: $(EXECUTABLE)"

$(EXECUTABLE): $(EXECUTABLE_MODE)
	@:

clean:
	@rm -rf $(BUILDDIR)
	@rm -f $(EXECUTABLE_MODE) $(EXECUTABLE)
	@echo "✓ Cleaned artifacts for mode tag: $(MODE_TAG)"

clean-all:
	@rm -rf build
	@rm -f $(BINDIR)/rm_synthesis $(BINDIR)/rm_synthesis_*
	@echo "✓ Cleaned all build artifacts for every mode"

install: $(EXECUTABLE)
	@install -d /usr/local/bin
	@install -m 755 $(EXECUTABLE) /usr/local/bin/
	@install -d /usr/local/share/rm_synthesis
	@cp -r cfg /usr/local/share/rm_synthesis/
	@echo "✓ Installed to /usr/local/bin/rm_synthesis"

uninstall:
	@rm -f /usr/local/bin/rm_synthesis
	@rm -rf /usr/local/share/rm_synthesis
	@echo "✓ Uninstalled"

help:
	@echo "RM-Synthesis Build System"
	@echo "========================="
	@echo "Usage: make [target] [MODE=debug|release] [OMP=0|1]"
	@echo ""
	@echo "Targets:"
	@echo "  make                         - Build executable (default, release mode)"
	@echo "  make MODE=debug              - Build with debug symbols and checks"
	@echo "  make OMP=1                   - Build with OpenMP enabled"
	@echo "  make clean [MODE=.. OMP=..]  - Remove artifacts for selected mode/OMP"
	@echo "  make clean-all               - Remove all mode/OMP build artifacts"
	@echo "  make install      - Install to /usr/local/bin"
	@echo "  make uninstall    - Remove installation"
	@echo "  make help         - Show this message"
	@echo ""
	@echo "Note: Artifacts are mode-specific under build/<mode>_omp<0|1>."
	@echo "      Switching MODE/OMP does not require make clean."
	@echo ""
	@echo "Examples:"
	@echo "  make                            # Build release version"
	@echo "  make MODE=debug                 # Build debug version"
	@echo "  make MODE=release OMP=1         # Build OpenMP-enabled release version"
	@echo "  make MODE=debug OMP=1           # Build OpenMP-enabled debug version"
	@echo "  make clean MODE=debug OMP=1     # Clean only debug+OMP artifacts"
	@echo "  make clean-all                  # Clean everything"
	@echo "  make install                    # Install to system"
	@echo "  CFITSIO_LIB=-lcfitsio make  # Specify CFITSIO library"
