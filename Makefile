# Makefile for rm_synthesis (Fortran RM-synthesis package)
# Quick build without CMake - targets: make, make clean, make install

.PHONY: all clean install help build_dir

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

# Directories
SRCDIR := src
BUILDDIR := build
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

# Target executable
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
$(EXECUTABLE): $(OBJFILES) | $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)
	@echo "✓ Executable created: $@"

clean:
	@rm -rf $(BUILDDIR) $(BINDIR)
	@echo "✓ Build artifacts cleaned"

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
	@echo "Usage: make [target] [MODE=debug|release]"
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build executable (default, release mode)"
	@echo "  make MODE=debug   - Build with debug symbols and checks"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make install      - Install to /usr/local/bin"
	@echo "  make uninstall    - Remove installation"
	@echo "  make help         - Show this message"
	@echo ""
	@echo "Examples:"
	@echo "  make                     # Build release version"
	@echo "  make MODE=debug          # Build debug version"
	@echo "  make install             # Install to system"
	@echo "  CFITSIO_LIB=-lcfitsio make  # Specify CFITSIO library"

.PHONY: help
