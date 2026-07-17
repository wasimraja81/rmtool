# Makefile for rm_synthesis (Fortran RM-synthesis package)
# Quick build without CMake - targets: make, make clean, make install

.PHONY: all clean clean-all install uninstall help build_dir

# Compiler and flags
FC := gfortran
GPU_FC ?= nvfortran
NVFORTRAN_PATH := $(shell command -v nvfortran 2>/dev/null)
GFORTRAN_PATH := $(shell command -v gfortran 2>/dev/null)
CPPFLAGS := -cpp
BASEFLAGS := $(CPPFLAGS) -std=gnu -fallow-argument-mismatch -ffree-line-length-none
CPU_OPTFLAGS := -O3 -march=native
CPU_DEBUGFLAGS := -g -fbacktrace -fbounds-check
CPU_PROFILEFLAGS := -O3 -march=native -g -fno-omit-frame-pointer
CPU_OMPFLAGS := -fopenmp
GPU_NVFLAGS := -cpp -O3 -mp=gpu -gpu=cc80,managed -DUSE_GPU
GPU_GNUFLAGS := $(BASEFLAGS) -O3 -fopenmp -foffload=nvptx-none -foffload="-lm" -ffast-math -fno-finite-math-only -DUSE_GPU

FFLAGS := $(BASEFLAGS)

# Build mode: release, profile, or debug
MODE ?= release
# Optional OpenMP support (set OMP=1 to enable)
OMP ?= 0
OMP_EFFECTIVE := $(OMP)
HOST_OMP_CPP := -DHOST_OMP=$(OMP_EFFECTIVE)

# Optional GPU/offload build (set GPU=1 to enable)
GPU ?= 0

# Auto-select GPU compiler only when user did not explicitly set GPU_FC.
# Preference order: nvfortran, then gfortran.
ifeq ($(GPU),1)
  ifeq ($(origin GPU_FC),file)
    ifneq ($(NVFORTRAN_PATH),)
      GPU_FC := nvfortran
    else ifneq ($(GFORTRAN_PATH),)
      GPU_FC := gfortran
    endif
  endif
endif

ifeq ($(GPU_FC),gfortran)
  GPUFLAGS := $(GPU_GNUFLAGS)
else
  GPUFLAGS := $(GPU_NVFLAGS)
endif

ifeq ($(GPU),1)
  FC := $(GPU_FC)
  FFLAGS := $(GPUFLAGS)
else
  ifeq ($(MODE),debug)
    FFLAGS += $(CPU_DEBUGFLAGS)
	else ifeq ($(MODE),profile)
		FFLAGS += $(CPU_PROFILEFLAGS)
  else
    FFLAGS += $(CPU_OPTFLAGS)
  endif
  ifeq ($(OMP_EFFECTIVE),1)
    FFLAGS += $(CPU_OMPFLAGS)
  endif
endif
FFLAGS += $(HOST_OMP_CPP)

# Human-readable build flavor naming
ifeq ($(GPU),0)
	ifeq ($(OMP_EFFECTIVE),1)
		FLAVOR := cpu_omp
	else
		FLAVOR := cpu_serial
	endif
else
	ifeq ($(OMP_EFFECTIVE),1)
		FLAVOR := gpu_offload_hostomp
	else
		FLAVOR := gpu_offload
	endif
endif

# Mode/flavor specific artifact tag so build outputs do not conflict
BUILD_TAG := $(MODE)_$(FLAVOR)

# Directories
SRCDIR := src
BUILDDIR := build/$(BUILD_TAG)
BINDIR ?= bin
MODDIR := $(BUILDDIR)/modules
PROFILE_BINDIR := scratch/profiles/bin

# Profile builds default to a dedicated scratch binary directory.
# Users can still override BINDIR explicitly on the command line.
ifeq ($(MODE),profile)
ifneq ($(origin BINDIR),command line)
BINDIR := $(PROFILE_BINDIR)
endif
endif

# Default bin directory (used to decide whether to update the convenience symlink)
DEFAULT_BINDIR := bin

# CFITSIO library
CFITSIO_LIB ?= -lcfitsio
# -lpthread: async tile-write (io_overlap) dispatches writes on a raw POSIX
# thread outside the OpenMP runtime. A no-op on glibc >= 2.34 (pthread is
# folded into libc there) but kept explicit for portability to older libc.
LIBS := $(CFITSIO_LIB) -lpthread

# Source files
MODSRC := $(SRCDIR)/rm_synthesis_mod.f90
MAINSRC := $(SRCDIR)/rm_synthesis.f90
INCSRC := $(SRCDIR)/myfits_info.f $(SRCDIR)/printerror.f

SOURCES := $(MODSRC) $(MAINSRC) $(INCSRC)
OBJFILES := $(BUILDDIR)/rm_synthesis_mod.o $(BUILDDIR)/rm_synthesis.o

# Target executable (mode-specific plus default convenience path)
EXECUTABLE_MODE := $(BINDIR)/rm_synthesis_$(BUILD_TAG)
EXECUTABLE := $(BINDIR)/rm_synthesis

# Default target
all: $(EXECUTABLE)

ifeq ($(GPU),1)
CHECK_GPU_COMPILER := check_gpu_compiler
else
CHECK_GPU_COMPILER :=
endif

check_gpu_compiler:
	@command -v $(FC) >/dev/null 2>&1 || \
	  { echo "ERROR: GPU compiler '$(FC)' not found in PATH."; \
	    echo "       Auto-select order (when GPU_FC not set): nvfortran -> gfortran"; \
	    echo "       Set GPU_FC=<compiler> explicitly, e.g. GPU_FC=gfortran or GPU_FC=nvfortran."; \
	    echo "       nvfortran uses flags: $(GPU_NVFLAGS)"; \
	    echo "       gfortran uses flags:  $(GPU_GNUFLAGS)"; \
	    exit 127; }

$(BUILDDIR):
	@mkdir -p $(BUILDDIR) $(MODDIR)

$(BINDIR):
	@mkdir -p $(BINDIR)

# Module compilation
$(BUILDDIR)/rm_synthesis_mod.o: $(MODSRC) | $(BUILDDIR) $(CHECK_GPU_COMPILER)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

# Main program compilation
$(BUILDDIR)/rm_synthesis.o: $(MAINSRC) $(BUILDDIR)/rm_synthesis_mod.o | $(BUILDDIR) $(CHECK_GPU_COMPILER)
	$(FC) $(FFLAGS) -I$(MODDIR) -J$(MODDIR) -c $< -o $@

# Linking
$(EXECUTABLE_MODE): $(OBJFILES) | $(BINDIR) $(CHECK_GPU_COMPILER)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)
	@echo "✓ Executable created: $@"
ifeq ($(BINDIR),$(DEFAULT_BINDIR))
	@cp -f $@ $(EXECUTABLE)
	@echo "✓ Updated default executable: $(EXECUTABLE)"
else
	@echo "  (default bin/rm_synthesis not updated; BINDIR=$(BINDIR))"
endif

$(EXECUTABLE): $(EXECUTABLE_MODE)
	@:

clean:
	@rm -rf $(BUILDDIR)
	@rm -f $(EXECUTABLE_MODE) $(EXECUTABLE)
	@echo "✓ Cleaned artifacts for build tag: $(BUILD_TAG)"

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
	@echo "Usage: make [target] [MODE=release|profile|debug] [OMP=0|1] [GPU=0|1]"
	@echo ""
	@echo "Targets:"
	@echo "  make                         - Build executable (default, release mode)"
	@echo "  make MODE=profile            - Build profiling binary (optimized + symbols + frame pointers)"
	@echo "  make MODE=debug              - Build with debug symbols and checks"
	@echo "  make OMP=1                   - Build with OpenMP enabled CPU backend"
	@echo "  make GPU=1                   - Build GPU/offload backend (auto: nvfortran -> gfortran)"
	@echo "  make GPU=1 GPU_FC=gfortran   - Build GPU/offload backend with GNU offload"
	@echo "  make clean [MODE=.. OMP=.. GPU=..] - Remove artifacts for selected mode/OMP/GPU"
	@echo "  make clean-all               - Remove all mode/OMP/GPU build artifacts"
	@echo "  make install      - Install to /usr/local/bin"
	@echo "  make uninstall    - Remove installation"
	@echo "  make help         - Show this message"
	@echo ""
	@echo "Note: Artifacts are mode-specific under build/<mode>_<flavor>."
	@echo "      Flavors: cpu_serial, cpu_omp, gpu_offload, gpu_offload_hostomp."
	@echo "      MODE=profile defaults binaries to scratch/profiles/bin (unless BINDIR=... is provided)."
	@echo "      GPU=1 OMP=0/1 sets HOST_OMP=0/1; host OpenMP regions are gated accordingly."
	@echo "      GPU and OMP can be enabled together (e.g., GPU=1 OMP=1)."
	@echo "      Switching MODE/OMP/GPU does not require make clean."
	@echo ""
	@echo "Examples:"
	@echo "  make                            # Build release version"
	@echo "  make MODE=profile               # Build CPU profiling-friendly version"
	@echo "  make MODE=debug                 # Build debug version"
	@echo "  make MODE=release OMP=1         # Build OpenMP-enabled CPU release version"
	@echo "  make MODE=profile OMP=1         # Build OpenMP-enabled CPU profiling version"
	@echo "  make MODE=debug OMP=1           # Build OpenMP-enabled CPU debug version"
	@echo "  make GPU=1                      # Build GPU/offload binary (auto compiler)"
	@echo "  make GPU=1 GPU_FC=nvfortran     # Select GPU compiler explicitly"
	@echo "  make GPU=1 GPU_FC=gfortran      # Use GNU OpenMP offload backend"
	@echo "  make clean MODE=debug OMP=1 GPU=0 # Clean only debug+OMP CPU artifacts"
	@echo "  make clean-all                  # Clean everything"
	@echo "  make install                    # Install to system"
	@echo "  CFITSIO_LIB=-lcfitsio make  # Specify CFITSIO library"
