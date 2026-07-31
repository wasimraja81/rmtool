#!/bin/bash -l
#SBATCH --partition=askaprt
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=128
#SBATCH --account=askaprt
#SBATCH --time=12:00:00
#SBATCH --mem=0
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --job-name=rm-synthesis
#SBATCH --output="rm-synthesis_%j.out"
#SBATCH --export=NONE

module use /askapbuffer/payne/raj030/askaprtModules
module load singularity/4.1.0-askap 
module load rmtool/4.0

# Create the cfg file: 
CFG_PATH=/scratch/ja3/raj030/for_Jennifer/rm-synthesis/cfg
outDir=/scratch/ja3/raj030/for_Jennifer/rm-synthesis/scratch
cfgFile=${CFG_PATH}/rmsymthesis-jennifer.cfg

# Cleanup old files: 

echo "INFO - Cleaning up old files: "
rm -rf ${outDir}/JENNIFER_TOO_FULLIM_*

cat >${cfgFile} <<EOFcfg
# RM-synthesis test config: full-image run
# Same limited RM settings as subimage test.
# Adjust path/infile names if your data location differs.

#path=/scratch/ja3/raj030/for_Jennifer/rm-synthesis/data/
path=/scratch/ja3/westjl/80775/fullcubes/

#infileQ=cutout-stokesQ.fits
#infileU=cutout-stokesU.fits
infileQ=image.restored.q.EMU_1615-64.SB80775.contcube.conv.fits
infileU=image.restored.u.EMU_1615-64.SB80775.contcube.conv.fits
outfile=${outDir}/JENNIFER_TOO_FULLIM_CONV

remove_badchan=y
global_badchan_file=/scratch/ja3/raj030/for_Jennifer/rm-synthesis/cfg/askap_nan_channels.burdies

# Optional mask controls:
# mask_input_cube_file: input frequency-mask cube (0=bad, >0 good).
# mask_cube_file: override output mask filename (default: <outfile>.MASK.CUBE.FITS).
# mask_trust_mode: safe|strict. safe keeps NaN checks enabled.
mask_input_cube_file=
mask_cube_file=
mask_trust_mode=safe
write_mask_output=n
write_nvalid_output=y
use_gpu=y

subim=n
subim_ra_blc=1
subim_ra_trc=0
subim_ra_inc=1
subim_dec_blc=1
subim_dec_trc=0
subim_dec_inc=1
subim_chan_blc=1
subim_chan_trc=0
subim_chan_inc=1

rem_mean=0
remove_qu_bias=n
resiQ=0.0
slopeQ=0.0
resiU=0.0
slopeU=0.0

#path_I=/scratch/ja3/raj030/for_Jennifer/rm-synthesis/data/
path_I=/scratch/ja3/westjl/80775/fullcubes/image.restored.i.EMU_1615-64.SB80775.contcube.conv.fits
infileI="" 

ofac=1
fac=3.14159265358979
use_auto_rm_range=0
output_mode=ap
ap_angle_mode=phase
beg_rm=-500.0
end_rm=500.0
nrm=101
dry_run=n


tile_auto=y
tile_ra=4501
tile_dec=1315
mem_frac_ram=0.50
mem_frac_vram=0.70
gpu_vram_mib=6144
io_overlap=y
io_read_threads=4
io_write_threads=4
cubestat=y

# Logging and timing
log_level=debug
timing_enabled=y
timing_tile_enabled=y
timing_io_enabled=y
log_output_file=./rmsynthesis.largeCube.cpu.setonix.4.0.log
timing_csv_file=./rmsynthesis.largeCube.cpu.setonix.4.0.csv

EOFcfg

#rm_synthesis_OMP1_GPU0 ${cfgFile}

# Set OpenMP threads explicitly (override for scaling sweeps as needed).
export OMP_NUM_THREADS=8
export OMP_THREAD_LIMIT=${OMP_NUM_THREADS}
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1

echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT}"
echo "OMP_DYNAMIC=${OMP_DYNAMIC}"

# Bind OpenMP threads to physical cores (Prevents thread migration)
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_DISPLAY_ENV=VERBOSE

echo "OMP_PLACES=${OMP_PLACES}"
echo "OMP_PROC_BIND=${OMP_PROC_BIND}"

# Forward OMP controls into container runtime environments.
export SINGULARITYENV_OMP_NUM_THREADS=${OMP_NUM_THREADS}
export SINGULARITYENV_OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT}
export SINGULARITYENV_OMP_DYNAMIC=${OMP_DYNAMIC}
export SINGULARITYENV_OMP_MAX_ACTIVE_LEVELS=${OMP_MAX_ACTIVE_LEVELS}
export SINGULARITYENV_OMP_PLACES=${OMP_PLACES}
export SINGULARITYENV_OMP_PROC_BIND=${OMP_PROC_BIND}
export SINGULARITYENV_OMP_DISPLAY_ENV=${OMP_DISPLAY_ENV}

export APPTAINERENV_OMP_NUM_THREADS=${OMP_NUM_THREADS}
export APPTAINERENV_OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT}
export APPTAINERENV_OMP_DYNAMIC=${OMP_DYNAMIC}
export APPTAINERENV_OMP_MAX_ACTIVE_LEVELS=${OMP_MAX_ACTIVE_LEVELS}
export APPTAINERENV_OMP_PLACES=${OMP_PLACES}
export APPTAINERENV_OMP_PROC_BIND=${OMP_PROC_BIND}
export APPTAINERENV_OMP_DISPLAY_ENV=${OMP_DISPLAY_ENV}

# Use container-provided wrapper binary (rmtool module / SIF entry).
RM_BIN=rm_synthesis_OMP1_GPU0
echo "RM_BIN=${RM_BIN}"

# Keep execution in Slurm cpuset and enforce core binding.
# Even on an exclusive node, constrain this step to OMP_NUM_THREADS CPUs.
srun --cpu-bind=cores --hint=nomultithread \
		 --cpus-per-task=${OMP_NUM_THREADS} --exact --export=ALL \
	env OMP_NUM_THREADS=${OMP_NUM_THREADS} \
			OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT} \
			OMP_DYNAMIC=${OMP_DYNAMIC} \
			OMP_MAX_ACTIVE_LEVELS=${OMP_MAX_ACTIVE_LEVELS} \
			OMP_PLACES=${OMP_PLACES} \
			OMP_PROC_BIND=${OMP_PROC_BIND} \
			OMP_DISPLAY_ENV=${OMP_DISPLAY_ENV} \
			SINGULARITYENV_OMP_NUM_THREADS=${SINGULARITYENV_OMP_NUM_THREADS} \
			SINGULARITYENV_OMP_THREAD_LIMIT=${SINGULARITYENV_OMP_THREAD_LIMIT} \
			SINGULARITYENV_OMP_DYNAMIC=${SINGULARITYENV_OMP_DYNAMIC} \
			SINGULARITYENV_OMP_MAX_ACTIVE_LEVELS=${SINGULARITYENV_OMP_MAX_ACTIVE_LEVELS} \
			SINGULARITYENV_OMP_PLACES=${SINGULARITYENV_OMP_PLACES} \
			SINGULARITYENV_OMP_PROC_BIND=${SINGULARITYENV_OMP_PROC_BIND} \
			SINGULARITYENV_OMP_DISPLAY_ENV=${SINGULARITYENV_OMP_DISPLAY_ENV} \
			APPTAINERENV_OMP_NUM_THREADS=${APPTAINERENV_OMP_NUM_THREADS} \
			APPTAINERENV_OMP_THREAD_LIMIT=${APPTAINERENV_OMP_THREAD_LIMIT} \
			APPTAINERENV_OMP_DYNAMIC=${APPTAINERENV_OMP_DYNAMIC} \
			APPTAINERENV_OMP_MAX_ACTIVE_LEVELS=${APPTAINERENV_OMP_MAX_ACTIVE_LEVELS} \
			APPTAINERENV_OMP_PLACES=${APPTAINERENV_OMP_PLACES} \
			APPTAINERENV_OMP_PROC_BIND=${APPTAINERENV_OMP_PROC_BIND} \
			APPTAINERENV_OMP_DISPLAY_ENV=${APPTAINERENV_OMP_DISPLAY_ENV} \
			${RM_BIN} ${cfgFile}

# Legacy direct launch kept for reference:
# env OMP_NUM_THREADS=${OMP_NUM_THREADS} OMP_PLACES=cores OMP_PROC_BIND=spread rm_synthesis_OMP1_GPU0 ${cfgFile}
