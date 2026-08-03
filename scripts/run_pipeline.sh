#!/usr/bin/env bash
set -euo pipefail

# scripts/run_pipeline.sh -- orchestrator for match_cubes -> rm_synthesis
# -> rmclean_cubes, driven by one small pipeline cfg (cfg/pipeline-
# example.cfg is a documented template). Ticket: the user's own item-4
# ask, alongside item 3 ("run a full test on moderately big real data")
# -- see docs/dev/RMCLEAN_INTEGRATION_PLAN.md's own T4 entry and this
# script's own git history for the design discussion.
#
# ANY non-empty subset of {match, rmsynth, rmclean} is valid (T18,
# docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md) -- not just the full chain.
# Found needed for real: a real multi-band run wanted match_cubes'
# output on one disk (a slow, high-capacity one, alongside the raw
# input) and rm_synthesis'/rmclean_cubes' own output on a different,
# faster one (this script's own single outdir= can't express two
# different disks for one invocation) -- splitting into a match-only
# invocation (outdir on the first disk) followed by a rmsynth[,rmclean]-
# only invocation (outdir on the second, rmsynth_cfg_template's own
# path=/infileQ=/infileU= pointed at the first invocation's real
# output) is exactly this case. Whichever stage would normally receive
# a chained path from an earlier stage uses its own cfg template's own
# value unmodified instead, whenever that earlier stage is not ALSO
# present in this same invocation.
#
# This script invents NO new algorithmic options. Every stage tool
# (match_cubes, rm_synthesis, rmclean_cubes) already has its own complete,
# tested cfg/CLI interface; this script's only job is to chain
# path-shaped values between stages (via sed substitution into a COPY of
# each tool's own native cfg template -- never a new cfg-parsing DSL) and
# drive execution in the right order with the right CPU thread pinning.
#
# --- match_cubes output location (a real constraint, not a gap here) ---
# match_cubes writes each output ALONGSIDE its own input
# (<infile><suffix>, confirmed directly in src/match_cubes.f90's own
# process_one_file_*/outsuffix logic) -- it has no output-directory
# option of its own. This script does NOT relocate its multi-GB outputs
# into outdir afterward (that would double disk I/O and defeat
# match_cubes' own "no intermediate copy" design goal) -- it reads
# match_cubes' own manifest (see below) to learn the resulting paths and
# feeds those straight into the rmsynth stage's own infileQ/infileU.
#
# --- Why Q and U are combined into ONE match_cubes call ---
# match_cubes derives ONE common footprint (reproject) and ONE common
# target beam (convolve, auto-derived as "the smallest beam every good
# channel of every input can be deconvolved from") from whatever set of
# files it's given in a single call. Running it separately for the Q list
# and the U list risks each call auto-deriving a DIFFERENT common target
# beam if Q's and U's own per-channel restoring beams for the same band
# aren't bit-identical (not guaranteed by every imaging pipeline, even
# though the sky footprint always matches). Combining Q+U into one
# infiles= list guarantees both polarizations land on the IDENTICAL
# output grid and beam by construction, at zero extra cost (redundant
# per-band WCS/beam entries are harmless to the auto-derivation).
#
# --- Never infer outcome from filesystem state ---
# Whether a given match_cubes input was skipped (already matched the
# target) or processed is read from match_cubes' own explicit
# manifest=<path> output (one line per input, tab-separated
# "<infile> SKIPPED|PROCESSED <effective_path>"), never guessed from
# whether a file happens to exist on disk. A stray/stale file left by an
# unrelated earlier run could sit at the exact path match_cubes would
# have written and be silently (and wrongly) mistaken for this run's own
# result -- the same class of trap as the historical rmclean_cubes
# SIGSEGV (unchecked FTINIT status on a pre-existing file). Standing rule
# throughout this whole pipeline: a DATA output file this script or any
# stage tool is about to write that already exists is always refused,
# never silently reused or overwritten. This script's own provenance
# directory (see below) is metadata scaffolding, not a data output --
# it is shared across repeated runs, not subject to this rule.
#
# --- rmclean_cubes cfg override order (different from reproject_cubes!) ---
# rmclean_cubes parses --config AFTER any CLI key=value args and applies
# each key unconditionally (src/rmclean_cubes.f90's own parse_args), so
# the FILE wins over CLI for any key it also sets -- this script therefore
# always drives rmclean_cubes via a fully-substituted --config file, never
# CLI overrides alongside --config.
#
# --- Provenance ---
# Every run leaves records in <outdir>/<run_name>.provenance/ (mkdir -p,
# shared across repeated runs of the same outdir/run_name -- the
# directory itself is scaffolding, not an output, so its pre-existence
# is never a reason to abort): a verbatim copy of the pipeline cfg
# actually used, the fully-substituted per-stage cfgs actually fed to
# each tool, the match_cubes manifest (if the match stage ran), a log of
# every external command actually executed, and each stage's own full
# stdout+stderr (captured via tee, so it's both visible live and kept
# permanently), and a run_summary.txt with per-stage timing and the
# final output paths. Every filename in this directory is tagged with
# this run's own start timestamp, so repeated runs never collide or
# overwrite each other's provenance -- only the actual DATA outputs
# (match_cubes/rm_synthesis/rmclean_cubes cubes) are refused-if-exists,
# each by its own stage tool.
#
# Usage: scripts/run_pipeline.sh <pipeline_cfg>
#    or: scripts/run_pipeline.sh --help | -h
#
# See cfg/pipeline-example.cfg for the full, documented key list.

usage() {
  cat <<'EOF'
Usage: run_pipeline.sh <pipeline_cfg>

<pipeline_cfg>: path to a pipeline cfg file (see cfg/pipeline-example.cfg
for the full, documented key list: stages, outdir, run_name, match_*,
rmsynth_cfg_template, rmsynth_backend, rmsynth_omp_threads,
rmclean_cfg_template).

Runs, in this fixed relative order, whichever of these are named in the
cfg's own stages= key -- ANY non-empty subset is valid, not just the
full match,rmsynth,rmclean sequence (e.g. stages=match alone to prepare
data now and run rmsynth separately/later on different storage;
stages=rmclean alone against an existing dirty cube from a prior
rmsynth-only run):
  match    bin/match_cubes (convolve+reproject), Q and U combined in one
           call so both land on the identical output grid/beam
  rmsynth  bin/rm_synthesis_release_cpu_omp|_gpu_offload_hostomp
  rmclean  bin/rmclean_cubes

Whenever a stage's own inputs would normally come from an earlier stage
that is NOT present in this same invocation, that stage's own cfg
template's own path-shaped keys (path=/infileQ=/infileU= for rmsynth;
ampfile=/phafile=/maskfile=/outfile= for rmclean) are used exactly as
you wrote them, instead of being overridden with a chained value.

Each stage's own native cfg template is copied into
<outdir>/<run_name>.provenance/ with just the path-shaped keys
substituted; every other setting in your own template is used
unchanged. Leaves a full provenance record (cfgs actually used, exact
commands run, per-stage logs, run summary) in that same directory --
see this script's own top-of-file comment for the full design
rationale (match's output-location constraint, why Q/U are combined,
the manifest mechanism, rmclean_cubes' cfg-vs-CLI override order).
EOF
}

if [[ $# -eq 0 || "${1}" == "-h" || "${1}" == "--help" ]]; then
  usage
  [[ $# -eq 0 ]] && exit 1
  exit 0
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PIPELINE_CFG="${1}"
if [[ ! -f "${PIPELINE_CFG}" ]]; then
  if [[ -f "${ROOT_DIR}/cfg/$(basename "${PIPELINE_CFG}")" ]]; then
    PIPELINE_CFG="${ROOT_DIR}/cfg/$(basename "${PIPELINE_CFG}")"
  else
    echo "[pipeline] ERROR: cfg not found: ${PIPELINE_CFG}" >&2
    exit 1
  fi
fi
PIPELINE_CFG="$(realpath "${PIPELINE_CFG}")"

# --- cfg_get <file> <key>: first non-comment line matching
# ^[space]*<key>[space]*=[space]*<value>, trimmed. Tolerates both
# "key=value" and "key = value" styles (both already in use across this
# project's own cfg files) -- never `source`s a cfg, only ever reads one
# named field at a time, same discipline scratch/run_rmsynthesis_test.sh
# already uses. ---
cfg_get() {
  local file="$1" key="$2"
  awk -v k="${key}" '
    { line = $0 }
    { t = line; sub(/^[[:space:]]+/, "", t); if (t ~ /^#/ || t == "") next }
    line ~ ("^[[:space:]]*" k "[[:space:]]*=") {
      sub("^[^=]*=", "", line)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
      print line
      exit
    }
  ' "${file}"
}

# --- cfg_set_inplace <file> <key> <value>: replace an existing
# "<key>=..." line (any whitespace style) with "<key>=<value>" exactly
# once. Errors if the key isn't present -- this script only ever
# substitutes keys it already confirmed exist in the template (via
# cfg_get), so a miss here means the template changed shape unexpectedly. ---
cfg_set_inplace() {
  local file="$1" key="$2" value="$3"
  if ! grep -qE "^[[:space:]]*${key}[[:space:]]*=" "${file}"; then
    echo "[pipeline] ERROR: expected key '${key}' not found in ${file}" >&2
    exit 1
  fi
  # '|' delimiter: values here are always paths, which always contain '/'.
  sed -i -E "s|^[[:space:]]*${key}[[:space:]]*=.*|${key}=${value}|" "${file}"
}

flag_is_true() {
  # Same convention as rm_synthesis_mod.f90's own flag_from_value:
  # first non-blank char is 1/y/t (case-insensitive) -> true.
  local v
  v="$(echo "${1}" | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')"
  [[ "${v:0:1}" == "1" || "${v:0:1}" == "y" || "${v:0:1}" == "t" ]]
}

log() { echo "[pipeline] $*"; }
die() { echo "[pipeline] ERROR: $*" >&2; exit 1; }

# --- run_logged <log_file> <cmd...>: record the exact command in
# commands.log (timestamped) BEFORE running it, then run it with stdout
# +stderr captured to <log_file> via tee (still visible live). ---
run_logged() {
  local logfile="$1"; shift
  echo "[$(date -Iseconds)] $*" >> "${COMMANDS_LOG}"
  "$@" 2>&1 | tee "${logfile}"
  return "${PIPESTATUS[0]}"
}

# ============================================================
# Parse the top-level pipeline cfg
# ============================================================
STAGES_RAW="$(cfg_get "${PIPELINE_CFG}" stages)"
[[ -z "${STAGES_RAW}" ]] && die "stages= is required in ${PIPELINE_CFG}"
IFS=',' read -ra STAGE_LIST <<< "${STAGES_RAW}"
DO_MATCH=0; DO_RMSYNTH=0; DO_RMCLEAN=0
for s in "${STAGE_LIST[@]}"; do
  case "${s}" in
    match) DO_MATCH=1 ;;
    rmsynth) DO_RMSYNTH=1 ;;
    rmclean) DO_RMCLEAN=1 ;;
    *) die "unrecognised stage '${s}' in stages= (expected match, rmsynth, rmclean)" ;;
  esac
done
# Any non-empty subset is allowed -- match/rmsynth/rmclean are each
# independently useful (e.g. match alone to prepare data for a later,
# separate rmsynth run on different storage; rmclean alone against an
# existing dirty cube from a prior rmsynth-only run). When a later
# stage's own inputs would normally come from an earlier stage that
# ISN'T present in THIS invocation, that stage's own cfg template's
# own path-shaped keys are used unmodified instead (see the rmsynth/
# rmclean sections below) -- the same convention rmsynth already used
# for path/infileQ/infileU whenever match wasn't requested.

OUTDIR="$(cfg_get "${PIPELINE_CFG}" outdir)"
RUN_NAME="$(cfg_get "${PIPELINE_CFG}" run_name)"
[[ -z "${OUTDIR}" ]] && die "outdir= is required in ${PIPELINE_CFG}"
[[ -z "${RUN_NAME}" ]] && die "run_name= is required in ${PIPELINE_CFG}"
mkdir -p "${OUTDIR}"
OUTDIR="$(realpath "${OUTDIR}")"

# ============================================================
# Provenance directory -- shared across every run of this outdir/
# run_name (mkdir -p, never refused: it's a scaffolding directory, not
# an output). Each run's own records inside it are tagged with this
# run's own start timestamp so repeated runs never collide or
# overwrite each other's provenance -- the "refuse to overwrite" rule
# applies to actual DATA outputs (match_cubes/rm_synthesis/rmclean_cubes
# cubes, each already self-guarded by its own abort-if-exists check),
# not to this metadata directory.
# ============================================================
PROVDIR="${OUTDIR}/${RUN_NAME}.provenance"
mkdir -p "${PROVDIR}"
RUN_TAG="$(date +%Y%m%dT%H%M%S)"
COMMANDS_LOG="${PROVDIR}/commands.${RUN_TAG}.log"
: > "${COMMANDS_LOG}"
cp "${PIPELINE_CFG}" "${PROVDIR}/pipeline.${RUN_TAG}.cfg"

RUN_START="$(date -Iseconds)"
STAGE_TIMES=()

log "Config: ${PIPELINE_CFG}"
log "Stages: ${STAGES_RAW}"
log "Output: ${OUTDIR}/${RUN_NAME}.*"
log "Provenance: ${PROVDIR}"

RMSYNTH_INFILEQ=""
RMSYNTH_INFILEU=""

# ============================================================
# Stage: match (optional)
# ============================================================
if [[ "${DO_MATCH}" -eq 1 ]]; then
  log "=== match stage ==="
  MATCH_T0="$(date +%s)"
  MATCH_INPUT_PATH="$(cfg_get "${PIPELINE_CFG}" match_input_path)"
  MATCH_INFILEQ_RAW="$(cfg_get "${PIPELINE_CFG}" match_infileQ)"
  MATCH_INFILEU_RAW="$(cfg_get "${PIPELINE_CFG}" match_infileU)"
  MATCH_ARGS="$(cfg_get "${PIPELINE_CFG}" match_args)"
  [[ -z "${MATCH_INPUT_PATH}" ]] && die "match_input_path= is required when stages includes match"
  [[ -z "${MATCH_INFILEQ_RAW}" ]] && die "match_infileQ= is required when stages includes match"
  [[ -z "${MATCH_INFILEU_RAW}" ]] && die "match_infileU= is required when stages includes match"

  IFS=',' read -ra Q_BANDS <<< "${MATCH_INFILEQ_RAW}"
  IFS=',' read -ra U_BANDS <<< "${MATCH_INFILEU_RAW}"
  [[ "${#Q_BANDS[@]}" -eq "${#U_BANDS[@]}" ]] || die "match_infileQ and match_infileU must list the same number of bands"

  MATCH_EXE="${ROOT_DIR}/bin/match_cubes"
  [[ -x "${MATCH_EXE}" ]] || die "match_cubes not built -- run: make match_cubes"

  # match_cubes always writes its own output ALONGSIDE its input (no
  # output-directory option of its own) -- to get that output onto NVMe
  # (outdir, see this cfg's own comment) rather than match_input_path's
  # own disk, symlink each real input file into a fresh directory under
  # outdir and pass THOSE symlink paths as infiles=. Reads still follow
  # the symlink back to the real file wherever match_input_path
  # physically is (unchanged, unavoidable without copying the source
  # data itself); only the WRITE side (the new <symlink_name><outsuffix>
  # file) lands in the symlink's own directory, on NVMe. A pre-existing
  # symlink pointing at the SAME real file is reused across runs
  # (harmless -- it's just a directory entry, not a data output, so this
  # doesn't relax the abort-if-exists guarantee on the actual FITS
  # output next to it); one pointing elsewhere aborts rather than being
  # silently repointed.
  MATCH_SYMLINK_DIR="${OUTDIR}/match_input_symlinks"
  mkdir -p "${MATCH_SYMLINK_DIR}"

  ALL_INFILES=""
  for b in "${Q_BANDS[@]}" "${U_BANDS[@]}"; do
    f="${MATCH_INPUT_PATH%/}/${b}"
    [[ -f "${f}" ]] || die "match input not found: ${f}"
    link="${MATCH_SYMLINK_DIR}/${b}"
    if [[ -e "${link}" || -L "${link}" ]]; then
      [[ "$(readlink -f "${link}")" == "$(readlink -f "${f}")" ]] || \
        die "stale symlink at ${link} does not point at ${f} -- remove it first"
    else
      ln -s "${f}" "${link}"
    fi
    ALL_INFILES="${ALL_INFILES:+${ALL_INFILES},}${link}"
  done
  MATCH_INPUT_PATH="${MATCH_SYMLINK_DIR}"

  MATCH_MANIFEST="${PROVDIR}/match_manifest.${RUN_TAG}.txt"
  log "Running: bin/match_cubes ${MATCH_ARGS} infiles=<${#Q_BANDS[@]} Q + ${#U_BANDS[@]} U bands> manifest=${MATCH_MANIFEST}"
  # shellcheck disable=SC2086
  run_logged "${PROVDIR}/match.${RUN_TAG}.stdout.log" \
    "${MATCH_EXE}" ${MATCH_ARGS} "infiles=${ALL_INFILES}" "manifest=${MATCH_MANIFEST}"

  # Effective path per file comes from match_cubes' own manifest --
  # never guessed from filesystem state (see this script's own top
  # comment). One tab-separated line per input:
  # "<infile>\t(SKIPPED|PROCESSED)\t<effective_path>". match_cubes
  # always writes its own output ALONGSIDE its input (<infile><suffix>,
  # no separate output-directory option), so every effective path --
  # SKIPPED (the original input) or PROCESSED (<infile><outsuffix>) --
  # is guaranteed to still live directly under match_input_path. This
  # lets rmsynth's own cfg use the SAME path=/infileQ=/infileU=
  # convention every other cfg in this project already uses (a single
  # directory prefix + bare filenames) rather than inventing an
  # absolute-path special case rm_synthesis' own cfg parser doesn't
  # support (it hard-rejects a blank path= value).
  manifest_lookup() {
    local infile="$1"
    awk -F'\t' -v f="${infile}" '$1==f {print $3; exit}' "${MATCH_MANIFEST}"
  }

  RMSYNTH_INFILEQ=""
  for b in "${Q_BANDS[@]}"; do
    f="${MATCH_INPUT_PATH%/}/${b}"
    eff="$(manifest_lookup "${f}")"
    [[ -z "${eff}" ]] && die "match_cubes manifest has no entry for: ${f}"
    [[ "$(dirname "${eff}")" == "${MATCH_INPUT_PATH%/}" ]] || die "match_cubes manifest entry for ${f} (${eff}) is not directly under match_input_path (${MATCH_INPUT_PATH}) -- cannot express as a single path=+infileQ= pair"
    RMSYNTH_INFILEQ="${RMSYNTH_INFILEQ:+${RMSYNTH_INFILEQ},}$(basename "${eff}")"
  done
  RMSYNTH_INFILEU=""
  for b in "${U_BANDS[@]}"; do
    f="${MATCH_INPUT_PATH%/}/${b}"
    eff="$(manifest_lookup "${f}")"
    [[ -z "${eff}" ]] && die "match_cubes manifest has no entry for: ${f}"
    [[ "$(dirname "${eff}")" == "${MATCH_INPUT_PATH%/}" ]] || die "match_cubes manifest entry for ${f} (${eff}) is not directly under match_input_path (${MATCH_INPUT_PATH}) -- cannot express as a single path=+infileQ= pair"
    RMSYNTH_INFILEU="${RMSYNTH_INFILEU:+${RMSYNTH_INFILEU},}$(basename "${eff}")"
  done
  STAGE_TIMES+=("match:$(( $(date +%s) - MATCH_T0 ))s")
  log "match stage done."
fi

# ============================================================
# Stage: rmsynth (optional)
# ============================================================
if [[ "${DO_RMSYNTH}" -eq 1 ]]; then
  log "=== rmsynth stage ==="
  RMSYNTH_T0="$(date +%s)"
  RMSYNTH_TEMPLATE="$(cfg_get "${PIPELINE_CFG}" rmsynth_cfg_template)"
  RMSYNTH_BACKEND="$(cfg_get "${PIPELINE_CFG}" rmsynth_backend)"
  RMSYNTH_THREADS="$(cfg_get "${PIPELINE_CFG}" rmsynth_omp_threads)"
  [[ -z "${RMSYNTH_TEMPLATE}" ]] && die "rmsynth_cfg_template= is required when stages includes rmsynth"
  [[ -z "${RMSYNTH_BACKEND}" ]] && RMSYNTH_BACKEND="auto"
  [[ -z "${RMSYNTH_THREADS}" ]] && RMSYNTH_THREADS="6"
  if [[ ! -f "${RMSYNTH_TEMPLATE}" ]]; then
    [[ -f "${ROOT_DIR}/${RMSYNTH_TEMPLATE}" ]] && RMSYNTH_TEMPLATE="${ROOT_DIR}/${RMSYNTH_TEMPLATE}"
  fi
  [[ -f "${RMSYNTH_TEMPLATE}" ]] || die "rmsynth_cfg_template not found: ${RMSYNTH_TEMPLATE}"

  if [[ "${DO_RMCLEAN}" -eq 1 ]]; then
    WMO="$(cfg_get "${RMSYNTH_TEMPLATE}" write_mask_output)"
    # write_mask_output defaults to TRUE when absent (rm_synthesis_mod.f90's
    # own cfg%write_mask_output = .true. default) -- only an EXPLICIT
    # falsy value is a problem.
    if [[ -n "${WMO}" ]] && ! flag_is_true "${WMO}"; then
      die "rmsynth_cfg_template (${RMSYNTH_TEMPLATE}) has write_mask_output=${WMO}, but rmclean_cubes requires the MASK cube -- set write_mask_output=y (or remove the key, it defaults to true) in your template, or drop 'rmclean' from stages="
    fi
  fi

  RMSYNTH_CFG="${PROVDIR}/rmsynth.${RUN_TAG}.cfg"
  cp "${RMSYNTH_TEMPLATE}" "${RMSYNTH_CFG}"
  RMSYNTH_OUTFILE="${OUTDIR}/${RUN_NAME}"
  cfg_set_inplace "${RMSYNTH_CFG}" outfile "${RMSYNTH_OUTFILE}"
  if [[ "${DO_MATCH}" -eq 1 ]]; then
    cfg_set_inplace "${RMSYNTH_CFG}" path "${MATCH_INPUT_PATH%/}/"
    cfg_set_inplace "${RMSYNTH_CFG}" infileQ "${RMSYNTH_INFILEQ}"
    cfg_set_inplace "${RMSYNTH_CFG}" infileU "${RMSYNTH_INFILEU}"
  fi
  # else: match did not run in THIS invocation -- use the template's own
  # path=/infileQ=/infileU= unmodified (e.g. pointing at a previous,
  # separate match-only run's own real output).

  OUTPUT_MODE="$(cfg_get "${RMSYNTH_CFG}" output_mode)"
  AP_MODE="$(cfg_get "${RMSYNTH_CFG}" ap_angle_mode)"
  [[ -z "${OUTPUT_MODE}" ]] && OUTPUT_MODE="ap"
  [[ -z "${AP_MODE}" ]] && AP_MODE="phase"
  if [[ "${OUTPUT_MODE}" == "ri" ]]; then
    RMSYNTH_OUT1="${RMSYNTH_OUTFILE}.REAL.RMCUBE.FITS"
    RMSYNTH_OUT2="${RMSYNTH_OUTFILE}.IMAG.RMCUBE.FITS"
  else
    RMSYNTH_OUT1="${RMSYNTH_OUTFILE}.AMP.RMCUBE.FITS"
    RMSYNTH_OUT2="${RMSYNTH_OUTFILE}.PHA.RMCUBE.FITS"
    [[ "${AP_MODE}" == "pol" ]] && RMSYNTH_OUT2="${RMSYNTH_OUTFILE}.POLA.RMCUBE.FITS"
  fi

  # scratch/run_rmsynthesis_test.sh unconditionally `rm -f`'s these two
  # cubes before every run (a deliberate dev-test convenience for its own
  # direct, interactive use) -- check for pre-existing output HERE, before
  # calling it, so this pipeline's own "never silently delete/overwrite an
  # existing output" guarantee holds regardless of that runner's own
  # behaviour.
  if [[ -e "${RMSYNTH_OUT1}" || -e "${RMSYNTH_OUT2}" ]]; then
    die "rmsynth output already exists, refusing to proceed (would be deleted and regenerated): ${RMSYNTH_OUT1} / ${RMSYNTH_OUT2}"
  fi

  RUNNER="${ROOT_DIR}/scratch/run_rmsynthesis_test.sh"
  [[ -x "${RUNNER}" ]] || die "runner not found or not executable: ${RUNNER}"
  log "Running: ${RUNNER} ${RMSYNTH_CFG} ${RMSYNTH_THREADS} ${RMSYNTH_BACKEND}"
  run_logged "${PROVDIR}/rmsynth.${RUN_TAG}.stdout.log" \
    "${RUNNER}" "${RMSYNTH_CFG}" "${RMSYNTH_THREADS}" "${RMSYNTH_BACKEND}"

  [[ -s "${RMSYNTH_OUT1}" && -s "${RMSYNTH_OUT2}" ]] || die "rmsynth stage did not produce expected output cubes (${RMSYNTH_OUT1}, ${RMSYNTH_OUT2})"
  STAGE_TIMES+=("rmsynth:$(( $(date +%s) - RMSYNTH_T0 ))s")
  log "rmsynth stage done: ${RMSYNTH_OUT1}, ${RMSYNTH_OUT2}"
fi

# ============================================================
# Stage: rmclean (optional)
# ============================================================
if [[ "${DO_RMCLEAN}" -eq 1 ]]; then
  log "=== rmclean stage ==="
  RMCLEAN_T0="$(date +%s)"
  RMCLEAN_TEMPLATE="$(cfg_get "${PIPELINE_CFG}" rmclean_cfg_template)"
  [[ -z "${RMCLEAN_TEMPLATE}" ]] && die "rmclean_cfg_template= is required when stages includes rmclean"
  if [[ ! -f "${RMCLEAN_TEMPLATE}" ]]; then
    [[ -f "${ROOT_DIR}/${RMCLEAN_TEMPLATE}" ]] && RMCLEAN_TEMPLATE="${ROOT_DIR}/${RMCLEAN_TEMPLATE}"
  fi
  [[ -f "${RMCLEAN_TEMPLATE}" ]] || die "rmclean_cfg_template not found: ${RMCLEAN_TEMPLATE}"

  RMCLEAN_CFG="${PROVDIR}/rmclean.${RUN_TAG}.cfg"
  cp "${RMCLEAN_TEMPLATE}" "${RMCLEAN_CFG}"

  if [[ "${DO_RMSYNTH}" -eq 1 ]]; then
    # rmsynth ran in THIS invocation -- chain its own real output
    # straight in, overriding whatever the template says.
    MASK_FILE="$(cfg_get "${RMSYNTH_CFG}" mask_cube_file)"
    [[ -z "${MASK_FILE}" ]] && MASK_FILE="${RMSYNTH_OUTFILE}.MASK.CUBE.FITS"
    [[ -s "${MASK_FILE}" ]] || die "expected MASK cube missing: ${MASK_FILE} (is write_mask_output=y in ${RMSYNTH_TEMPLATE}?)"
    RMCLEAN_OUTFILE="${OUTDIR}/${RUN_NAME}_cleaned"
    cfg_set_inplace "${RMCLEAN_CFG}" ampfile "${RMSYNTH_OUT1}"
    cfg_set_inplace "${RMCLEAN_CFG}" phafile "${RMSYNTH_OUT2}"
    cfg_set_inplace "${RMCLEAN_CFG}" maskfile "${MASK_FILE}"
    cfg_set_inplace "${RMCLEAN_CFG}" outfile "${RMCLEAN_OUTFILE}"
  else
    # rmsynth did NOT run in this invocation -- use the template's own
    # ampfile=/phafile=/maskfile=/outfile= unmodified (e.g. pointing at
    # a previous, separate rmsynth-only run's own real output).
    AMP_FILE="$(cfg_get "${RMCLEAN_CFG}" ampfile)"
    PHA_FILE="$(cfg_get "${RMCLEAN_CFG}" phafile)"
    MASK_FILE="$(cfg_get "${RMCLEAN_CFG}" maskfile)"
    RMCLEAN_OUTFILE="$(cfg_get "${RMCLEAN_CFG}" outfile)"
    [[ -z "${AMP_FILE}" || -z "${PHA_FILE}" || -z "${MASK_FILE}" || -z "${RMCLEAN_OUTFILE}" ]] && \
      die "rmclean_cfg_template (${RMCLEAN_TEMPLATE}) must set ampfile=/phafile=/maskfile=/outfile= itself when rmsynth is not also in stages="
    [[ -s "${AMP_FILE}" ]] || die "rmclean_cfg_template's own ampfile does not exist: ${AMP_FILE}"
    [[ -s "${PHA_FILE}" ]] || die "rmclean_cfg_template's own phafile does not exist: ${PHA_FILE}"
    [[ -s "${MASK_FILE}" ]] || die "rmclean_cfg_template's own maskfile does not exist: ${MASK_FILE}"
  fi

  RMCLEAN_EXE="${ROOT_DIR}/bin/rmclean_cubes"
  [[ -x "${RMCLEAN_EXE}" ]] || die "rmclean_cubes not built -- run: make rmclean_cubes"

  log "Running: bin/rmclean_cubes --config ${RMCLEAN_CFG}"
  run_logged "${PROVDIR}/rmclean.${RUN_TAG}.stdout.log" \
    "${RMCLEAN_EXE}" --config "${RMCLEAN_CFG}"

  for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA RESTORED.AMP RESTORED.PHA; do
    f="${RMCLEAN_OUTFILE}.${suffix}.RMCUBE.FITS"
    [[ -s "${f}" ]] || die "rmclean stage did not produce expected output: ${f}"
  done
  STAGE_TIMES+=("rmclean:$(( $(date +%s) - RMCLEAN_T0 ))s")
  log "rmclean stage done: ${RMCLEAN_OUTFILE}.{CLEAN,RESID,RESTORED}.{AMP,PHA}.RMCUBE.FITS"
fi

# ============================================================
# run_summary.txt
# ============================================================
{
  echo "pipeline cfg: ${PIPELINE_CFG}"
  echo "started:      ${RUN_START}"
  echo "finished:     $(date -Iseconds)"
  echo "stages:       ${STAGES_RAW}"
  echo "stage times:  ${STAGE_TIMES[*]}"
  echo "status:       OK"
  echo ""
  echo "outputs:"
  [[ "${DO_MATCH}" -eq 1 ]] && echo "  match: manifest ${MATCH_MANIFEST} (per-file effective paths -- SKIPPED entries are the original input, PROCESSED entries <input><outsuffix>)"
  [[ "${DO_RMSYNTH}" -eq 1 ]] && echo "  rmsynth: ${RMSYNTH_OUT1}, ${RMSYNTH_OUT2}"
  if [[ "${DO_RMCLEAN}" -eq 1 ]]; then
    for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA RESTORED.AMP RESTORED.PHA; do
      echo "  rmclean: ${RMCLEAN_OUTFILE}.${suffix}.RMCUBE.FITS"
    done
  fi
} > "${PROVDIR}/run_summary.${RUN_TAG}.txt"

log "Provenance: ${PROVDIR}"
log "Pipeline finished."
