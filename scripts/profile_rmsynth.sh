#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/profile_rmsynth.sh --cfg <config> [options]

Required:
  --cfg <path>              Config file for rm_synthesis

Options:
  --out <dir>               Output directory (default: scratch/profiles/<timestamp>)
  --profile-bin-dir <dir>   Directory for dedicated profiling binaries (default: scratch/profiles/bin)
  --omp-threads <N>         OMP threads for OMP run (default: 6)
  --run-omp                 Run OMP profiling only
  --run-gpu                 Run GPU profiling only
  --vtune                   Use Intel VTune for OMP profiling (required when --run-omp is used)
  --gpu-offload-info        Enable LIBOMPTARGET_INFO=16 for GPU run logs
  --cleanup-existing        Remove existing science outputs for this profile prefix before run
  --keep-outputs            Keep generated science FITS outputs (default: delete after profiling stage)
  --skip-build              Skip building binaries
  -h, --help                Show this help

Examples:
  bash scripts/profile_rmsynth.sh --cfg cfg/rmsynth-casa.fullim.cfg
  bash scripts/profile_rmsynth.sh --cfg cfg/rmsynth-casa.fullim.cfg --omp-threads 8 --run-gpu
EOF
}

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

abs_path() {
  local p="$1"
  local d=""
  local b=""

  if [[ "$p" == /* ]]; then
    printf '%s\n' "$p"
    return 0
  fi

  d="$(dirname "$p")"
  b="$(basename "$p")"
  if cd "$d" >/dev/null 2>&1; then
    printf '%s/%s\n' "$(pwd -P)" "$b"
    return 0
  fi

  return 1
}

cfg_value() {
  local file="$1"
  local key="$2"
  awk -F= -v k="$key" '
    $1==k {
      v=$0
      sub(/^[^=]*=/, "", v)
      gsub(/^[ \t]+|[ \t]+$/, "", v)
      print v
      exit
    }
  ' "$file" 2>/dev/null || true
}

prepare_profile_cfg_and_outprefix() {
  local src_cfg="$1"
  local outdir="$2"
  local suffix="${3:-}"
  local prof_cfg="$outdir/profile_run${suffix}.cfg"
  local cfg_out=""
  local out_base=""
  local out_prefix=""

  cfg_out="$(cfg_value "$src_cfg" outfile)"
  if [[ -z "$cfg_out" ]]; then
    echo "ERROR: config does not define 'outfile=': $src_cfg" >&2
    exit 2
  fi

  out_base="$(basename "$cfg_out")"
  if [[ -z "$out_base" || "$out_base" == "." || "$out_base" == "/" ]]; then
    echo "ERROR: invalid outfile value in config: $cfg_out" >&2
    exit 2
  fi

  out_prefix="$outdir/$out_base$suffix"

  awk -v v="$out_prefix" '
    BEGIN{done=0}
    {
      if($0 ~ /^outfile=/ && done==0){
        print "outfile=" v
        done=1
      } else {
        print $0
      }
    }
    END{
      if(done==0) print "outfile=" v
    }
  ' "$src_cfg" >"$prof_cfg"

  printf '%s\n' "$prof_cfg|$out_prefix"
}

abort_if_output_exists() {
  local out_prefix="$1"
  local -a hits=()
  local f

  for f in \
    "$out_prefix.AMP.RMCUBE.FITS" \
    "$out_prefix.REAL.RMCUBE.FITS" \
    "$out_prefix.PHA.RMCUBE.FITS" \
    "$out_prefix.POLA.RMCUBE.FITS" \
    "$out_prefix.IMAG.RMCUBE.FITS" \
    "$out_prefix.NVALID.MAP.FITS" \
    "$out_prefix.MASK.CUBE.FITS"
  do
    if [[ -e "$f" ]]; then
      hits+=("$f")
    fi
  done

  if [[ ${#hits[@]} -gt 0 ]]; then
    echo "" >&2
    echo "ALARM: profiling output files already exist for this run prefix:" >&2
    printf '  %s\n' "${hits[@]}" >&2
    echo "" >&2
    echo "Refusing to run profiling to avoid mixing old/new science outputs." >&2
    echo "Clean up these files and rerun:" >&2
    printf '  rm -f %q*\n' "$out_prefix" >&2
    echo "" >&2
    exit 3
  fi
}

cleanup_existing_outputs() {
  local out_prefix="$1"
  local f

  for f in \
    "$out_prefix.AMP.RMCUBE.FITS" \
    "$out_prefix.REAL.RMCUBE.FITS" \
    "$out_prefix.PHA.RMCUBE.FITS" \
    "$out_prefix.POLA.RMCUBE.FITS" \
    "$out_prefix.IMAG.RMCUBE.FITS" \
    "$out_prefix.NVALID.MAP.FITS" \
    "$out_prefix.MASK.CUBE.FITS"
  do
    if [[ -e "$f" ]]; then
      rm -f "$f"
    fi
  done
}

resolve_tool() {
  local tool="$1"
  local p=""

  if p="$(command -v "$tool" 2>/dev/null)" && [[ -n "$p" ]]; then
    printf '%s\n' "$p"
    return 0
  fi

  case "$tool" in
    nsys)
      p="$(ls -1d /opt/nvidia/nsight-systems/*/bin/nsys 2>/dev/null | sort -V | tail -n 1 || true)"
      if [[ -n "$p" ]]; then
        printf '%s\n' "$p"
        return 0
      fi
      p="$(ls -1d /opt/nvidia/nsight-systems/*/target-linux-x64/nsys 2>/dev/null | sort -V | tail -n 1 || true)"
      ;;
    ncu)
      p="$(ls -1d /opt/nvidia/nsight-compute/*/ncu 2>/dev/null | sort -V | tail -n 1 || true)"
      if [[ -n "$p" ]]; then
        printf '%s\n' "$p"
        return 0
      fi
      p="$(ls -1d /opt/nvidia/nsight-compute/*/target/linux-desktop-glibc_2_11_3-x64/ncu 2>/dev/null | sort -V | tail -n 1 || true)"
      ;;
    vtune)
      # oneAPI default paths (Intel oneAPI toolkit)
      for d in \
        /opt/intel/oneapi/vtune/latest/bin64 \
        /opt/intel/vtune_profiler/bin64 \
        /opt/intel/oneapi/vtune/*/bin64
      do
        if [[ -x "$d/vtune" ]]; then
          printf '%s/vtune\n' "$d"
          return 0
        fi
      done

      # If VTune was installed but not exported into PATH yet, load oneAPI env
      # and try again. This is what a fresh shell needs after installation.
      if [[ -f /opt/intel/oneapi/setvars.sh ]]; then
        # shellcheck disable=SC1091
        source /opt/intel/oneapi/setvars.sh --force >/dev/null 2>&1 || true
        if p="$(command -v vtune 2>/dev/null)" && [[ -n "$p" ]]; then
          printf '%s\n' "$p"
          return 0
        fi
      fi
      ;;
    *)
      p=""
      ;;
  esac

  if [[ -n "$p" ]]; then
    printf '%s\n' "$p"
    return 0
  fi
  return 1
}

now_ts() {
  date +%Y%m%d_%H%M%S
}

log() {
  printf '[profile] %s\n' "$*"
}

start_sampler() {
  local label="$1"
  local outdir="$2"
  local pid=""

  case "$label" in
    iostat)
      if have_cmd iostat; then
        iostat -xz 1 >"$outdir/iostat.log" 2>&1 &
        pid=$!
      fi
      ;;
    nvidia)
      if have_cmd nvidia-smi; then
        nvidia-smi \
          --query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used,memory.total,power.draw,temperature.gpu,clocks.gr,clocks.mem \
          --format=csv -l 1 >"$outdir/nvidia_smi.log" 2>&1 &
        pid=$!
      fi
      ;;
    pidstat_omp)
      if have_cmd pidstat; then
        pidstat -durh -C rm_synthesis 1 >"$outdir/pidstat_omp.log" 2>&1 &
        pid=$!
      fi
      ;;
    pidstat_gpu)
      if have_cmd pidstat; then
        pidstat -durh -C rm_synthesis 1 >"$outdir/pidstat_gpu.log" 2>&1 &
        pid=$!
      fi
      ;;
  esac

  echo "$pid"
}

stop_sampler() {
  local pid="$1"
  if [[ -n "$pid" ]]; then
    kill "$pid" >/dev/null 2>&1 || true
    wait "$pid" >/dev/null 2>&1 || true
  fi
}

extract_time_summary() {
  local time_log="$1"
  local out_txt="$2"
  {
    grep -E 'Elapsed \(wall clock\) time|Percent of CPU this job got|Maximum resident set size|File system inputs|File system outputs' "$time_log" || true
  } >>"$out_txt"
}

get_time_value() {
  local time_log="$1"
  local key="$2"
  awk -v k="$key" '
    index($0, k) > 0 {
      p=index($0, k)
      line=substr($0, p + length(k))
      d=index(line, ": ")
      if(d>0) line=substr(line, d+2)
      else sub(/^[[:space:]]*:[[:space:]]*/, "", line)
      gsub(/^[ \t]+/, "", line)
      print line
      exit
    }
  ' "$time_log" 2>/dev/null || true
}

avg_pidstat_cpu() {
  local pid_log="$1"
  local field="$2"
  # field choices: usr sys iowait
  awk -v f="$field" '
    BEGIN {col=0}
    /Average:.*UID.*usr.*sys.*iowait/ {
      for(i=1;i<=NF;i++) {
        if($i=="usr" && f=="usr") col=i
        if($i=="sys" && f=="sys") col=i
        if($i=="iowait" && f=="iowait") col=i
      }
      next
    }
    /^Average:/ && col>0 {
      if($col ~ /^-?[0-9]+(\.[0-9]+)?$/) {s+=$col; n+=1}
    }
    END {if(n>0) printf("%.2f", s/n)}
  ' "$pid_log" 2>/dev/null || true
}

avg_pidstat_sample_field() {
  local pid_log="$1"
  local field="$2"
  # field choices from sampled rows: usr sys wait cpu
  awk -v f="$field" '
    $0 ~ /^#/ {next}
    NF < 9 {next}
    $NF !~ /^rm_synthesis/ {next}
    {
      if(f=="usr")  v=$4
      else if(f=="sys")  v=$5
      else if(f=="wait") v=$7
      else if(f=="cpu")  v=$8
      else v=""
      if(v ~ /^-?[0-9]+(\.[0-9]+)?$/) {s+=v; n+=1}
    }
    END {if(n>0) printf("%.2f", s/n)}
  ' "$pid_log" 2>/dev/null || true
}

pct_of() {
  local part="$1"
  local total="$2"
  awk -v p="$part" -v t="$total" 'BEGIN {if(t+0>0) printf("%.2f", (100.0*p)/t)}' 2>/dev/null || true
}

max_iostat_util() {
  local io_log="$1"
  awk '
    BEGIN {mx=0}
    /%util/ {next}
    NF>=6 {
      val=$(NF)
      if(val ~ /^-?[0-9]+(\.[0-9]+)?$/) {
        if(val+0 > mx) mx=val+0
      }
    }
    END {if(mx>0) printf("%.2f", mx)}
  ' "$io_log" 2>/dev/null || true
}

avg_iostat_await() {
  local io_log="$1"
  awk '
    BEGIN {col=0}
    /await/ {
      for(i=1;i<=NF;i++) if($i=="await") col=i
      next
    }
    col>0 && NF>=col {
      v=$col
      if(v ~ /^-?[0-9]+(\.[0-9]+)?$/){s+=v;n+=1}
    }
    END {if(n>0) printf("%.2f", s/n)}
  ' "$io_log" 2>/dev/null || true
}

extract_nsys_launch_count() {
  local nsys_stats="$1"
  # best-effort count from textual stats table mentions
  grep -Eic 'cudaLaunchKernel|Kernel Launch|\bkernel\b' "$nsys_stats" 2>/dev/null || true
}

extract_offload_events() {
  local run_log="$1"
  grep -Eic 'Libomptarget|libomptarget|TARGET|omptarget' "$run_log" 2>/dev/null || true
}

is_perf_report_valid() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  grep -Eqi 'SKIPPED:|zero-sized data|No permission to enable|failed to open|Permission denied' "$f" && return 1
  return 0
}

is_nsys_stats_valid() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  # If CUDA summaries are present, treat as valid even with importer warnings.
  if grep -q 'CUDA API Summary (cuda_api_sum)' "$f" 2>/dev/null ||
     grep -q 'CUDA GPU Kernel Summary (cuda_gpu_kern_sum)' "$f" 2>/dev/null; then
    return 0
  fi
  grep -Eqi 'SKIPPED:|usage: nsys stats|Try .nsys stats --help.|No such file|ERROR' "$f" && return 1
  return 0
}

nsys_section_pct_sum() {
  local stats="$1"
  local section_pat="$2"
  local stop_pat="$3"
  local name_regex="$4"
  awk -v sec="$section_pat" -v stop="$stop_pat" -v r="$name_regex" '
    BEGIN {in=0; s=0}
    $0 ~ sec {in=1; next}
    in && $0 ~ stop {in=0}
    in {
      if($1 ~ /^[0-9]+(\.[0-9]+)?$/) {
        name=$NF
        if(name ~ r) s+=$1
      }
    }
    END {if(s>0) printf("%.2f", s)}
  ' "$stats" 2>/dev/null || true
}

is_ncu_csv_valid() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  grep -Eqi 'SKIPPED:|==ERROR==|No metrics to collect' "$f" && return 1
  return 0
}

extract_ncu_pct_metric() {
  local csv="$1"
  local metric_regex="$2"
  awk -F, -v r="$metric_regex" '
    BEGIN {IGNORECASE=1}
    $0 ~ r {
      v=$NF
      gsub(/"/,"",v)
      gsub(/%/,"",v)
      gsub(/^[ \t]+|[ \t]+$/, "", v)
      if(v ~ /^-?[0-9]+(\.[0-9]+)?$/) {print v; exit}
    }
  ' "$csv" 2>/dev/null || true
}

avg_nvidia_col() {
  local nvsmi_log="$1"
  local col="$2"
  awk -F, -v c="$col" 'NR>1 {
    gsub(/^[ \t]+/,"",$c)
    gsub(/ %/,"",$c)
    gsub(/ W/,"",$c)
    gsub(/ MiB/,"",$c)
    if($c ~ /^-?[0-9]+(\.[0-9]+)?$/){s+=$c;n+=1}
  }
  END {if(n>0) printf("%.2f", s/n)}' "$nvsmi_log" 2>/dev/null || true
}

generate_consolidated_report() {
  local outdir="$1"
  local report_md="$outdir/consolidated_report.md"
  local report_html="$outdir/consolidated_report.html"

  local om_elapsed="N/A"
  local om_cpu="N/A"
  local om_rss="N/A"
  local om_fsin="N/A"
  local om_fsout="N/A"

  local gp_elapsed="N/A"
  local gp_cpu="N/A"
  local gp_rss="N/A"
  local gp_fsin="N/A"
  local gp_fsout="N/A"
  local gp_util="N/A"
  local gp_mem_util="N/A"
  local gp_mem_used="N/A"
  local gp_power="N/A"
  local gp_temp="N/A"
  local om_usr="N/A"
  local om_sys="N/A"
  local om_iow="N/A"
  local gp_usr="N/A"
  local gp_sys="N/A"
  local gp_iow="N/A"
  local io_util_max="N/A"
  local io_await_avg="N/A"
  local nsys_launch_count="N/A"
  local offload_evt_count="N/A"
  local ncu_sm_pct="N/A"
  local ncu_mem_pct="N/A"
  local om_usr_s=""
  local om_sys_s=""
  local om_wait_s=""
  local om_cpu_s=""
  local gp_usr_s=""
  local gp_sys_s=""
  local gp_wait_s=""
  local gp_cpu_s=""
  local om_compute_pct="N/A"
  local om_iowait_pct="N/A"
  local om_other_pct="N/A"
  local gp_host_compute_pct="N/A"
  local gp_host_iowait_pct="N/A"
  local gp_host_other_pct="N/A"
  local gp_cuda_sync_pct="N/A"
  local gp_cuda_xfer_pct="N/A"
  local gp_cuda_launch_pct="N/A"
  local gp_cuda_other_pct="N/A"

  if [[ -f "$outdir/time_omp.log" ]]; then
    om_elapsed=$(get_time_value "$outdir/time_omp.log" 'Elapsed (wall clock) time')
    om_cpu=$(get_time_value "$outdir/time_omp.log" 'Percent of CPU this job got')
    om_rss=$(get_time_value "$outdir/time_omp.log" 'Maximum resident set size')
    om_fsin=$(get_time_value "$outdir/time_omp.log" 'File system inputs')
    om_fsout=$(get_time_value "$outdir/time_omp.log" 'File system outputs')
  fi

  if [[ -f "$outdir/time_gpu.log" ]]; then
    gp_elapsed=$(get_time_value "$outdir/time_gpu.log" 'Elapsed (wall clock) time')
    gp_cpu=$(get_time_value "$outdir/time_gpu.log" 'Percent of CPU this job got')
    gp_rss=$(get_time_value "$outdir/time_gpu.log" 'Maximum resident set size')
    gp_fsin=$(get_time_value "$outdir/time_gpu.log" 'File system inputs')
    gp_fsout=$(get_time_value "$outdir/time_gpu.log" 'File system outputs')
  fi

  if [[ -f "$outdir/nvidia_smi.log" ]]; then
    gp_util=$(avg_nvidia_col "$outdir/nvidia_smi.log" 2)
    gp_mem_util=$(avg_nvidia_col "$outdir/nvidia_smi.log" 3)
    gp_mem_used=$(avg_nvidia_col "$outdir/nvidia_smi.log" 4)
    gp_power=$(avg_nvidia_col "$outdir/nvidia_smi.log" 6)
    gp_temp=$(avg_nvidia_col "$outdir/nvidia_smi.log" 7)
  fi

  if [[ -f "$outdir/pidstat_omp.log" ]]; then
    om_usr=$(avg_pidstat_cpu "$outdir/pidstat_omp.log" usr)
    om_sys=$(avg_pidstat_cpu "$outdir/pidstat_omp.log" sys)
    om_iow=$(avg_pidstat_cpu "$outdir/pidstat_omp.log" iowait)
    om_usr_s=$(avg_pidstat_sample_field "$outdir/pidstat_omp.log" usr)
    om_sys_s=$(avg_pidstat_sample_field "$outdir/pidstat_omp.log" sys)
    om_wait_s=$(avg_pidstat_sample_field "$outdir/pidstat_omp.log" wait)
    om_cpu_s=$(avg_pidstat_sample_field "$outdir/pidstat_omp.log" cpu)
  fi

  if [[ -f "$outdir/pidstat_gpu.log" ]]; then
    gp_usr=$(avg_pidstat_cpu "$outdir/pidstat_gpu.log" usr)
    gp_sys=$(avg_pidstat_cpu "$outdir/pidstat_gpu.log" sys)
    gp_iow=$(avg_pidstat_cpu "$outdir/pidstat_gpu.log" iowait)
    gp_usr_s=$(avg_pidstat_sample_field "$outdir/pidstat_gpu.log" usr)
    gp_sys_s=$(avg_pidstat_sample_field "$outdir/pidstat_gpu.log" sys)
    gp_wait_s=$(avg_pidstat_sample_field "$outdir/pidstat_gpu.log" wait)
    gp_cpu_s=$(avg_pidstat_sample_field "$outdir/pidstat_gpu.log" cpu)
  fi

  if [[ -f "$outdir/iostat.log" ]]; then
    io_util_max=$(max_iostat_util "$outdir/iostat.log")
    io_await_avg=$(avg_iostat_await "$outdir/iostat.log")
  fi

  [[ -z "$om_usr" ]] && om_usr="N/A"
  [[ -z "$om_sys" ]] && om_sys="N/A"
  [[ -z "$om_iow" ]] && om_iow="N/A"
  [[ -z "$gp_usr" ]] && gp_usr="N/A"
  [[ -z "$gp_sys" ]] && gp_sys="N/A"
  [[ -z "$gp_iow" ]] && gp_iow="N/A"
  [[ -z "$io_util_max" ]] && io_util_max="N/A"
  [[ -z "$io_await_avg" ]] && io_await_avg="N/A"

  if [[ -f "$outdir/nsys_gpu_stats.txt" ]]; then
    nsys_launch_count=$(extract_nsys_launch_count "$outdir/nsys_gpu_stats.txt")
    gp_cuda_sync_pct=$(nsys_section_pct_sum "$outdir/nsys_gpu_stats.txt" 'CUDA API Summary \(cuda_api_sum\)' '^Processing ' 'cuCtxSynchronize')
    gp_cuda_xfer_pct=$(nsys_section_pct_sum "$outdir/nsys_gpu_stats.txt" 'CUDA API Summary \(cuda_api_sum\)' '^Processing ' 'cuMemcpy')
    gp_cuda_launch_pct=$(nsys_section_pct_sum "$outdir/nsys_gpu_stats.txt" 'CUDA API Summary \(cuda_api_sum\)' '^Processing ' 'cuLaunchKernel')
  fi

  if [[ -f "$outdir/run_gpu.log" ]]; then
    offload_evt_count=$(extract_offload_events "$outdir/run_gpu.log")
  fi

  if is_ncu_csv_valid "$outdir/ncu_gpu.csv"; then
    ncu_sm_pct=$(extract_ncu_pct_metric "$outdir/ncu_gpu.csv" 'sm__throughput.avg.pct_of_peak_sustained_elapsed|SM Throughput')
    ncu_mem_pct=$(extract_ncu_pct_metric "$outdir/ncu_gpu.csv" 'dram__throughput.avg.pct_of_peak_sustained_elapsed|Memory Throughput')
  fi

  [[ -z "$ncu_sm_pct" ]] && ncu_sm_pct="N/A"
  [[ -z "$ncu_mem_pct" ]] && ncu_mem_pct="N/A"

  if [[ -n "$om_cpu_s" && "$om_cpu_s" != "0.00" ]]; then
    om_compute_pct=$(pct_of "${om_usr_s:-0}" "$om_cpu_s")
    om_iowait_pct=$(pct_of "${om_wait_s:-0}" "$om_cpu_s")
    tmp_other=$(awk -v c="${om_compute_pct:-0}" -v i="${om_iowait_pct:-0}" 'BEGIN{v=100-c-i; if(v<0)v=0; printf("%.2f", v)}')
    om_other_pct="$tmp_other"
  fi

  if [[ -n "$gp_cpu_s" && "$gp_cpu_s" != "0.00" ]]; then
    gp_host_compute_pct=$(pct_of "${gp_usr_s:-0}" "$gp_cpu_s")
    gp_host_iowait_pct=$(pct_of "${gp_wait_s:-0}" "$gp_cpu_s")
    tmp_other=$(awk -v c="${gp_host_compute_pct:-0}" -v i="${gp_host_iowait_pct:-0}" 'BEGIN{v=100-c-i; if(v<0)v=0; printf("%.2f", v)}')
    gp_host_other_pct="$tmp_other"
  fi

  if [[ "$gp_cuda_sync_pct" != "N/A" || "$gp_cuda_xfer_pct" != "N/A" || "$gp_cuda_launch_pct" != "N/A" ]]; then
    c_sync=${gp_cuda_sync_pct:-0}
    c_xfer=${gp_cuda_xfer_pct:-0}
    c_launch=${gp_cuda_launch_pct:-0}
    gp_cuda_other_pct=$(awk -v s="$c_sync" -v x="$c_xfer" -v l="$c_launch" 'BEGIN{v=100-s-x-l; if(v<0)v=0; printf("%.2f", v)}')
  fi

  [[ -z "$gp_cuda_sync_pct" ]] && gp_cuda_sync_pct="N/A"
  [[ -z "$gp_cuda_xfer_pct" ]] && gp_cuda_xfer_pct="N/A"
  [[ -z "$gp_cuda_launch_pct" ]] && gp_cuda_launch_pct="N/A"
  [[ -z "$gp_cuda_other_pct" ]] && gp_cuda_other_pct="N/A"

  {
    echo "# RM-Synthesis Profiling Report"
    echo
    echo "Generated: $(date -Is)"
    echo
    if [[ -f "$outdir/metadata.txt" ]]; then
      echo "## Run Metadata"
      echo
      sed 's/^/- /' "$outdir/metadata.txt"
      echo
    fi

    echo "## Executive Summary"
    echo
    if [[ -f "$outdir/classify_report.txt" ]]; then
      sed '1,2d' "$outdir/classify_report.txt"
    else
      echo "No classification report available."
    fi
    echo

    echo "## OMP vs GPU Key Metrics"
    echo
    echo "| Metric | OMP | GPU |"
    echo "|---|---:|---:|"
    echo "| Wall time | ${om_elapsed} | ${gp_elapsed} |"
    echo "| CPU utilization | ${om_cpu} | ${gp_cpu} |"
    echo "| Max RSS (kB) | ${om_rss} | ${gp_rss} |"
    echo "| Filesystem inputs | ${om_fsin} | ${gp_fsin} |"
    echo "| Filesystem outputs | ${om_fsout} | ${gp_fsout} |"
    echo "| Avg GPU util (%) | N/A | ${gp_util} |"
    echo "| Avg GPU mem util (%) | N/A | ${gp_mem_util} |"
    echo "| Avg GPU mem used (MiB) | N/A | ${gp_mem_used} |"
    echo "| Avg GPU power (W) | N/A | ${gp_power} |"
    echo "| Avg GPU temp (C) | N/A | ${gp_temp} |"
    echo "| OMP avg usr/sys/iowait (%) | ${om_usr}/${om_sys}/${om_iow} | N/A |"
    echo "| GPU avg usr/sys/iowait (%) | N/A | ${gp_usr}/${gp_sys}/${gp_iow} |"
    echo "| iostat max %util | ${io_util_max} | ${io_util_max} |"
    echo "| iostat avg await (ms) | ${io_await_avg} | ${io_await_avg} |"
    echo "| nsys kernel-launch indicator count | N/A | ${nsys_launch_count} |"
    echo "| Offload runtime event count (log) | N/A | ${offload_evt_count} |"
    echo "| NCU SM throughput (% peak) | N/A | ${ncu_sm_pct} |"
    echo "| NCU DRAM throughput (% peak) | N/A | ${ncu_mem_pct} |"
    echo

    echo "## Activity Breakdown (Best-Effort)"
    echo
    echo "All percentages below are derived from profiler sub-views, not a single unified timeline denominator."
    echo "Use these as directional signals to decide what to optimise next, not as strict wall-time totals."
    echo
    echo "| Activity Metric | OMP | GPU | Plain meaning |"
    echo "|---|---:|---:|---|"
    echo "| Host compute share (% of process CPU time) | ${om_compute_pct} | ${gp_host_compute_pct} | How much CPU effort was actual math/work. Higher = good. |"
    echo "| Host iowait share (% of process CPU time) | ${om_iowait_pct} | ${gp_host_iowait_pct} | CPU time stalled waiting on storage. High = disk I/O is hurting. |"
    echo "| Host other/system share (% of process CPU time) | ${om_other_pct} | ${gp_host_other_pct} | Runtime overhead, sync, driver calls, thread mgmt. High = orchestration cost. |"
    echo "| CUDA API sync share (% of CUDA API time) | N/A | ${gp_cuda_sync_pct} | GPU API time spent waiting for device to finish. High = lots of stall/sync. |"
    echo "| CUDA API transfer share (% of CUDA API time) | N/A | ${gp_cuda_xfer_pct} | GPU API time moving data between host and device. High = transfer overhead. |"
    echo "| CUDA API kernel-launch share (% of CUDA API time) | N/A | ${gp_cuda_launch_pct} | GPU API time spent just launching kernels. High = many small/fragmented launches. |"
    echo "| CUDA API other share (% of CUDA API time) | N/A | ${gp_cuda_other_pct} | Remaining GPU API activities not in the three above. |"
    echo

    echo "## Tool Availability and Artifacts"
    echo
    echo "| Tool/Artifact | Status | File |"
    echo "|---|---|---|"
    if is_perf_report_valid "$outdir/perf_omp.txt"; then
      echo "| perf (OMP) | Present | perf_omp.txt |"
    else
      echo "| perf (OMP) | Missing/Skipped/Invalid | perf_omp.txt |"
    fi
    if is_nsys_stats_valid "$outdir/nsys_gpu_stats.txt"; then
      echo "| nsys (GPU) | Present | nsys_gpu_stats.txt |"
    else
      echo "| nsys (GPU) | Missing/Skipped/Invalid | nsys_gpu_stats.txt |"
    fi
    if is_ncu_csv_valid "$outdir/ncu_gpu.csv"; then
      echo "| ncu (GPU) | Present | ncu_gpu.csv |"
    else
      echo "| ncu (GPU) | Missing/Skipped/Invalid | ncu_gpu.csv |"
    fi
    [[ -f "$outdir/pidstat_omp.log" ]] && echo "| pidstat OMP | Present | pidstat_omp.log |" || echo "| pidstat OMP | Missing/Skipped | pidstat_omp.log |"
    [[ -f "$outdir/pidstat_gpu.log" ]] && echo "| pidstat GPU | Present | pidstat_gpu.log |" || echo "| pidstat GPU | Missing/Skipped | pidstat_gpu.log |"
    [[ -f "$outdir/iostat.log" ]] && echo "| iostat | Present | iostat.log |" || echo "| iostat | Missing/Skipped | iostat.log |"
    [[ -f "$outdir/nvidia_smi.log" ]] && echo "| nvidia-smi sampler | Present | nvidia_smi.log |" || echo "| nvidia-smi sampler | Missing/Skipped | nvidia_smi.log |"
    echo

    echo "## Missing Tool Hints"
    echo
    if [[ ! -f "$outdir/nsys_gpu_stats.txt" ]] || grep -q 'SKIPPED: nsys not found' "$outdir/nsys_gpu_stats.txt" 2>/dev/null; then
      echo "- Nsight Systems (nsys) missing: install from NVIDIA CUDA toolkit/repo to get kernel launch timeline."
    fi
    if [[ -f "$outdir/nsys_gpu_stats.txt" ]] && ! is_nsys_stats_valid "$outdir/nsys_gpu_stats.txt"; then
      echo "- nsys stats output is invalid (parse/usage error). Try rerun; script now uses --force-export=true to refresh sqlite export."
    fi
    if [[ ! -f "$outdir/ncu_gpu.csv" ]] || grep -q 'SKIPPED: ncu not found' "$outdir/ncu_gpu.csv" 2>/dev/null; then
      echo "- Nsight Compute (ncu) missing: install from NVIDIA CUDA toolkit/repo for occupancy/memory-bound metrics."
    fi
    if [[ -f "$outdir/ncu_gpu.csv" ]] && ! is_ncu_csv_valid "$outdir/ncu_gpu.csv"; then
      echo "- ncu output invalid (target app error or no metrics). Check ncu_run.log and GPU command success."
    fi
    if [[ -f "$outdir/perf_omp.txt" ]] && ! is_perf_report_valid "$outdir/perf_omp.txt"; then
      echo "- perf output invalid/empty. Common causes: perf_event_paranoid permissions or short/failed OMP run."
    fi
    if ! command -v pandoc >/dev/null 2>&1; then
      echo "- pandoc missing: install with apt to generate HTML report from markdown."
    fi
    echo

    if [[ -f "$outdir/perf_omp.txt" ]]; then
      echo "## Top OMP perf symbols (first lines)"
      echo
      echo '```text'
      sed -n '1,40p' "$outdir/perf_omp.txt"
      echo '```'
      echo
    fi

    if [[ -f "$outdir/nsys_gpu_stats.txt" ]]; then
      echo "## nsys GPU stats excerpt"
      echo
      echo '```text'
      sed -n '1,80p' "$outdir/nsys_gpu_stats.txt"
      echo '```'
      echo
    fi

    echo "## Raw Logs"
    echo
    for f in metadata.txt classify_report.txt time_omp.log time_gpu.log run_omp.log run_gpu.log \
             perf_omp.txt nsys_gpu_stats.txt ncu_gpu.csv pidstat_omp.log pidstat_gpu.log iostat.log nvidia_smi.log; do
      if [[ -f "$outdir/$f" ]]; then
        echo "- $f"
      fi
    done
  } >"$report_md"

  if have_cmd pandoc; then
    pandoc "$report_md" -o "$report_html" >/dev/null 2>&1 || true
  fi
}

classify_omp() {
  local outdir="$1"
  local rep="$outdir/classify_report.txt"
  local time_log="$outdir/time_omp.log"

  local cpu_pct=""
  cpu_pct=$(awk -F: '/Percent of CPU this job got/ {gsub(/^[ \t]+/,"",$2); gsub(/%/,"",$2); print int($2)}' "$time_log" 2>/dev/null || true)

  local fs_in=""
  local fs_out=""
  fs_in=$(awk -F: '/File system inputs/ {gsub(/^[ \t]+/,"",$2); print int($2)}' "$time_log" 2>/dev/null || true)
  fs_out=$(awk -F: '/File system outputs/ {gsub(/^[ \t]+/,"",$2); print int($2)}' "$time_log" 2>/dev/null || true)

  echo "" >>"$rep"
  echo "[OMP Classification]" >>"$rep"

  if [[ -n "$cpu_pct" && "$cpu_pct" -ge 450 ]]; then
    echo "Primary: compute-limited (high aggregate CPU utilization)." >>"$rep"
  else
    echo "Primary: mixed/under-utilized CPU (CPU% not near expected thread aggregate)." >>"$rep"
  fi

  if [[ -n "$fs_in" && -n "$fs_out" ]]; then
    if [[ "$fs_in" -gt 5000000 || "$fs_out" -gt 5000000 ]]; then
      echo "Secondary: possible disk-I/O pressure (high filesystem input/output counters)." >>"$rep"
    else
      echo "Secondary: disk-I/O likely not dominant (modest filesystem I/O counters)." >>"$rep"
    fi
  fi

  echo "Evidence:" >>"$rep"
  extract_time_summary "$time_log" "$rep"

  if [[ -f "$outdir/perf_omp.txt" ]]; then
    echo "Top perf symbols:" >>"$rep"
    sed -n '1,30p' "$outdir/perf_omp.txt" >>"$rep" || true
  fi
}

classify_gpu() {
  local outdir="$1"
  local rep="$outdir/classify_report.txt"
  local time_log="$outdir/time_gpu.log"

  local gpu_util_avg=""
  if [[ -f "$outdir/nvidia_smi.log" ]]; then
    gpu_util_avg=$(awk -F, 'NR>1 {gsub(/ %/,"",$2); gsub(/^[ \t]+/,"",$2); s+=$2; n+=1} END {if(n>0) printf("%.1f",s/n)}' "$outdir/nvidia_smi.log" 2>/dev/null || true)
  fi

  local memcpy_hits=0
  local kernel_hits=0
  if [[ -f "$outdir/nsys_gpu_stats.txt" ]]; then
    memcpy_hits=$(grep -Eic 'memcpy|cudaMemcpy|\bHtoD\b|\bDtoH\b' "$outdir/nsys_gpu_stats.txt" || true)
    kernel_hits=$(grep -Eic 'kernel|cudaLaunchKernel|gpu__time_duration' "$outdir/nsys_gpu_stats.txt" || true)
  fi

  echo "" >>"$rep"
  echo "[GPU Classification]" >>"$rep"

  if [[ -n "$gpu_util_avg" ]]; then
    if awk "BEGIN{exit !($gpu_util_avg >= 85.0)}"; then
      echo "Primary: compute-active GPU (high average GPU utilization)." >>"$rep"
    else
      echo "Primary: potential idle/gap-limited GPU (average utilization not high)." >>"$rep"
    fi
  else
    echo "Primary: unknown (no nvidia-smi samples)." >>"$rep"
  fi

  if [[ -f "$outdir/nsys_gpu_stats.txt" ]]; then
    if [[ "$memcpy_hits" -gt "$kernel_hits" ]]; then
      echo "Secondary: transfer/sync-overhead likely significant (nsys text mentions memcpy/sync frequently)." >>"$rep"
    else
      echo "Secondary: kernel execution likely dominates over transfer mentions in nsys text." >>"$rep"
    fi
  else
    echo "Secondary: nsys unavailable, cannot separate kernel vs transfer timeline precisely." >>"$rep"
  fi

  echo "Evidence:" >>"$rep"
  extract_time_summary "$time_log" "$rep"

  if [[ -n "$gpu_util_avg" ]]; then
    echo "Average GPU utilization from nvidia-smi samples: ${gpu_util_avg}%" >>"$rep"
  fi

  if [[ -f "$outdir/ncu_gpu.csv" ]]; then
    if is_ncu_csv_valid "$outdir/ncu_gpu.csv"; then
      echo "ncu collected: see ncu_gpu.csv for memory-vs-compute metrics." >>"$rep"
    else
      echo "ncu output invalid: check ncu_run.log and ncu_gpu.csv for target-app errors." >>"$rep"
    fi
  fi
}

run_with_samplers() {
  local label="$1"
  local cmd="$2"
  local outdir="$3"

  local pid_iostat=""
  local pid_pidstat=""
  local pid_nvidia=""

  pid_iostat=$(start_sampler iostat "$outdir")

  if [[ "$label" == "omp" ]]; then
    pid_pidstat=$(start_sampler pidstat_omp "$outdir")
  else
    pid_pidstat=$(start_sampler pidstat_gpu "$outdir")
    pid_nvidia=$(start_sampler nvidia "$outdir")
  fi

  if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v -o "$outdir/time_${label}.log" bash -lc "$cmd" >"$outdir/run_${label}.log" 2>&1
  else
    bash -lc "$cmd" >"$outdir/run_${label}.log" 2>&1
    echo "NOTE: /usr/bin/time not found; detailed timing unavailable" >"$outdir/time_${label}.log"
  fi

  stop_sampler "$pid_nvidia"
  stop_sampler "$pid_pidstat"
  stop_sampler "$pid_iostat"
}

CFG=""
OUT_DIR=""
PROFILE_BIN_DIR=""
OMP_THREADS=6
RUN_OMP=1
RUN_GPU=1
VTUNE_REQUESTED=0
SKIP_BUILD=0
GPU_OFFLOAD_INFO=0
CLEANUP_EXISTING=0
KEEP_OUTPUTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cfg)
      CFG="${2:-}"
      shift 2
      ;;
    --out)
      OUT_DIR="${2:-}"
      shift 2
      ;;
    --profile-bin-dir)
      PROFILE_BIN_DIR="${2:-}"
      shift 2
      ;;
    --omp-threads)
      OMP_THREADS="${2:-}"
      shift 2
      ;;
    --run-omp)
      RUN_GPU=0
      RUN_OMP=1
      shift
      ;;
    --run-gpu)
      RUN_OMP=0
      RUN_GPU=1
      shift
      ;;
    --vtune)
      VTUNE_REQUESTED=1
      shift
      ;;
    --skip-build)
      SKIP_BUILD=1
      shift
      ;;
    --gpu-offload-info)
      GPU_OFFLOAD_INFO=1
      shift
      ;;
    --cleanup-existing)
      CLEANUP_EXISTING=1
      shift
      ;;
    --keep-outputs)
      KEEP_OUTPUTS=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 2
      ;;
  esac
done

if [[ -z "$CFG" ]]; then
  echo "ERROR: --cfg is required"
  usage
  exit 2
fi

if [[ ! -f "$CFG" ]]; then
  echo "ERROR: config file not found: $CFG"
  exit 2
fi

INVOKE_CWD="$(pwd -P)"
CFG="$(abs_path "$CFG")"
if [[ -z "$CFG" || ! -f "$CFG" ]]; then
  echo "ERROR: could not resolve config file path"
  exit 2
fi

if [[ -n "$OUT_DIR" && "$OUT_DIR" != /* ]]; then
  OUT_DIR="$INVOKE_CWD/$OUT_DIR"
fi

if [[ -n "$PROFILE_BIN_DIR" && "$PROFILE_BIN_DIR" != /* ]]; then
  PROFILE_BIN_DIR="$INVOKE_CWD/$PROFILE_BIN_DIR"
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

NSYS_BIN="$(resolve_tool nsys || true)"
NCU_BIN="$(resolve_tool ncu || true)"
VTUNE_BIN="$(resolve_tool vtune || true)"

if [[ "$RUN_OMP" -eq 1 && -z "$VTUNE_BIN" ]]; then
  echo "ERROR: VTune is required for OMP profiling but was not found on PATH or in /opt/intel/oneapi" >&2
  exit 1
fi

if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="$REPO_ROOT/scratch/profiles/$(now_ts)"
fi

if [[ -z "$PROFILE_BIN_DIR" ]]; then
  PROFILE_BIN_DIR="$REPO_ROOT/scratch/profiles/bin"
fi

if [[ "$VTUNE_REQUESTED" -eq 1 ]]; then
  log "VTune requested explicitly; OMP profiling will run through VTune"
fi

mkdir -p "$OUT_DIR"
mkdir -p "$PROFILE_BIN_DIR"

PROFILE_CFG_OMP=""
PROFILE_CFG_GPU=""
PROFILE_CFG_GPU_NSYS=""
PROFILE_CFG_GPU_NCU=""
PROFILE_OUT_PREFIX_OMP=""
PROFILE_OUT_PREFIX_GPU=""
PROFILE_OUT_PREFIX_GPU_NSYS=""
PROFILE_OUT_PREFIX_GPU_NCU=""

if [[ "$RUN_OMP" -eq 1 ]]; then
  cfg_prep_omp="$(prepare_profile_cfg_and_outprefix "$CFG" "$OUT_DIR" "_omp")"
  PROFILE_CFG_OMP="${cfg_prep_omp%%|*}"
  PROFILE_OUT_PREFIX_OMP="${cfg_prep_omp##*|}"
fi

if [[ "$RUN_GPU" -eq 1 ]]; then
  cfg_prep_gpu="$(prepare_profile_cfg_and_outprefix "$CFG" "$OUT_DIR" "_gpu")"
  PROFILE_CFG_GPU="${cfg_prep_gpu%%|*}"
  PROFILE_OUT_PREFIX_GPU="${cfg_prep_gpu##*|}"

  cfg_prep_gpu_nsys="$(prepare_profile_cfg_and_outprefix "$CFG" "$OUT_DIR" "_gpu_nsys")"
  PROFILE_CFG_GPU_NSYS="${cfg_prep_gpu_nsys%%|*}"
  PROFILE_OUT_PREFIX_GPU_NSYS="${cfg_prep_gpu_nsys##*|}"

  cfg_prep_gpu_ncu="$(prepare_profile_cfg_and_outprefix "$CFG" "$OUT_DIR" "_gpu_ncu")"
  PROFILE_CFG_GPU_NCU="${cfg_prep_gpu_ncu%%|*}"
  PROFILE_OUT_PREFIX_GPU_NCU="${cfg_prep_gpu_ncu##*|}"
fi

if [[ "$CLEANUP_EXISTING" -eq 1 ]]; then
  log "Cleaning existing science outputs for profile prefix(es)"
  [[ -n "$PROFILE_OUT_PREFIX_OMP" ]] && cleanup_existing_outputs "$PROFILE_OUT_PREFIX_OMP"
  [[ -n "$PROFILE_OUT_PREFIX_GPU" ]] && cleanup_existing_outputs "$PROFILE_OUT_PREFIX_GPU"
  [[ -n "$PROFILE_OUT_PREFIX_GPU_NSYS" ]] && cleanup_existing_outputs "$PROFILE_OUT_PREFIX_GPU_NSYS"
  [[ -n "$PROFILE_OUT_PREFIX_GPU_NCU" ]] && cleanup_existing_outputs "$PROFILE_OUT_PREFIX_GPU_NCU"
fi

[[ -n "$PROFILE_OUT_PREFIX_OMP" ]] && abort_if_output_exists "$PROFILE_OUT_PREFIX_OMP"
[[ -n "$PROFILE_OUT_PREFIX_GPU" ]] && abort_if_output_exists "$PROFILE_OUT_PREFIX_GPU"
[[ -n "$PROFILE_OUT_PREFIX_GPU_NSYS" ]] && abort_if_output_exists "$PROFILE_OUT_PREFIX_GPU_NSYS"
[[ -n "$PROFILE_OUT_PREFIX_GPU_NCU" ]] && abort_if_output_exists "$PROFILE_OUT_PREFIX_GPU_NCU"

{
  echo "repo_root=$REPO_ROOT"
  echo "cfg=$CFG"
  echo "profile_cfg_omp=${PROFILE_CFG_OMP:-n/a}"
  echo "profile_cfg_gpu=${PROFILE_CFG_GPU:-n/a}"
  echo "profile_cfg_gpu_nsys=${PROFILE_CFG_GPU_NSYS:-n/a}"
  echo "profile_cfg_gpu_ncu=${PROFILE_CFG_GPU_NCU:-n/a}"
  echo "out_dir=$OUT_DIR"
  echo "profile_bin_dir=$PROFILE_BIN_DIR"
  echo "profile_out_prefix_omp=${PROFILE_OUT_PREFIX_OMP:-n/a}"
  echo "profile_out_prefix_gpu=${PROFILE_OUT_PREFIX_GPU:-n/a}"
  echo "profile_out_prefix_gpu_nsys=${PROFILE_OUT_PREFIX_GPU_NSYS:-n/a}"
  echo "profile_out_prefix_gpu_ncu=${PROFILE_OUT_PREFIX_GPU_NCU:-n/a}"
  echo "omp_threads=$OMP_THREADS"
  echo "run_omp=$RUN_OMP"
  echo "run_gpu=$RUN_GPU"
  echo "skip_build=$SKIP_BUILD"
  echo "gpu_offload_info=$GPU_OFFLOAD_INFO"
  echo "omp_target_offload=MANDATORY"
  echo "cleanup_existing=$CLEANUP_EXISTING"
  echo "keep_outputs=$KEEP_OUTPUTS"
  echo "nsys_bin=${NSYS_BIN:-not_found}"
  echo "ncu_bin=${NCU_BIN:-not_found}"
  echo "vtune_bin=${VTUNE_BIN:-not_found}"
  echo "timestamp=$(date -Is)"
  git rev-parse --short HEAD 2>/dev/null | sed 's/^/git_commit=/' || true
  git rev-parse --abbrev-ref HEAD 2>/dev/null | sed 's/^/git_branch=/' || true
} >"$OUT_DIR/metadata.txt"

log "Output directory: $OUT_DIR"

if [[ "$SKIP_BUILD" -eq 0 ]]; then
  if [[ "$RUN_OMP" -eq 1 ]]; then
      log "Building OMP binary (MODE=profile) -> $PROFILE_BIN_DIR"
      make MODE=profile OMP=1 GPU=0 BINDIR="$PROFILE_BIN_DIR" >"$OUT_DIR/build_omp.log" 2>&1
  fi
  if [[ "$RUN_GPU" -eq 1 ]]; then
    log "Building GPU binary -> $PROFILE_BIN_DIR"
    make GPU=1 OMP=0 BINDIR="$PROFILE_BIN_DIR" >"$OUT_DIR/build_gpu.log" 2>&1
  fi
fi

OMP_BIN="$PROFILE_BIN_DIR/rm_synthesis_profile_omp1_gpu0"
GPU_BIN="$PROFILE_BIN_DIR/rm_synthesis_release_omp0_gpu1"

if [[ ! -x "$OMP_BIN" ]]; then
  OMP_BIN="$REPO_ROOT/bin/rm_synthesis_release_omp1_gpu0"
fi

if [[ ! -x "$GPU_BIN" ]]; then
  GPU_BIN="$REPO_ROOT/bin/rm_synthesis_release_omp0_gpu1"
fi

if [[ "$RUN_OMP" -eq 1 ]]; then
  if [[ ! -x "$OMP_BIN" ]]; then
    echo "ERROR: OMP binary missing: $OMP_BIN"
    exit 1
  fi
  OMP_CMD="OMP_NUM_THREADS=$OMP_THREADS '$OMP_BIN' '$PROFILE_CFG_OMP'"
  echo "$OMP_CMD" >"$OUT_DIR/cmd_omp.txt"
  log "Running OMP baseline"
  run_with_samplers omp "$OMP_CMD" "$OUT_DIR"
  if [[ "$KEEP_OUTPUTS" -eq 0 ]]; then
    cleanup_existing_outputs "$PROFILE_OUT_PREFIX_OMP"
  fi

  echo "SKIPPED: perf disabled in VTune profiling mode" >"$OUT_DIR/perf_omp.txt"

  # ---- VTune collections ------------------------------------------------
  # VTune is mandatory for OMP profiling. If it is unavailable, the script
  # exits earlier before reaching this block.
  cfg_prep_vtune="$(prepare_profile_cfg_and_outprefix "$CFG" "$OUT_DIR" "_vtune")"
  PROFILE_CFG_VTUNE="${cfg_prep_vtune%%|*}"
  PROFILE_OUT_PREFIX_VTUNE="${cfg_prep_vtune##*|}"
  abort_if_output_exists "$PROFILE_OUT_PREFIX_VTUNE"
  VTUNE_APP_CMD=(env OMP_NUM_THREADS="$OMP_THREADS" "$OMP_BIN" "$PROFILE_CFG_VTUNE")

  # 1. Hotspots: which functions use the most CPU time
  log "VTune: collecting hotspots"
  rm -rf "$OUT_DIR/vtune_hotspots"
  "$VTUNE_BIN" -collect hotspots \
    -result-dir "$OUT_DIR/vtune_hotspots" \
    -app-working-dir "$REPO_ROOT" \
    -knob sampling-mode=hw \
    -- "${VTUNE_APP_CMD[@]}" >"$OUT_DIR/vtune_hotspots_run.log" 2>&1
  if [[ -d "$OUT_DIR/vtune_hotspots" ]]; then
    "$VTUNE_BIN" -report hotspots \
      -result-dir "$OUT_DIR/vtune_hotspots" \
      -format text -report-output "$OUT_DIR/vtune_hotspots.txt" 2>&1
  else
    echo "ERROR: VTune hotspots did not create a result directory" >&2
    exit 1
  fi
  if [[ "$KEEP_OUTPUTS" -eq 0 ]]; then
    cleanup_existing_outputs "$PROFILE_OUT_PREFIX_VTUNE"
  fi

  # 2. Threading: OpenMP imbalance, wait time, concurrency
  log "VTune: collecting threading analysis"
  rm -rf "$OUT_DIR/vtune_threading"
  "$VTUNE_BIN" -collect threading \
    -result-dir "$OUT_DIR/vtune_threading" \
    -app-working-dir "$REPO_ROOT" \
    -- "${VTUNE_APP_CMD[@]}" >"$OUT_DIR/vtune_threading_run.log" 2>&1
  if [[ -d "$OUT_DIR/vtune_threading" ]]; then
    "$VTUNE_BIN" -report summary \
      -result-dir "$OUT_DIR/vtune_threading" \
      -format text -report-output "$OUT_DIR/vtune_threading.txt" 2>&1
  else
    echo "ERROR: VTune threading did not create a result directory" >&2
    exit 1
  fi
  if [[ "$KEEP_OUTPUTS" -eq 0 ]]; then
    cleanup_existing_outputs "$PROFILE_OUT_PREFIX_VTUNE"
  fi

  # 3. Memory access: cache misses, bandwidth, NUMA locality
  log "VTune: collecting memory-access analysis"
  rm -rf "$OUT_DIR/vtune_memory"
  "$VTUNE_BIN" -collect memory-access \
    -result-dir "$OUT_DIR/vtune_memory" \
    -app-working-dir "$REPO_ROOT" \
    -- "${VTUNE_APP_CMD[@]}" >"$OUT_DIR/vtune_memory_run.log" 2>&1
  if [[ -d "$OUT_DIR/vtune_memory" ]]; then
    "$VTUNE_BIN" -report summary \
      -result-dir "$OUT_DIR/vtune_memory" \
      -format text -report-output "$OUT_DIR/vtune_memory.txt" 2>&1
  else
    echo "ERROR: VTune memory-access did not create a result directory" >&2
    exit 1
  fi
  if [[ "$KEEP_OUTPUTS" -eq 0 ]]; then
    cleanup_existing_outputs "$PROFILE_OUT_PREFIX_VTUNE"
  fi

  # Consolidated vtune summary
  {
    echo "# VTune Summary"
    echo "Generated: $(date -Is)"
    echo ""
    echo "## Hotspots (CPU time by function)"
    head -n 60 "$OUT_DIR/vtune_hotspots.txt" 2>/dev/null || echo "no data"
    echo ""
    echo "## Threading (OpenMP concurrency/wait)"
    head -n 60 "$OUT_DIR/vtune_threading.txt" 2>/dev/null || echo "no data"
    echo ""
    echo "## Memory Access (cache/NUMA/bandwidth)"
    head -n 60 "$OUT_DIR/vtune_memory.txt" 2>/dev/null || echo "no data"
  } >"$OUT_DIR/vtune_summary.txt"

  log "VTune collections done. Raw result dirs:"
  log "  $OUT_DIR/vtune_hotspots  (open in vtune-gui for interactive call tree)"
  log "  $OUT_DIR/vtune_threading"
  log "  $OUT_DIR/vtune_memory"
fi

if [[ "$RUN_GPU" -eq 1 ]]; then
  if [[ ! -x "$GPU_BIN" ]]; then
    echo "ERROR: GPU binary missing: $GPU_BIN"
    exit 1
  fi
  GPU_ENV_PREFIX="OMP_TARGET_OFFLOAD=MANDATORY"
  if [[ "$GPU_OFFLOAD_INFO" -eq 1 ]]; then
    GPU_CMD="$GPU_ENV_PREFIX LIBOMPTARGET_INFO=16 '$GPU_BIN' '$PROFILE_CFG_GPU'"
  else
    GPU_CMD="$GPU_ENV_PREFIX '$GPU_BIN' '$PROFILE_CFG_GPU'"
  fi
  echo "$GPU_CMD" >"$OUT_DIR/cmd_gpu.txt"
  log "Running GPU baseline"
  run_with_samplers gpu "$GPU_CMD" "$OUT_DIR"
  if [[ "$KEEP_OUTPUTS" -eq 0 ]]; then
    cleanup_existing_outputs "$PROFILE_OUT_PREFIX_GPU"
  fi

  if [[ -n "$NSYS_BIN" ]]; then
    log "Collecting nsys timeline"
    if [[ "$GPU_OFFLOAD_INFO" -eq 1 ]]; then
      GPU_CMD_NSYS="$GPU_ENV_PREFIX LIBOMPTARGET_INFO=16 '$GPU_BIN' '$PROFILE_CFG_GPU_NSYS'"
    else
      GPU_CMD_NSYS="$GPU_ENV_PREFIX '$GPU_BIN' '$PROFILE_CFG_GPU_NSYS'"
    fi
    "$NSYS_BIN" profile --stats=true -o "$OUT_DIR/nsys_gpu" bash -lc "$GPU_CMD_NSYS" >"$OUT_DIR/nsys_run.log" 2>&1 || true
    "$NSYS_BIN" stats --force-export=true "$OUT_DIR/nsys_gpu.nsys-rep" >"$OUT_DIR/nsys_gpu_stats.txt" 2>&1 || true
    if [[ "$KEEP_OUTPUTS" -eq 0 ]]; then
      cleanup_existing_outputs "$PROFILE_OUT_PREFIX_GPU_NSYS"
    fi
  else
    echo "SKIPPED: nsys not found" >"$OUT_DIR/nsys_gpu_stats.txt"
  fi

  if [[ -n "$NCU_BIN" ]]; then
    log "Collecting ncu kernel metrics"
    if [[ "$GPU_OFFLOAD_INFO" -eq 1 ]]; then
      GPU_CMD_NCU="$GPU_ENV_PREFIX LIBOMPTARGET_INFO=16 '$GPU_BIN' '$PROFILE_CFG_GPU_NCU'"
    else
      GPU_CMD_NCU="$GPU_ENV_PREFIX '$GPU_BIN' '$PROFILE_CFG_GPU_NCU'"
    fi
    "$NCU_BIN" --set speedOfLight --target-processes all --csv --log-file "$OUT_DIR/ncu_gpu.csv" \
      bash -lc "$GPU_CMD_NCU" >"$OUT_DIR/ncu_run.log" 2>&1 || true
    if [[ "$KEEP_OUTPUTS" -eq 0 ]]; then
      cleanup_existing_outputs "$PROFILE_OUT_PREFIX_GPU_NCU"
    fi
  else
    echo "SKIPPED: ncu not found" >"$OUT_DIR/ncu_gpu.csv"
  fi
fi

{
  echo "RM-Synthesis Bottleneck Classification"
  echo "Generated: $(date -Is)"
  echo ""
} >"$OUT_DIR/classify_report.txt"

if [[ "$RUN_OMP" -eq 1 ]]; then
  classify_omp "$OUT_DIR"
fi
if [[ "$RUN_GPU" -eq 1 ]]; then
  classify_gpu "$OUT_DIR"
fi

generate_consolidated_report "$OUT_DIR"

log "Done. See:"
log "  $OUT_DIR/classify_report.txt"
log "  $OUT_DIR/consolidated_report.md"
if [[ -f "$OUT_DIR/consolidated_report.html" ]]; then
  log "  $OUT_DIR/consolidated_report.html"
fi
if [[ -f "$OUT_DIR/vtune_summary.txt" ]]; then
  log "  $OUT_DIR/vtune_summary.txt"
fi
if [[ -d "$OUT_DIR/vtune_hotspots" ]]; then
  log "  vtune-gui $OUT_DIR/vtune_hotspots   # open hotspot result in GUI"
  log "  vtune-gui $OUT_DIR/vtune_threading  # open threading result in GUI"
  log "  vtune-gui $OUT_DIR/vtune_memory     # open memory result in GUI"
fi
log "  $OUT_DIR (all raw logs)"
[[ -n "$PROFILE_OUT_PREFIX_OMP" ]] && log "  OMP science outputs at prefix: $PROFILE_OUT_PREFIX_OMP"
[[ -n "$PROFILE_OUT_PREFIX_GPU" ]] && log "  GPU science outputs at prefix: $PROFILE_OUT_PREFIX_GPU"
[[ -n "$PROFILE_OUT_PREFIX_GPU_NSYS" ]] && log "  GPU nsys outputs at prefix: $PROFILE_OUT_PREFIX_GPU_NSYS"
[[ -n "$PROFILE_OUT_PREFIX_GPU_NCU" ]] && log "  GPU ncu outputs at prefix: $PROFILE_OUT_PREFIX_GPU_NCU"
