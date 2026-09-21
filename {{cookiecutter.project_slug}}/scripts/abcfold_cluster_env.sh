#!/usr/bin/env bash
# scripts/abcfold_cluster_env.sh — SOURCE this from the run_abcfold rule.
# ============================================================================
# Mined verbatim from workflows/processing/_reference/submit_abcfold.sh
# (the DRB2 pipeline's proven SLURM submitter). A Snakemake `shell:` block on
# a compute node is its own non-login process: `module` is undefined and no
# conda env is active, exactly the conditions that cost real debugging time
# on IFB (2026-08-21). This script fixes all of that in one `source`.
#
# Usage:   source scripts/abcfold_cluster_env.sh
# Reads (env, with the defaults below):
#   AF3_MODULE_VERSION      default 3.0.2
#   CUDA_MODULE_VERSION     default 12.9.1
#   ABCFOLD_ENV             default /shared/projects/npf_abinitio/conda/envs/abcfold-npf-pipeline
#   AF3_SIF_PATH            optional explicit override (empty = auto-discover)
# Exports:
#   ABCFOLD_BIN   AF3_SIF_PATH   CUDA_HOME
# and runs `module load singularity` + `module load cuda-toolkit/$VER`.
# ============================================================================

AF3_MODULE_VERSION="${AF3_MODULE_VERSION:-3.0.2}"
CUDA_MODULE_VERSION="${CUDA_MODULE_VERSION:-12.9.1}"
ABCFOLD_ENV="${ABCFOLD_ENV:-/shared/projects/npf_abinitio/conda/envs/abcfold-npf-pipeline}"
AF3_SIF_PATH="${AF3_SIF_PATH:-}"

# neither a non-login remote shell nor a SLURM batch script sources this
[[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh

discover_af3_sif() {
    local wrapper_dir="/shared/software/singularity/wrappers/alphafold/$AF3_MODULE_VERSION"
    local wrapper="$wrapper_dir/run_alphafold.py"
    [[ -f "$wrapper" ]] || return 1
    local sif
    sif=$(grep -oE '[^[:space:]]+\.sif' "$wrapper" | head -n1)
    [[ -z "$sif" ]] && return 1
    [[ "$sif" != /* ]] && sif="$wrapper_dir/$sif"
    [[ -f "$sif" ]] || return 1
    echo "$sif"
}

discover_cuda_home() {
    local modulefile="/shared/software/modulefiles/cuda-toolkit/$CUDA_MODULE_VERSION"
    [[ -f "$modulefile" ]] || return 1
    local env_bin env_root targets_root
    env_bin=$(grep -oE '[^[:space:]]+/envs/cuda-toolkit-[^[:space:]]+/bin' "$modulefile" | head -n1)
    [[ -z "$env_bin" ]] && return 1
    env_root=$(dirname "$env_bin")
    targets_root=$(find "$env_root/targets" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | head -n1)
    if [[ -z "$targets_root" || ! -f "$targets_root/include/cuda_runtime_api.h" ]]; then
        echo "$env_root"; return 0
    fi
    local shim="/shared/projects/npf_abinitio/conda/cuda_home_shim"
    mkdir -p "$shim"
    ln -sfn "$env_root/bin" "$shim/bin"
    ln -sfn "$env_root/nvvm" "$shim/nvvm"
    ln -sfn "$targets_root/include" "$shim/include"
    ln -sfn "$targets_root/lib" "$shim/lib"
    ln -sfn "$targets_root/lib" "$shim/lib64"
    echo "$shim"
}

if [[ -z "$AF3_SIF_PATH" ]]; then
    if _sif=$(discover_af3_sif); then
        AF3_SIF_PATH="$_sif"
        echo "[abcfold_cluster_env] AF3 .sif: $AF3_SIF_PATH"
    else
        echo "[abcfold_cluster_env] WARNING: could not auto-discover the AF3 .sif for" \
             "alphafold/$AF3_MODULE_VERSION — AF3 (-a) will fall back to 'docker run'" >&2
    fi
fi

if _cuda=$(discover_cuda_home); then
    export CUDA_HOME="$_cuda"
    echo "[abcfold_cluster_env] CUDA_HOME: $CUDA_HOME"
else
    echo "[abcfold_cluster_env] WARNING: could not auto-discover CUDA_HOME from" \
         "cuda-toolkit/$CUDA_MODULE_VERSION — Protenix (-p) may fail to build" >&2
fi

command -v module >/dev/null 2>&1 && {
    module load singularity 2>/dev/null || true
    module load "cuda-toolkit/$CUDA_MODULE_VERSION" 2>/dev/null || true
}

ABCFOLD_BIN="$ABCFOLD_ENV/bin/abcfold"
[[ -x "$ABCFOLD_BIN" ]] || ABCFOLD_BIN="abcfold"
export ABCFOLD_BIN AF3_SIF_PATH
echo "[abcfold_cluster_env] abcfold: $ABCFOLD_BIN"
