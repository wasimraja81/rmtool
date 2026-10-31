#!/usr/bin/env bash
# Every command .github/workflows/build-and-verify.yml runs against the
# Docker images, in one place: CI calls this script, and so does a local
# run before tagging a cut, so the two cannot drift apart.
#
#   docker/ci_run.sh suite <image-tag>
#       Build docker/dockerfile, check the tools are on PATH, run
#       tests/run_tests.sh inside the image.
#   docker/ci_run.sh fresh-install <image-tag>
#       Build docker/dockerfile.docs-verify, run
#       docker/verify_fresh_install.sh (the documented install path).
#
# Kept identical between CI and a local run, by construction:
# - Build context: only git-tracked files (with their working-tree
#   contents), streamed in as a tar -- what CI's fresh checkout holds.
#   Never local leftovers: test outputs, fixtures, FFTW wisdom, editor
#   settings.
# - CPUs: pinned to the first 2, so a local machine with more CPUs, or
#   a larger runner, does not change the conditions the tests run under
#   (CI's runner log shows 1 physical core, i.e. 2 CPUs with SMT).
# - Suite user: the caller's own UID (--user), so files written to the
#   bind-mounted test folders belong to whoever ran it. The image is
#   built with its owner at UID/GID 4242, so the caller never owns
#   /app/rmtool -- the condition every `docker run --user` user meets --
#   whatever UID runs this.
# - Suite start state: the bind-mounted test folders must be empty, as
#   in CI's fresh checkout. They default to tests/output and tests/data;
#   set CI_RUN_OUTPUT_DIR/CI_RUN_DATA_DIR to use other (empty) folders.
set -euo pipefail

IMAGE_OWNER_ID=4242
CPUSET=0-1

usage() {
    echo "Usage: $0 suite|fresh-install <image-tag>" >&2
    exit 2
}

[[ $# -eq 2 ]] || usage
mode="$1"
image_tag="$2"

cd "$(git rev-parse --show-toplevel)"

build_from_tracked_files() {
    local dockerfile="$1" tag="$2"
    shift 2
    echo "################ docker build -f $dockerfile $* (git-tracked files only) ################"
    git ls-files -z | tar --null -T - -cf - | docker build -f "$dockerfile" "$@" -t "$tag" -
}

require_empty_dir() {
    local dir="$1"
    mkdir -p "$dir"
    if [[ -n "$(ls -A "$dir")" ]]; then
        echo "ERROR: $dir is not empty -- CI starts from an empty folder; empty it or set CI_RUN_OUTPUT_DIR/CI_RUN_DATA_DIR to an empty one." >&2
        exit 1
    fi
}

case "$mode" in
    suite)
        caller_uid="$(id -u)"
        caller_gid="$(id -g)"
        if [[ "$caller_uid" -eq 0 || "$caller_uid" -eq "$IMAGE_OWNER_ID" ]]; then
            echo "ERROR: run as a non-root user other than UID $IMAGE_OWNER_ID (the image's owner): the suite must run as a user who does not own /app/rmtool." >&2
            exit 1
        fi
        output_dir="${CI_RUN_OUTPUT_DIR:-tests/output}"
        data_dir="${CI_RUN_DATA_DIR:-tests/data}"
        require_empty_dir "$output_dir"
        require_empty_dir "$data_dir"
        output_dir="$(cd "$output_dir" && pwd)"
        data_dir="$(cd "$data_dir" && pwd)"

        build_from_tracked_files docker/dockerfile "$image_tag" \
            --build-arg RMTOOL_UID="$IMAGE_OWNER_ID" \
            --build-arg RMTOOL_GID="$IMAGE_OWNER_ID"

        echo "################ casda_fetch.py runs by bare name ################"
        docker run --rm "$image_tag" casda_fetch.py -h
        echo "################ astroquery/pyvo/requests import ################"
        docker run --rm "$image_tag" python3 -c "import astroquery, pyvo, requests"
        echo "################ run_pipeline.sh and every compiled binary are on PATH ################"
        docker run --rm "$image_tag" bash -c "
            which run_pipeline.sh casda_fetch.py &&
            which rm_synthesis match_cubes rmclean_cubes convolve_cubes reproject_cubes
        "

        # tests/output and tests/data are bind-mounted (not the image's own
        # copies) so dry_run=y's disk-type check resolves against a host
        # block device via df, as it would for a user's own bind-mounted
        # output folder (see src/convolve_cubes.f90's run_dry_run_check).
        echo "################ tests/run_tests.sh as UID $caller_uid, CPUs $CPUSET ################"
        docker run --rm \
            --cpuset-cpus="$CPUSET" \
            --user "$caller_uid:$caller_gid" \
            -v "$output_dir:/app/rmtool/tests/output" \
            -v "$data_dir:/app/rmtool/tests/data" \
            "$image_tag" \
            bash tests/run_tests.sh
        ;;
    fresh-install)
        build_from_tracked_files docker/dockerfile.docs-verify "$image_tag"
        # .dockerignore keeps docker/ out of every image, so the script is
        # bind-mounted in rather than COPY'd.
        echo "################ docker/verify_fresh_install.sh, CPUs $CPUSET ################"
        docker run --rm \
            --cpuset-cpus="$CPUSET" \
            -v "$PWD/docker/verify_fresh_install.sh:/verify_fresh_install.sh:ro" \
            "$image_tag" \
            bash /verify_fresh_install.sh
        ;;
    *)
        usage
        ;;
esac
