#!/usr/bin/env bash
set -euo pipefail

# Build and push rmtool image from a specific git tag snapshot.
# Workflow:
#   1) create a temporary clone under tmp/
#   2) git pull (with tags)
#   3) checkout the requested git tag
#   4) build using docker/dockerfile from that checked-out source
#   5) push only if target Docker Hub tag does not already exist
# Usage:
#   docker/build_push_from_tag.sh <git_tag> [docker_image_tag]
# Example:
#   docker/build_push_from_tag.sh rmcube-stat rmcube-stat
#   docker/build_push_from_tag.sh v1.2.0

REPO="wasimraja81/rmtool"
GIT_REMOTE="${GIT_REMOTE:-origin}"
DEFAULT_REMOTE_URL="https://github.com/wasimraja81/rmtool.git"

usage() {
  cat <<'EOF'
Usage:
  docker/build_push_from_tag.sh <git_tag> [docker_image_tag]

Builds and pushes wasimraja81/rmtool:<docker_image_tag> from a checked-out git tag.

Arguments:
  <git_tag>           Required source git tag to checkout in temporary clone
  [docker_image_tag]  Optional Docker image tag (defaults to <git_tag>)

Environment:
  GIT_REMOTE       Git remote name when discovering URL from local repo (default: origin)
  GIT_REMOTE_URL   Explicit git clone URL override
  TMP_BASE_DIR     Base directory for temporary clone work
EOF
}

if [[ $# -eq 1 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
  usage
  exit 0
fi

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage >&2
  exit 1
fi

GIT_TAG="$1"
IMAGE_TAG="${2:-$GIT_TAG}"
IMAGE_REF="${REPO}:${IMAGE_TAG}"

if [[ "$GIT_TAG" == -* ]]; then
  echo "Error: invalid git tag '${GIT_TAG}'." >&2
  echo "Hint: use -h or --help for usage." >&2
  exit 1
fi

for cmd in git docker mktemp; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: required command '$cmd' is not available." >&2
    exit 1
  fi
done

ROOT_DIR="$(git rev-parse --show-toplevel 2>/dev/null || true)"
IN_REPO=0
if [[ -n "$ROOT_DIR" ]]; then
  IN_REPO=1
fi

echo "Checking Docker Hub tag availability for ${IMAGE_REF} ..."
if docker manifest inspect "$IMAGE_REF" >/dev/null 2>&1; then
  echo ""
  echo "Refusing to overwrite existing remote image tag: ${IMAGE_REF}" >&2
  echo "Delete the existing tag on Docker Hub first, then rebuild and push." >&2
  echo ""
  exit 2
fi

if [[ -n "${GIT_REMOTE_URL:-}" ]]; then
  remote_url="$GIT_REMOTE_URL"
elif [[ "$IN_REPO" -eq 1 ]]; then
  remote_url="$(git -C "$ROOT_DIR" remote get-url "$GIT_REMOTE" 2>/dev/null || true)"
else
  remote_url="$DEFAULT_REMOTE_URL"
fi

if [[ -z "$remote_url" ]]; then
  echo "Error: unable to resolve remote URL for '$GIT_REMOTE'." >&2
  echo "Set GIT_REMOTE_URL explicitly, e.g. GIT_REMOTE_URL=https://github.com/<user>/<repo>.git" >&2
  exit 1
fi

if [[ -n "${TMP_BASE_DIR:-}" ]]; then
  tmp_base="$TMP_BASE_DIR"
elif [[ "$IN_REPO" -eq 1 ]]; then
  tmp_base="$ROOT_DIR/tmp"
else
  tmp_base="$(pwd)/tmp"
fi

mkdir -p "$tmp_base"
safe_tag="${GIT_TAG//[^A-Za-z0-9._-]/_}"
tmpdir="$(mktemp -d "$tmp_base/rmtool-build-${safe_tag}-XXXXXX")"
clone_dir="$tmpdir/repo"

cleanup() {
  rm -rf "$tmpdir"
}
trap cleanup EXIT

echo "Cloning '${remote_url}' into temporary directory: $clone_dir"
if ! git clone "$remote_url" "$clone_dir" >/dev/null; then
  echo "Error: unable to clone remote repository '${remote_url}'." >&2
  exit 1
fi

pushd "$clone_dir" >/dev/null

echo "Running git pull (with tags) in temporary clone ..."
if ! git pull --tags --ff-only "$GIT_REMOTE" >/dev/null; then
  echo "Error: git pull failed in temporary clone." >&2
  popd >/dev/null
  exit 1
fi

echo "Checking out tag '${GIT_TAG}' in temporary clone ..."
if ! git checkout "tags/${GIT_TAG}" >/dev/null; then
  echo "Error: unable to checkout tag '${GIT_TAG}' in temporary clone." >&2
  popd >/dev/null
  exit 1
fi

if [[ ! -f docker/dockerfile ]]; then
  echo "Error: docker/dockerfile not found in tag '${GIT_TAG}'." >&2
  popd >/dev/null
  exit 1
fi

echo "Building image ${IMAGE_REF} from git tag ${GIT_TAG} ..."
docker build \
  --pull \
  --build-arg RMTOOL_GIT_TAG="$GIT_TAG" \
  -f docker/dockerfile \
  -t "$IMAGE_REF" \
  .

popd >/dev/null

echo "Pushing ${IMAGE_REF} ..."
docker push "$IMAGE_REF"

echo "Done: ${IMAGE_REF}"
