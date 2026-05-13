#!/usr/bin/env bash
set -euo pipefail

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "${script_dir}/init-containerization.sh"

cd "${project_dir}"

target="biopsyk/${deploy_image_tag}"

echo ">> Pushing multi-arch deployment docker image to Docker Hub"

# Tag the architecture-specific images for Docker Hub
echo "Tagging images for Docker Hub..."
docker tag "${deploy_image_tag}-amd64" "${target}-amd64"
docker tag "${deploy_image_tag}-arm64" "${target}-arm64"

platform_digest() {
  local image_ref="$1"
  local platform="$2"

  docker buildx imagetools inspect "${image_ref}" | awk -v platform="${platform}" '
    $1 == "Name:" && $2 ~ /@sha256:/ { name = $2 }
    $1 == "Platform:" && $2 == platform {
      sub(/^.*@/, "", name)
      print name
      exit
    }
  '
}

# Push the architecture-specific images and capture their platform manifest digests.
# Docker may push these tags as OCI indexes when attestations are enabled, so the
# top-level digest is not always the arch-specific manifest digest.
echo "Pushing architecture-specific images..."
docker push "${target}-amd64"
docker push "${target}-arm64"

amd64_digest=$(platform_digest "docker.io/${target}-amd64" "linux/amd64")
arm64_digest=$(platform_digest "docker.io/${target}-arm64" "linux/arm64")

[[ -z "${amd64_digest}" ]] && { echo "Error: could not resolve amd64 manifest digest"; exit 1; }
[[ -z "${arm64_digest}" ]] && { echo "Error: could not resolve arm64 manifest digest"; exit 1; }

echo "AMD64 digest: ${amd64_digest}"
echo "ARM64 digest: ${arm64_digest}"

# Create and push the multi-arch manifest
echo "Creating and pushing multi-arch manifest..."

# Create new manifest with explicit platform manifest digests
docker buildx imagetools create -t "${target}" \
  "docker.io/${target}-amd64@${amd64_digest}" \
  "docker.io/${target}-arm64@${arm64_digest}"

# Wait a moment for Docker Hub to process the manifest
sleep 5

echo "Multi-arch image successfully pushed to Docker Hub as ${target}"
echo "Verifying manifest..."
docker buildx imagetools inspect "${target}"
