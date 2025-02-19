#!/usr/bin/env bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "${script_dir}/init-containerization.sh"

cd "${project_dir}"

echo ">> Setting up docker buildx for multi-arch support"
docker buildx create --name multiarch --driver docker-container --use || true
docker buildx inspect --bootstrap

echo ">> Building multi-arch docker image"

docker buildx build \
  --platform linux/amd64,linux/arm64 \
  --progress=plain \
  --tag "${image_tag}" \
  --load \
  ./docker \
  "$@"
