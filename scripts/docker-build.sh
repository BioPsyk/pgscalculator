#!/usr/bin/env bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "${script_dir}/init-containerization.sh"

cd "${project_dir}"

echo ">> Building docker image for local architecture"

docker buildx build \
  --platform linux/$(uname -m | sed 's/x86_64/amd64/;s/aarch64/arm64/') \
  --progress=plain \
  --load \
  ./docker \
  -t "${image_tag}" \
  "$@"
