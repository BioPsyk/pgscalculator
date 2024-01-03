#!/usr/bin/env bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${script_dir}/init-containerization.sh"

mkdir -p tmp_singularity_build
SINGULARITY_TMPDIR="tmp_singularity_build"

exec singularity build tmp/${singularity_image_tag} \
     docker-daemon:"${image_tag}"

