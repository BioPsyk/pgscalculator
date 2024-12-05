#!/usr/bin/env bash

# Docker run configuration
docker_run_args="--memory=20g \
                 --cpus=6 \
                 --ulimit nofile=1024:1024 \
                 --security-opt=no-new-privileges" 