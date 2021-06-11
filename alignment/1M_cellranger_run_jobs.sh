#!/usr/bin/env bash

set -e

mkdir -p ./logs

for script in ./jobs/run_*.sh; do
    lane_id=$(basename $script)
    sh $script 2>&1 | tee ./logs/${lane_id}.txt
done
