#!/bin/bash

# Concatenate model output from different Isca experiments.

DATA="/work/garywu/isca_frierson_byrne_18_15y/data"
OUTPUT="/data92/b11209013/Isca_Gary/DATA"

shopt -s nullglob

# first document the chunks for experiments
# the output file include: ("atmos_daily.nc", "pbl_daily.nc", "surface_daily.nc", "surface_snapshot_daily.nc")
# document the file name through loop first
for exp_dir in "$DATA"/*lock; do
    [[ -d "$exp_dir" ]] || continue
    exp="${exp_dir##*/}"

    # create output directory for the experiment
    mkdir -p "$OUTPUT/$exp"

    for file_name in "atmos_daily.nc" "pbl_daily.nc" "surface_daily.nc" "surface_snapshot_daily.nc"; do

        echo "Processing $file_name for experiment $exp"

        runs=("$exp_dir"/run????/"$file_name")
        if (( ${#runs[@]} == 0 )); then
            printf 'No %s files found for experiment %s; skipping.\n' "$file_name" "$exp" >&2
            continue
        fi
        
        # concatenate the runs for the current file_name
        cdo -P 8 mergetime "${runs[@]}" "$OUTPUT/$exp/$file_name"
        
    done

done
