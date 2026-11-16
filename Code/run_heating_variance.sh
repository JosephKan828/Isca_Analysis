#!/usr/bin/env bash
# Run the selected q/rh/radiation combinations sequentially using the active Python environment.
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
variance_python="${PYTHON:-python}"

# Inputs: edit the defaults here, or override them when launching this script.
read -r -a variance_qflux_values <<< "${HEATING_QFLUX_VALUES:-0 25 50}"
read -r -a variance_rh_values <<< "${HEATING_RH_VALUES:-60 70 80}"
read -r -a variance_rad_values <<< "${HEATING_RAD_VALUES:-${HEATING_RAD:-original _lock}}"
variance_spinup="${HEATING_SPINUP:-1800}"

for variance_qflux in "${variance_qflux_values[@]}"; do
    for variance_rh in "${variance_rh_values[@]}"; do
        for variance_rad in "${variance_rad_values[@]}"; do
            printf 'Running heating variance: q%s_rh%s_b%s\n' \
                "$variance_qflux" "$variance_rh" "$variance_rad"
            HEATING_QFLUX="$variance_qflux" HEATING_RH="$variance_rh" \
                HEATING_RAD="$variance_rad" HEATING_SPINUP="$variance_spinup" \
                "$variance_python" -u "$script_dir/Heating_variance.py"
        done
    done
done

printf 'All %s heating-variance analyses completed.\n' \
    "$(( ${#variance_qflux_values[@]} * ${#variance_rh_values[@]} * ${#variance_rad_values[@]} ))"
