#!/usr/bin/env bash
# Run the selected q/rh/radiation combinations sequentially using the active Python environment.
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
zonal_python="${PYTHON:-python}"

# Edit these defaults, or override them when launching this script.
read -r -a zonal_qflux_values <<< "${ZONAL_QFLUX_VALUES:-0 25 50}"
read -r -a zonal_rh_values <<< "${ZONAL_RH_VALUES:-60 70 80}"
read -r -a zonal_rad_values <<< "${ZONAL_RAD_VALUES:-original _lock}"
zonal_spinup="${ZONAL_SPINUP:-${HEATING_SPINUP:-1800}}"
zonal_num_threads="${ZONAL_NUM_THREADS:-${HEATING_NUM_THREADS:-1}}"
zonal_data_root="${ZONAL_DATA_ROOT:-${HEATING_DATA_ROOT:-/data92/b11209013/Isca_Gary/DATA}}"
zonal_low_period="${ZONAL_LOW_PERIOD_DAYS:-30}"
zonal_high_period="${ZONAL_HIGH_PERIOD_DAYS:-10}"

# Optional ZONAL_OUTPUT_DIR and ZONAL_EARTH_RADIUS_M
# environment settings are inherited directly by the Python script.
for zonal_qflux in "${zonal_qflux_values[@]}"; do
    for zonal_rh in "${zonal_rh_values[@]}"; do
        for zonal_rad in "${zonal_rad_values[@]}"; do
            printf 'Running zonal DSE transport: q%s_rh%s_b%s\n' \
                "$zonal_qflux" "$zonal_rh" "$zonal_rad"
            HEATING_QFLUX="$zonal_qflux" HEATING_RH="$zonal_rh" \
                HEATING_RAD="$zonal_rad" HEATING_SPINUP="$zonal_spinup" \
                HEATING_NUM_THREADS="$zonal_num_threads" HEATING_DATA_ROOT="$zonal_data_root" \
                ZONAL_LOW_PERIOD_DAYS="$zonal_low_period" \
                ZONAL_HIGH_PERIOD_DAYS="$zonal_high_period" \
                "$zonal_python" -u "$script_dir/zonal_energy_transpose.py"
        done
    done
done

printf 'All %s zonal DSE transport analyses completed.\n' \
    "$(( ${#zonal_qflux_values[@]} * ${#zonal_rh_values[@]} * ${#zonal_rad_values[@]} ))"
