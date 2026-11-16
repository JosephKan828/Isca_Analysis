#!/usr/bin/env bash
# Run the 18 q/rh/radiation combinations sequentially using the active Python environment.
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
budget_python="${PYTHON:-python}"

for budget_qflux in 0 25 50; do
    for budget_rh in 60 70 80; do
        for budget_rad in original _lock; do
            printf 'Running EAPE budget: q%s_rh%s_b%s\n' \
                "$budget_qflux" "$budget_rh" "$budget_rad"
            EAPE_QFLUX="$budget_qflux" EAPE_RH="$budget_rh" \
                EAPE_RAD="$budget_rad" \
                "$budget_python" -u "$script_dir/EAPE_budget.py"
        done
    done
done

printf 'All 18 EAPE budget analyses completed.\n'
