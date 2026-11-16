"""
    Calculate vertical variance profiles of convective and condensational heating.
    Pool time, latitude, and longitude samples separately at each vertical level.
    Save both variance profiles and pressure coordinates to a compressed NPZ file.
    Isca heating is stored in K/s; variances below are calculated in K^2/day^2.
"""

# %% Configure CPU threading before importing numerical libraries
import os
import json

# Override at launch, e.g. HEATING_NUM_THREADS=2 python Code/Heating_variance.py
NUM_THREADS: int = int(os.environ.get("HEATING_NUM_THREADS", "1"))
if NUM_THREADS < 1:
    raise ValueError("HEATING_NUM_THREADS must be a positive integer")

for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[variable] = str(NUM_THREADS)

# %% Import necessary libraries
import numpy as np
import netCDF4 as nc

from pathlib import Path

# %% Load data
# Select parameters through run_heating_variance.sh or environment variables.
qflux: str = os.environ.get("HEATING_QFLUX", "0")  # intensity of ocean heat uptake
rh   : str = os.environ.get("HEATING_RH", "60")  # relative humidity for simple Betts-Miller scheme
rad  : str = os.environ.get("HEATING_RAD", "original")  # radiation intensity; ["0", "original"]

# cutoff first 1800 days (spin-up period)
spinup: int = int(os.environ.get("HEATING_SPINUP", "1800"))
if spinup < 0:
    raise ValueError("HEATING_SPINUP must be nonnegative")

# data path
DATA_PATH: Path = Path(os.environ.get(
    "HEATING_DATA_ROOT", "/data92/b11209013/Isca_Gary/DATA",
)) / f"q{qflux}_rh{rh}_b{rad}"

# %% Calculate variance across time, latitude, and longitude at each level
variance: dict[str, np.ndarray] = {}
with nc.Dataset(DATA_PATH / "atmos_daily.nc") as ds:
    time_units = getattr(ds.variables["time"], "units", "days")
    time_calendar = getattr(ds.variables["time"], "calendar", "standard")

    coords: dict[str, np.ndarray] = {
        key: ds.variables[key][:]
        for key in ("time", "pfull", "lat", "lon")
    }

    # Match the tropical troposphere selected in EAPE_budget.py (aquaplanet).
    lat_mask = np.abs(coords["lat"]) <= 5.0
    level_mask = (coords["pfull"] >= 100.0) & (coords["pfull"] <= 1000.0)
    if not np.any(lat_mask) or not np.any(level_mask):
        raise ValueError("No grid points fall within 5 degrees latitude and 100–1000 hPa.")
    if coords["time"].size - spinup < 2:
        raise ValueError("At least two time samples are required after spinup.")

    # Read one heating component at a time to limit memory use.
    for key in ("dt_tg_convection", "dt_tg_condensation"):
        print(f"Calculating variance: {key}", flush=True)
        data: np.ndarray = np.ma.asarray(
            ds.variables[key][spinup:, level_mask, lat_mask, :], dtype=np.float64,
        ).filled(np.nan)  # shape: (time, level, lat, lon)
        data *= 86400.0  # K/s to K/day; variance therefore has units K^2/day^2

        # Each level has its own mean over all time/latitude/longitude samples.
        # Population variance (ddof=0); missing values are excluded.
        variance[key] = np.nanvar(data, axis=(0, 2, 3), ddof=0)
        del data

coords["time"] = coords["time"][spinup:]
coords["lat"] = coords["lat"][lat_mask]
coords["pfull"] = coords["pfull"][level_mask]
pressure_order = np.argsort(coords["pfull"])

# %% Save numerical diagnostics for separate visualization
script_dir = (
    Path(__file__).resolve().parent if "__file__" in globals()
    else Path("/home/b11209013/Isca_analysis/Code")  # IDE cell execution
)
output_dir = Path(os.environ.get(
    "HEATING_OUTPUT_DIR", str(script_dir.parent / "Output" / "Heating_variance"),
))
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / f"q{qflux}_rh{rh}_b{rad}_variance.npz"

metadata = {
    "title": "Vertical convective and condensational heating variance",
    "analysis_script": "Heating_variance.py",
    "source_file": str(DATA_PATH / "atmos_daily.nc"),
    "qflux": qflux, "rh": rh, "rad": rad,
    "spinup_records": spinup,
    "retained_time_samples": len(coords["time"]),
    "latitude_samples": len(coords["lat"]),
    "longitude_samples": len(coords["lon"]),
    "latitude_min": float(np.min(coords["lat"])),
    "latitude_max": float(np.max(coords["lat"])),
    "pressure_min_hpa": float(np.min(coords["pfull"])),
    "pressure_max_hpa": float(np.max(coords["pfull"])),
    "time_start": float(coords["time"][0]),
    "time_end": float(coords["time"][-1]),
    "time_units": time_units, "time_calendar": time_calendar,
    "variance_ddof": 0,
    "variance_sample_axes": "time latitude longitude",
    "sample_weighting": "Equal weight for each valid sample; missing values excluded",
    "heating_conversion": "Input K/s multiplied by 86400 before calculating variance",
}

# Named numerical arrays and JSON metadata can be loaded with allow_pickle=False.
output_arrays: dict[str, np.ndarray] = {
    "pfull": np.asarray(coords["pfull"])[pressure_order],
    "convection_variance": variance["dt_tg_convection"][pressure_order],
    "condensation_variance": variance["dt_tg_condensation"][pressure_order],
}
metadata["variables"] = {
    "pfull": {
        "dimensions": ["pfull"], "units": "hPa", "standard_name": "air_pressure",
        "positive": "down", "long_name": "Retained native full-level pressure in increasing order",
    },
    "convection_variance": {
        "dimensions": ["pfull"], "units": "K2 day-2",
        "long_name": "Population variance of convective heating", "source_variable": "dt_tg_convection",
    },
    "condensation_variance": {
        "dimensions": ["pfull"], "units": "K2 day-2",
        "long_name": "Population variance of condensational heating", "source_variable": "dt_tg_condensation",
    },
}

np.savez_compressed(output_path, **output_arrays, metadata=np.array(json.dumps(metadata)))
print(f"Saved {output_path}", flush=True)
# %%
