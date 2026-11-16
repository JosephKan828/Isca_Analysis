"""
Calculate and save full-record zonal dry static energy transport diagnostics.

Filter u and s=cp*T+g*z separately along the complete retained time axis.
Low: periods >=30 days; middle: 10<period<30 days; high: period<=10 days.
The zero-frequency term is the full-record mean. All 16 products of
(mean, low, middle, high) are retained before meridional integration.
Latitude strips are summed using dy=a*dphi over 20S-20N, without
normalization. A subsequent column-mass integral gives transport in watts.

Save coordinates, frequency masks, component maps, zonal anomalies, closure
diagnostics, column transports, signed bar contributions, and JSON metadata
to a compressed NPZ. Plot the archive with Code/plot/Zonal_DSE_transport.ipynb.

All selected latitudes and retained times are loaded together. Height is
referenced to zero at the lowest full model level.
"""

# %% Configure CPU threading before importing numerical libraries
import os
import json
from pathlib import Path

NUM_THREADS = int(os.environ.get("HEATING_NUM_THREADS", "1"))
if NUM_THREADS < 1:
    raise ValueError("HEATING_NUM_THREADS must be a positive integer")
for variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[variable] = str(NUM_THREADS)

import numpy as np
import netCDF4 as nc

# %% Analysis and output settings
# Retain the existing HEATING_* experiment settings for compatibility.
qflux = os.environ.get("HEATING_QFLUX", "0")
rh = os.environ.get("HEATING_RH", "60")
rad = os.environ.get("HEATING_RAD", "original")
spinup = int(os.environ.get("HEATING_SPINUP", "1800"))
DATA_PATH = Path(os.environ.get(
    "HEATING_DATA_ROOT", "/data92/b11209013/Isca_Gary/DATA",
)) / f"q{qflux}_rh{rh}_b{rad}"

LOW_PERIOD_DAYS = float(os.environ.get("ZONAL_LOW_PERIOD_DAYS", "30"))
HIGH_PERIOD_DAYS = float(os.environ.get("ZONAL_HIGH_PERIOD_DAYS", "10"))
LATITUDE_LIMIT = 20.0
EARTH_RADIUS = float(os.environ.get("ZONAL_EARTH_RADIUS_M", "6371000"))

SCRIPT_DIR = (
    Path(__file__).resolve().parent if "__file__" in globals()
    else Path("/home/b11209013/Isca_analysis/Code")
)
OUTPUT_DIR = Path(os.environ.get(
    "ZONAL_OUTPUT_DIR", str(SCRIPT_DIR.parent / "Output" / "Zonal_energy_transpose"),
))
GRAVITY = 9.80665
RD = 287.05
CP = 1004.5
COMPONENT_NAMES = ("mean", "low", "mid", "high")
CONTRIBUTION_NAMES = ("total", "mean_flow", "low", "mid", "high", "cross", "transient")


def meridional_strip_widths(latitude, latitude_limit, radius=EARTH_RADIUS):
    """Meridian lengths of global latitude cells clipped to the requested belt.

    Cell edges are midpoints of neighboring centers, with polar outer edges.
    Partially intersecting edge cells are included even if their centers are
    outside the belt. Return widths in meters in the original latitude order.
    """
    latitude = np.asarray(latitude, dtype=np.float64)
    if (latitude.ndim != 1 or latitude.size < 2
            or not np.all(np.isfinite(latitude)) or np.any(np.abs(latitude) > 90)
            or np.unique(latitude).size != latitude.size):
        raise ValueError("Expected distinct global latitude centers in [-90, 90] degrees")
    if (not np.isfinite(latitude_limit) or not 0 < latitude_limit <= 90
            or not np.isfinite(radius) or radius <= 0):
        raise ValueError("Latitude limit must be in (0, 90] and planet radius positive")
    order = np.argsort(latitude)
    centers = latitude[order]
    edges = np.r_[-90.0, 0.5 * (centers[:-1] + centers[1:]), 90.0]
    south = np.maximum(edges[:-1], -latitude_limit)
    north = np.minimum(edges[1:], latitude_limit)
    # A constant-longitude section has dy = a*dphi, with no cosine factor.
    widths = radius * np.deg2rad(np.maximum(north - south, 0.0))
    original_order = np.empty_like(widths)
    original_order[order] = widths
    return original_order


def geopotential_height(p, temperature, Zs=0.0):
    """Integrate upward from the lowest full level, preserving input ordering."""
    p = np.asarray(p, dtype=np.float64)
    temperature = np.asarray(temperature, dtype=np.float64)
    if (p.ndim != 1 or temperature.ndim < 2 or temperature.shape[1] != p.size
            or np.any(~np.isfinite(p)) or np.any(p <= 0)
            or np.unique(p).size != p.size):
        raise ValueError("Expected distinct positive pressures on temperature axis 1")
    order = np.argsort(p)[::-1]
    height = np.empty_like(temperature)
    height[:, order[0], ...] = Zs
    for lower, upper in zip(order[:-1], order[1:]):
        layer_temperature = 0.5 * (
            temperature[:, lower, ...] + temperature[:, upper, ...]
        )
        height[:, upper, ...] = (
            height[:, lower, ...]
            + RD / GRAVITY * layer_temperature * np.log(p[lower] / p[upper])
        )
    return height


def frequency_bands(sample_count, time_step_days, low_period, high_period):
    """Disjoint real-FFT masks covering every nonzero resolved frequency."""
    if (sample_count < 2 or not np.isfinite(time_step_days) or time_step_days <= 0
            or not np.isfinite(low_period) or not np.isfinite(high_period)
            or not 2 * time_step_days < high_period < low_period):
        raise ValueError("Require 2*sample spacing < high-period cutoff < low-period cutoff")
    frequency = np.fft.rfftfreq(sample_count, d=time_step_days)
    masks = (
        (frequency > 0) & (frequency <= 1.0 / low_period),
        (frequency > 1.0 / low_period) & (frequency < 1.0 / high_period),
        frequency >= 1.0 / high_period,
    )
    if any(not mask.any() for mask in masks):
        raise ValueError("The full record must resolve at least one frequency in every band")
    return frequency, masks


def split_time_bands(values, masks):
    """Return full-record mean and real band fields without tapering/detrending."""
    mean = np.mean(values, axis=0)
    coefficients = np.fft.rfft(values - mean, axis=0)
    mask_shape = (coefficients.shape[0],) + (1,) * (values.ndim - 1)
    bands = [
        np.fft.irfft(
            coefficients * mask.reshape(mask_shape), n=values.shape[0], axis=0,
        )
        for mask in masks
    ]
    return mean, bands


def transport_budget(u, s, masks):
    """Calculate all 16 component products averaged over the complete record."""
    if u.shape != s.shape or not (np.all(np.isfinite(u)) and np.all(np.isfinite(s))):
        raise ValueError("Wind and DSE must have equal shapes and finite, regularly sampled data")
    nt = u.shape[0]
    u_mean, u_bands = split_time_bands(u, masks)
    s_mean, s_bands = split_time_bands(s, masks)
    direct = np.einsum("t...,t...->...", u, s) / nt
    terms = np.empty((4, 4, *u.shape[1:]), dtype=np.float64)
    terms[0, 0] = u_mean * s_mean
    for i in range(3):
        terms[0, i + 1] = u_mean * s_bands[i].mean(axis=0)
        terms[i + 1, 0] = u_bands[i].mean(axis=0) * s_mean
        for j in range(3):
            terms[i + 1, j + 1] = (
                np.einsum("t...,t...->...", u_bands[i], s_bands[j]) / nt
            )

    reconstruction_errors = []
    for original, mean, bands in ((u, u_mean, u_bands), (s, s_mean, s_bands)):
        residual = original - mean
        for band in bands:
            residual -= band
        reconstruction_errors.append(float(np.max(np.abs(residual))))

    return {
        "mean_flux": direct, "transport_terms": terms,
        "reconstruction_errors": np.asarray(reconstruction_errors),
    }


def column_transport_contributions(mean_flux, terms, layer_mass):
    """Signed transports in watts, longitude-averaged after meridional integration."""
    # The budget fields already contain sum(F*dy). Apply dp/g once, with no
    # normalization by either meridional extent or total column mass.
    matrix = np.einsum("ijpl,p->ijl", terms, layer_mass).mean(axis=-1)
    total = float(np.einsum("pl,p->l", mean_flux, layer_mass).mean())
    summary = {
        "total": total,
        "mean_flow": float(matrix[0, 0]),
        "low": float(matrix[1, 1]),
        "mid": float(matrix[2, 2]),
        "high": float(matrix[3, 3]),
        "cross": float(sum(
            matrix[i, j] for i in range(4) for j in range(4) if i != j
        )),
        "transient": total - float(matrix[0, 0]),
    }
    reconstructed = sum(summary[name] for name in ("mean_flow", "low", "mid", "high", "cross"))
    if not np.isclose(reconstructed, total, rtol=1e-10, atol=1e-5):
        raise ValueError("Column transport contributions do not close the budget")
    return summary


def main():
    if spinup < 0:
        raise ValueError("Spinup must be nonnegative")
    source_path = DATA_PATH / "atmos_daily.nc"
    with nc.Dataset(source_path) as ds:
        time_variable = ds.variables["time"]
        time = np.asarray(time_variable[spinup:], dtype=np.float64)
        time_units = getattr(time_variable, "units", "days")
        calendar = getattr(time_variable, "calendar", "standard")
        if (time.size < 2 or not time_units.lower().startswith("day")
                or not np.allclose(np.diff(time), 1.0, rtol=0, atol=1e-6)):
            raise ValueError("Expected at least two regularly spaced daily samples after spinup")
        pressure = np.asarray(ds.variables["pfull"][:], dtype=np.float64)
        longitude = np.asarray(ds.variables["lon"][:], dtype=np.float64)
        latitude = np.asarray(ds.variables["lat"][:], dtype=np.float64)
        all_strip_widths = meridional_strip_widths(latitude, LATITUDE_LIMIT, EARTH_RADIUS)
        selected = np.flatnonzero(all_strip_widths > 0)
        if not selected.size:
            raise ValueError("No latitude cells intersect the selected domain")
        frequency, masks = frequency_bands(
            time.size, 1.0, LOW_PERIOD_DAYS, HIGH_PERIOD_DAYS,
        )
        latitude_strip_widths = all_strip_widths[selected]
        if "phalf" not in ds.variables:
            raise ValueError("Half-level pressure is required for column-mass weights")
        half_pressure = np.asarray(ds.variables["phalf"][:], dtype=np.float64)
        if half_pressure.size != pressure.size + 1 or np.any(np.diff(half_pressure) <= 0):
            raise ValueError("Expected increasing half-level pressures bracketing every level")
        layer_mass = np.diff(half_pressure) * 100.0 / GRAVITY

        print(
            f"Full-record filters: {time.size} daily samples; "
            f"all {selected.size} latitudes together; low ≥{LOW_PERIOD_DAYS:g} d, "
            f"middle {HIGH_PERIOD_DAYS:g}–{LOW_PERIOD_DAYS:g} d, "
            f"high ≤{HIGH_PERIOD_DAYS:g} d", flush=True,
        )
        print("Loading wind and temperature for all selected latitudes...", flush=True)
        temperature = np.ma.asarray(
            ds.variables["temp"][spinup:, :, selected, :], dtype=np.float64,
        ).filled(np.nan)
        u = np.ma.asarray(
            ds.variables["ucomp"][spinup:, :, selected, :], dtype=np.float64,
        ).filled(np.nan)
        s = geopotential_height(pressure, temperature)
        s *= GRAVITY
        s += CP * temperature
        del temperature

        print("Decomposing wind and DSE over the complete time axis...", flush=True)
        budget = transport_budget(u, s, masks)
        # Sum the native component products directly over all latitudes.
        terms = np.einsum(
            "ijzyl,y->ijzl", budget["transport_terms"], latitude_strip_widths,
        )
        mean_flux = np.einsum(
            "zyl,y->zl", budget["mean_flux"], latitude_strip_widths,
        )
        errors = budget["reconstruction_errors"]
        del u, s, budget

    pressure_order = np.argsort(pressure)
    longitude_order = np.argsort(longitude % 360.0)
    pressure = pressure[pressure_order]
    longitude = (longitude % 360.0)[longitude_order]
    terms = terms[:, :, pressure_order, :][..., longitude_order]
    mean_flux = mean_flux[pressure_order, :][:, longitude_order]
    layer_mass = layer_mass[pressure_order]
    pressure_weights = layer_mass / layer_mass.sum()
    transient = mean_flux - terms[0, 0]
    closure = mean_flux - terms.sum(axis=(0, 1))
    relative_closure = np.max(np.abs(closure)) / max(np.max(np.abs(mean_flux)), 1.0)
    if relative_closure > 1e-10:
        raise ValueError(f"Transport budget failed closure: relative error {relative_closure:.3e}")
    print(
        f"Field reconstruction max errors: u={errors[0]:.3e} m/s, "
        f"s={errors[1]:.3e} J/kg; relative transport closure={relative_closure:.3e}",
        flush=True,
    )

    prefix = f"q{qflux}_rh{rh}_b{rad}"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / f"{prefix}_frequency_budget.npz"
    metadata = {
        "source_file": str(source_path), "qflux": qflux, "rh": rh, "rad": rad,
        "spinup_records": spinup, "retained_time_samples": time.size,
        "time_start": float(time[0]), "time_end": float(time[-1]),
        "time_units": time_units, "time_calendar": calendar,
        "latitude_min": -LATITUDE_LIMIT,
        "latitude_max": LATITUDE_LIMIT,
        "sampled_latitude_center_min": float(latitude[selected].min()),
        "sampled_latitude_center_max": float(latitude[selected].max()),
        "meridional_reduction": (
            "Sum F_i * dy_i, with dy_i = radius * latitude-cell width in radians; "
            "products calculated before summation; no cosine factor or normalization"
        ),
        "spatial_processing": "All selected latitudes loaded together; direct meridional sum",
        "latitude_cell_bounds": (
            "Midpoints between global latitude centers, polar outer bounds; "
            "clipped to the requested belt, with partially intersecting edge cells retained"
        ),
        "planet_radius_m": EARTH_RADIUS,
        "meridional_extent_m": float(latitude_strip_widths.sum()),
        "frequency_filter": "Full-record real FFT/IRFFT, time axis only, disjoint masks",
        "preprocessing": "Full-record mean removed; no taper, segmenting, or additional detrending",
        "bands": {
            "low": f"0 < f <= 1/{LOW_PERIOD_DAYS:g} day-1",
            "mid": f"1/{LOW_PERIOD_DAYS:g} < f < 1/{HIGH_PERIOD_DAYS:g} day-1",
            "high": f"1/{HIGH_PERIOD_DAYS:g} <= f <= daily Nyquist",
        },
        "frequency_bin_counts": [int(mask.sum()) for mask in masks],
        "transport_terms_dimensions": ["u_component", "s_component", "pfull", "lon"],
        "native_specific_flux_units": "J m kg-1 s-1",
        "transport_units": "J m2 kg-1 s-1",
        "column_transport_units": "W",
        "column_transport_definition": "Sum of <u_i s_j> * dy * dp/g at each longitude",
        "column_mass": "Native half-level pressure differences in Pa divided by gravity",
        "height_reference": "z=0 at the lowest native full model level; dry hydrostatic integration upward",
        "rd_J_kg_K": RD, "cp_J_kg_K": CP, "gravity_m_s2": GRAVITY,
        "closure_max_relative": float(relative_closure),
        "field_reconstruction_max_errors": errors.tolist(),
    }
    # Keep the existing keys used by Zonal_DSE_transport.ipynb, and save
    # the individual diagnostics so visualization needs only this archive.
    summary = column_transport_contributions(mean_flux, terms, layer_mass)
    zonal_mean_flux = mean_flux.mean(axis=-1)
    band_flux = np.stack([terms[index, index] for index in (1, 2, 3)])
    cross_band_flux = sum(
        terms[i, j] for i in (1, 2, 3) for j in (1, 2, 3) if i != j
    )
    band_product_matrix = np.sum(
        terms[1:, 1:].mean(axis=-1) * pressure_weights, axis=-1,
    )
    metadata.update({
        "schema_version": 2,
        "analysis_script": "zonal_energy_transpose.py",
        "visualization_notebook": "Code/plot/Zonal_DSE_transport.ipynb",
        "low_period_days": LOW_PERIOD_DAYS,
        "high_period_days": HIGH_PERIOD_DAYS,
        "sample_interval_days": 1.0,
        "anomaly_definition": (
            "Time-mean meridionally integrated flux minus its longitude mean "
            "at each pressure level; calculated before any cyclic plotting column"
        ),
        "contribution_definition": (
            "Meridional integral, full-column mass integral, then longitude mean; "
            "cross includes all off-diagonal products; transient = total - mean_flow"
        ),
        "component_coordinates": {
            "u_component": "component", "s_component": "component",
            "u_band": "band", "s_band": "band",
        },
        "variables": {},
    })
    flux_units = "J m2 kg-1 s-1"
    output_arrays = {}
    for name, values, dimensions, units, long_name in (
        ("pfull", pressure, ("pfull",), "hPa", "Native full-level pressure, increasing"),
        ("lon", longitude, ("lon",), "degrees_east", "Native longitude, without cyclic duplication"),
        ("latitude", latitude[selected], ("latitude",), "degrees_north", "Selected latitude-cell centers"),
        ("latitude_strip_width_m", latitude_strip_widths, ("latitude",), "m",
         "Meridional cell widths clipped to the selected latitude belt"),
        ("time", time, ("time",), time_units, "Retained daily time coordinate"),
        ("frequency", frequency, ("frequency",), "day-1", "Real-FFT temporal frequencies"),
        ("component", np.array(COMPONENT_NAMES), ("component",), "1", "Wind and DSE component names"),
        ("band", np.array(COMPONENT_NAMES[1:]), ("band",), "1", "Nonzero-frequency band names"),
        ("band_masks", np.stack(masks), ("band", "frequency"), "1",
         "Disjoint nonzero-frequency masks applied to wind and DSE"),
        ("transport_terms", terms, ("u_component", "s_component", "pfull", "lon"), flux_units,
         "All 16 full-record component products, meridionally integrated"),
        ("mean_flux", mean_flux, ("pfull", "lon"), flux_units,
         "Direct time-mean meridionally integrated zonal DSE transport"),
        ("mean_flow_flux", terms[0, 0], ("pfull", "lon"), flux_units,
         "Transport by the full-record mean wind and mean DSE"),
        ("band_flux", band_flux, ("band", "pfull", "lon"), flux_units,
         "Same-band wind-DSE covariance contributions"),
        ("transient_flux", transient, ("pfull", "lon"), flux_units,
         "Total time-mean transport minus mean-flow transport"),
        ("zonal_mean_flux", zonal_mean_flux, ("pfull",), flux_units,
         "Longitude mean of the direct time-mean transport"),
        ("mean_flux_zonal_anomaly", mean_flux - zonal_mean_flux[:, None],
         ("pfull", "lon"), flux_units, "Longitude anomaly of the direct time-mean transport"),
        ("cross_band_flux", cross_band_flux, ("pfull", "lon"), flux_units,
         "Sum of off-diagonal nonzero-frequency band products"),
        ("band_product_matrix", band_product_matrix, ("u_band", "s_band"), flux_units,
         "Column-mass-weighted and longitude mean of meridionally integrated band products"),
        ("closure_residual", closure, ("pfull", "lon"), flux_units,
         "Direct mean transport minus the sum of all 16 component products"),
        ("closure_max_relative", relative_closure, (), "1",
         "Maximum absolute closure residual divided by the maximum absolute direct flux"),
        ("u_reconstruction_max_error", errors[0], (), "m s-1",
         "Maximum absolute wind reconstruction error"),
        ("s_reconstruction_max_error", errors[1], (), "J kg-1",
         "Maximum absolute DSE reconstruction error"),
        ("layer_mass", layer_mass, ("pfull",), "kg m-2", "Native layer mass from half-level dp/g"),
        ("pressure_weights", pressure_weights, ("pfull",), "1",
         "Layer mass normalized by total column mass for the diagnostic matrix"),
        ("column_transport_terms", np.einsum("ijpl,p->ijl", terms, layer_mass),
         ("u_component", "s_component", "lon"), "W",
         "All 16 meridional and column-integrated transport products"),
        ("column_mean_flux", np.einsum("pl,p->l", mean_flux, layer_mass), ("lon",), "W",
         "Direct meridional and column-integrated time-mean transport at each longitude"),
        ("contribution", np.array(CONTRIBUTION_NAMES), ("contribution",), "1",
         "Signed budget contribution names; reference totals overlap their components"),
        ("contribution_transport_W", np.array([summary[name] for name in CONTRIBUTION_NAMES]),
         ("contribution",), "W", "Signed longitude-mean column transport for each bar contribution"),
    ):
        output_arrays[name] = np.asarray(values)
        metadata["variables"][name] = {
            "dimensions": dimensions, "units": units, "long_name": long_name,
        }
    metadata["variables"]["pfull"].update({
        "standard_name": "air_pressure", "positive": "down",
    })
    metadata["variables"]["time"]["calendar"] = calendar
    np.savez_compressed(
        output_path, **output_arrays, metadata=np.array(json.dumps(metadata)),
    )
    print(f"Saved {output_path}", flush=True)


# %%
if __name__ == "__main__":
    main()
