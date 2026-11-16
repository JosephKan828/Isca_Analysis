"""
    Calculate spectral EAPE growth rates for the first two heating EOFs.
    Save spectra, growth rates, and vertical EOF structures to a compressed NPZ file.
    Isca heating is stored in K/s; rates below are calculated in day^-1.
"""

# %% Configure CPU threading before importing numerical libraries
import os
import json

# Override at launch, e.g. EAPE_NUM_THREADS=2 python Code/EAPE_budget.py
NUM_THREADS: int = int(os.environ.get("EAPE_NUM_THREADS", "1"))
if NUM_THREADS < 1:
    raise ValueError("EAPE_NUM_THREADS must be a positive integer")

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
# select parameters
qflux: str = os.environ.get("EAPE_QFLUX", "0")  # intensity of ocean heat uptake
rh   : str = os.environ.get("EAPE_RH", "60") # relative humidity for simple Betts-Miller scheme
rad  : str = os.environ.get("EAPE_RAD", "original")  # radiation setting; ["0", "original", "_lock"]

# cutoff first 1800 days (spin-up period)
spinup: int = 1800

# data path
DATA_PATH: Path = Path(os.environ.get(
    "EAPE_DATA_ROOT", "/data92/b11209013/Isca_Gary/DATA",
)) / f"q{qflux}_rh{rh}_b{rad}"

# load atmospheric data
with nc.Dataset(DATA_PATH / "atmos_daily.nc") as ds:
    time_units = getattr(ds.variables["time"], "units", "days")
    time_calendar = getattr(ds.variables["time"], "calendar", "standard")

    coords: dict[str, np.ndarray] = {
        key: ds.variables[key][:]
        for key in ds.variables.keys()
        if key in ds.dimensions.keys()  # include only coordinate variables
    }

    # Select the tropical troposphere before loading the fields (aquaplanet).
    lat_mask = np.abs(coords["lat"]) <= 5.0
    level_mask = (coords["pfull"] >= 100.0) & (coords["pfull"] <= 1000.0)

    data : dict[str, np.ndarray] = {
        key: np.asarray(ds.variables[key][spinup:, level_mask, lat_mask, :], dtype=np.float64)
        for key in [
            "temp", "dt_tg_convection", "dt_tg_condensation", "tdt_rad", "dt_tg_diffusion"
        ]
    } # shape: (time, level, lat, lon)

coords["time"] = coords["time"][spinup:]
coords["lat"] = coords["lat"][lat_mask]
coords["pfull"] = coords["pfull"][level_mask]
# Hydrostatic column mass is dp/g; use increasing pressure in Pa.
pressure_order = np.argsort(coords["pfull"])
pressure_pa = np.asarray(coords["pfull"])[pressure_order] * 100.0
gravity = 9.80665  # m/s^2
nt, nz, ny, nx = data["temp"].shape
if nt < 96 or not np.allclose(np.diff(coords["time"]), 1.0):
    raise ValueError("At least 96 regularly spaced daily samples are required after spinup.")

# Sum total heating and convert K/s to K/day; tdt_rad already includes solar heating.
data["tdt_diabatic"] = 86400.0 * (
    data["dt_tg_convection"] + data["dt_tg_condensation"]
    + data["tdt_rad"] + data["dt_tg_diffusion"]
)

# %% Processing data for EOF decomposition
# The model has no seasonal cycle; subtract only the time mean.
diabatic_anomaly: np.ndarray = np.ascontiguousarray(
    data["tdt_diabatic"] - np.nanmean(data["tdt_diabatic"], axis=0)
)
temp_anomaly: np.ndarray = np.ascontiguousarray(
    data["temp"] - np.nanmean(data["temp"], axis=0)
)

# transforming data into 2D array
diabatic_reshape: np.ndarray = diabatic_anomaly.transpose(1, 0, 2, 3).reshape(nz, nt * ny * nx)

# calculate EOF
U, s, VT = np.linalg.svd(diabatic_reshape, full_matrices=False)

# Retain the singular values so heating PCs have physical K/day amplitudes.
diabatic_mode1: np.ndarray = (s[0] * VT[0, :]).reshape(nt, ny, nx)
diabatic_mode2: np.ndarray = (s[1] * VT[1, :]).reshape(nt, ny, nx)

# calculate PCs of temperature perturbation with the same heating EOF basis
temp_reshape: np.ndarray = temp_anomaly.transpose(1, 0, 2, 3).reshape(nz, nt * ny * nx)

# project temperature anomaly onto diabatic EOFs
temp_mode1: np.ndarray = np.dot(U[:, 0], temp_reshape).reshape(nt, ny, nx)
temp_mode2: np.ndarray = np.dot(U[:, 1], temp_reshape).reshape(nt, ny, nx)

# %% Chunk data and taper only the five-day edges of each segment
CHUNK_SIZE: int = 96
OVERLAP   : int = 60
TAPER_DAYS: int = 5
window: np.ndarray = np.ones(CHUNK_SIZE)
edge = 0.5 * (1 - np.cos(np.pi * np.arange(TAPER_DAYS + 1) / TAPER_DAYS))
window[:TAPER_DAYS] = edge[:-1]
window[-TAPER_DAYS:] = edge[-2::-1]
window_energy = np.mean(window ** 2)

# chunking
Q1_1_chunks: list[np.ndarray] = [
    diabatic_mode1[i:i+CHUNK_SIZE, ...] * window[:, None, None]
    for i in range(0, nt - CHUNK_SIZE + 1, CHUNK_SIZE - OVERLAP)
]
Q1_2_chunks: list[np.ndarray] = [
    diabatic_mode2[i:i+CHUNK_SIZE, ...] * window[:, None, None]
    for i in range(0, nt - CHUNK_SIZE + 1, CHUNK_SIZE - OVERLAP)
]

T_1_chunks: list[np.ndarray] = [
    temp_mode1[i:i+CHUNK_SIZE, ...] * window[:, None, None]
    for i in range(0, nt - CHUNK_SIZE + 1, CHUNK_SIZE - OVERLAP)
]

T_2_chunks: list[np.ndarray] = [
    temp_mode2[i:i+CHUNK_SIZE, ...] * window[:, None, None]
    for i in range(0, nt - CHUNK_SIZE + 1, CHUNK_SIZE - OVERLAP)
]

T_chunks: list[np.ndarray] = [
    temp_anomaly[i:i+CHUNK_SIZE, ...] * window[:, None, None, None]
    for i in range(0, nt - CHUNK_SIZE + 1, CHUNK_SIZE - OVERLAP)
]

# %% 
def _calc_fft2(
        data: np.ndarray,
        lon_axis: int,
        time_axis: int
) -> np.ndarray:

    # perform FFT along time axis
    # and perform IFFT along the longitude axis
    data_fft: np.ndarray = np.fft.fft(data, axis=time_axis) / data.shape[time_axis]
    data_ifft: np.ndarray = np.fft.ifft(data_fft, axis=lon_axis)

    return data_ifft


def _smooth_121(
        data: np.ndarray,
        axes: tuple[int, ...] = (0, 1)
) -> np.ndarray:
    """Apply a 1-2-1 filter along each selected axis with zero-padded edges."""
    def smooth_line(line: np.ndarray) -> np.ndarray:
        # Trim the full convolution to preserve even one- or two-element axes.
        return np.convolve(line, [0.25, 0.5, 0.25], mode="full")[1:-1]

    smoothed = data
    for axis in axes:
        smoothed = np.apply_along_axis(smooth_line, axis=axis, arr=smoothed)
    return smoothed


cov_Q1T_1_list: list[np.ndarray] = []
cov_Q1T_2_list: list[np.ndarray] = []
cov_TT_1_list : list[np.ndarray] = []
cov_TT_2_list : list[np.ndarray] = []
cov_TT_list   : list[np.ndarray] = []

# calculate covariance matrix of Q1/T1 and Q2/T2 for each chunk
for i in range(len(Q1_1_chunks)):
    Q1_1_chunk: np.ndarray = Q1_1_chunks[i]
    T_1_chunk : np.ndarray = T_1_chunks[i]
    Q1_2_chunk: np.ndarray = Q1_2_chunks[i]
    T_2_chunk : np.ndarray = T_2_chunks[i]
    T_chunk   : np.ndarray = T_chunks[i]

    # apply FFT to chunked data
    Q1_1_FFT: np.ndarray = _calc_fft2(Q1_1_chunk, lon_axis=2, time_axis=0)
    T_1_FFT : np.ndarray = _calc_fft2(T_1_chunk, lon_axis=2, time_axis=0)
    Q1_2_FFT: np.ndarray = _calc_fft2(Q1_2_chunk, lon_axis=2, time_axis=0)
    T_2_FFT : np.ndarray = _calc_fft2(T_2_chunk, lon_axis=2, time_axis=0)

    T_FFT: np.ndarray = _calc_fft2(T_chunk, lon_axis=-1, time_axis=0)

    # calculate the covariance spectrum and average over latitude
    cov_Q1T_1: np.ndarray = np.nanmean(Q1_1_FFT * np.conj(T_1_FFT), axis=1)
    cov_Q1T_2: np.ndarray = np.nanmean(Q1_2_FFT * np.conj(T_2_FFT), axis=1)
    cov_TT_1 : np.ndarray = np.nanmean(T_1_FFT * np.conj(T_1_FFT), axis=1)
    cov_TT_2 : np.ndarray = np.nanmean(T_2_FFT * np.conj(T_2_FFT), axis=1)
    # Full temperature power retains the pressure axis: (frequency, level, wavenumber).
    cov_TT   : np.ndarray = np.nanmean(np.abs(T_FFT) ** 2, axis=-2)

    # apply a 1-2-1 filter along frequency (axis 0) and zonal wavenumber (axis 1)
    cov_Q1T_1: np.ndarray = _smooth_121(np.fft.fftshift(cov_Q1T_1), axes=(0, 1))
    cov_Q1T_2: np.ndarray = _smooth_121(np.fft.fftshift(cov_Q1T_2), axes=(0, 1))
    cov_TT_1 : np.ndarray = _smooth_121(np.fft.fftshift(cov_TT_1), axes=(0, 1))
    cov_TT_2 : np.ndarray = _smooth_121(np.fft.fftshift(cov_TT_2), axes=(0, 1))
    # Shift only frequency and wavenumber; preserve the physical level order.
    cov_TT   : np.ndarray = _smooth_121(
        np.fft.fftshift(cov_TT, axes=(0, -1)), axes=(0, -1),
    )

    # Integrate power, not Fourier amplitudes: integral P_T dp/g (K^2 kg/m^2).
    cov_TT_intg: np.ndarray = np.trapz(
        cov_TT[:, pressure_order, :], x=pressure_pa, axis=1,
    ) / gravity

    # store to a list for later use
    cov_Q1T_1_list.append(cov_Q1T_1)
    cov_Q1T_2_list.append(cov_Q1T_2)
    cov_TT_1_list.append(cov_TT_1)
    cov_TT_2_list.append(cov_TT_2)
    cov_TT_list.append(cov_TT_intg)

# Average segments and correct their common taper energy. This cancels in sigma.
cov_Q1T_1_avg: np.ndarray = np.mean(np.array(cov_Q1T_1_list), axis=0) / window_energy
cov_Q1T_2_avg: np.ndarray = np.mean(np.array(cov_Q1T_2_list), axis=0) / window_energy
cov_TT_1_avg: np.ndarray = np.mean(np.array(cov_TT_1_list), axis=0) / window_energy
cov_TT_2_avg: np.ndarray = np.mean(np.array(cov_TT_2_list), axis=0) / window_energy
cov_TT_avg: np.ndarray = np.mean(np.array(cov_TT_list), axis=0) / window_energy

# %% Calculate growth rates (day^-1)
sigma = []
for QT, TT in (
    (cov_Q1T_1_avg, cov_TT_1_avg),
    (cov_Q1T_2_avg, cov_TT_2_avg),
):
    # Exclude unresolved temperature power to avoid division by numerical zero.
    valid_T = TT.real > 1e-12 * TT.real.max()
    growth = np.full(TT.shape, np.nan)
    np.divide(2 * QT.real, TT.real, out=growth, where=valid_T)
    sigma.append(growth)
sigma_1, sigma_2 = sigma

# %% Save numerical diagnostics for separate visualization
frequency = np.fft.fftshift(np.fft.fftfreq(CHUNK_SIZE, d=1.0))
wavenumber = np.fft.fftshift(np.fft.fftfreq(nx, d=1.0 / nx))
temperature_power = cov_TT_avg.real
explained_variance = 100.0 * s[:2] ** 2 / np.sum(s ** 2)

script_dir = (
    Path(__file__).resolve().parent if "__file__" in globals()
    else Path("/home/b11209013/Isca_analysis/Code")  # IDE cell execution
)
output_dir = Path(os.environ.get(
    "EAPE_OUTPUT_DIR", str(script_dir.parent / "Output" / "EAPE_budget"),
))
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / f"q{qflux}_rh{rh}_b{rad}_eape.npz"

metadata = {
    "title": "Spectral diabatic growth rates and vertical heating EOFs",
    "analysis_script": "EAPE_budget.py",
    "source_file": str(DATA_PATH / "atmos_daily.nc"),
    "qflux": qflux, "rh": rh, "rad": rad,
    "spinup_records": spinup,
    "retained_time_samples": nt,
    "latitude_samples": ny, "longitude_samples": nx,
    "latitude_min": float(np.min(coords["lat"])),
    "latitude_max": float(np.max(coords["lat"])),
    "pressure_min_hpa": float(pressure_pa[0] / 100.0),
    "pressure_max_hpa": float(pressure_pa[-1] / 100.0),
    "time_start": float(coords["time"][0]),
    "time_end": float(coords["time"][-1]),
    "time_units": time_units, "time_calendar": time_calendar,
    "segment_length_days": CHUNK_SIZE,
    "segment_overlap_days": OVERLAP,
    "segment_count": len(Q1_1_chunks),
    "edge_taper_days": TAPER_DAYS,
    "window_energy": float(window_energy),
    "temperature_power_relative_cutoff": 1e-12,
    "gravity_m_s2": gravity,
    "heating_components": "dt_tg_convection dt_tg_condensation tdt_rad dt_tg_diffusion",
    "anomalies": "Full-record time mean removed independently at each grid point",
    "eof_metric": "Unweighted Euclidean inner product over native pressure levels",
    "fft_convention": "Time FFT / segment length, longitude IFFT; positive k at positive f is eastward",
    "spectral_normalization": "Two-sided power per discrete bin; segment mean divided by mean(window**2)",
    "spectral_smoothing": "One 1-2-1 pass each along frequency and wavenumber, with zero-padded edges",
    "cross_spectrum_definition": "heating_PC_hat * conjugate(temperature_PC_hat)",
    "growth_rate_definition": "2 * real(mean cross-spectrum) / mean temperature-PC power; units day-1",
    "column_integration": "Trapezoidal integral of full-field temperature power over increasing pressure in Pa, divided by gravity",
}

# Named numerical arrays and JSON metadata can be loaded with allow_pickle=False.
output_arrays: dict[str, np.ndarray] = {}
metadata["variables"] = {}
for name, values, dimensions, units, long_name in (
    ("frequency", frequency, ("frequency",), "day-1", "Temporal frequency in cycles per day"),
    ("wavenumber", wavenumber, ("wavenumber",), "1", "Global zonal wavenumber; positive is eastward at positive frequency"),
    ("pfull", pressure_pa / 100.0, ("pfull",), "hPa", "Retained native full-level pressure in increasing order"),
    ("mode", np.array([1, 2], dtype=np.int32), ("mode",), "1", "Heating EOF number"),
    ("temperature_power", temperature_power, ("frequency", "wavenumber"),
     "K2 kg m-2", "Column mass-integrated full-field temperature power per spectral bin"),
    ("sigma_1", sigma_1, ("frequency", "wavenumber"),
     "day-1", "Heating contribution to EOF 1 energy/variance growth rate"),
    ("sigma_2", sigma_2, ("frequency", "wavenumber"),
     "day-1", "Heating contribution to EOF 2 energy/variance growth rate"),
    ("heating_eof", U[pressure_order, :2].T, ("mode", "pfull"),
     "1", "Unit-norm vertical heating EOF; sign is arbitrary"),
    ("explained_variance", explained_variance, ("mode",),
     "%", "Percentage of total unweighted heating variance explained"),
    ("modal_temperature_power", np.stack((cov_TT_1_avg.real, cov_TT_2_avg.real)),
     ("mode", "frequency", "wavenumber"), "K2", "Temperature-PC power per spectral bin"),
    ("heating_temperature_cross_spectrum", np.stack((cov_Q1T_1_avg, cov_Q1T_2_avg)),
     ("mode", "frequency", "wavenumber"), "K2 day-1", "Complex heating-temperature PC cross-spectrum"),
):
    output_arrays[name] = np.asarray(values)
    metadata["variables"][name] = {
        "dimensions": dimensions, "units": units, "long_name": long_name,
    }
metadata["variables"]["pfull"].update({"standard_name": "air_pressure", "positive": "down"})

np.savez_compressed(output_path, **output_arrays, metadata=np.array(json.dumps(metadata)))
print(f"Saved {output_path}", flush=True)
# %%
