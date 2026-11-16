# Isca heating and spectral EAPE analysis

This guide combines the spectral EAPE procedure, its algorithm review, and the
vertical heating-variance documentation. Both calculation scripts use
sequential `# %%` workflows and save numerical diagnostics as compressed NumPy NPZ files.
Create visualizations in a separate script using the saved coordinates and units.

| Script | Diagnostic | Output directory |
|---|---|---|
| [EAPE_budget.py](EAPE_budget.py) | Heating contribution to EAPE growth for the first two heating EOFs, plus their vertical structures | `Output/EAPE_budget` |
| [Heating_variance.py](Heating_variance.py) | Separate convective and condensational heating-variance profiles | `Output/Heating_variance` |

The growth-rate definition follows
[Chien and Kim (2023), equations (10)–(12)](https://doi.org/10.1175/JAS-D-22-0067.1)
([accessible full text](https://www.osti.gov/servlets/purl/2420475)).
The review below distinguishes the current diagnostic from proposed scientific
and robustness improvements.

- [Spectral EAPE growth rates](#spectral-eape-growth-rates)
- [Running the EAPE analysis](#running-the-eape-analysis)
- [Vertical heating variance](#vertical-heating-variance)
- [Reading the saved results](#reading-the-saved-results)
- [Diagnostic review and proposed improvements](#diagnostic-review-and-proposed-improvements)
- [Historical validation and observations](#historical-validation-and-observations)

## Spectral EAPE growth rates

### Configuration

The script reads `qflux`, `rh`, and `rad` from `EAPE_QFLUX`, `EAPE_RH`,
and `EAPE_RAD`. Other calculation defaults are edited directly in the script.
The batch runner sets all three experiment parameters for each run.

| Setting | Current value |
|---|---|
| Ocean heat-uptake parameter, `qflux` | `"0"` |
| Betts–Miller relative-humidity parameter, `rh` | `"60"` |
| Radiation setting, `rad` | `"original"` |
| Spin-up discarded | First 1800 daily records |
| Latitude range | 5°S–5°N, inclusive |
| Pressure range | 100–1000 hPa, inclusive, using native model levels |
| Segment length | 96 days |
| Segment overlap | 60 days; successive segments start 36 days apart |
| Edge taper | Five-day half-Hann taper at each end |
| Numerical temperature-power cutoff | `1e-12` times each mode's maximum power |

With the default experiment parameters, the selected input is:

```text
/data92/b11209013/Isca_Gary/DATA/q0_rh60_boriginal/atmos_daily.nc
```

The script assumes an aquaplanet, daily time coordinates measured in days,
pressure coordinates in hPa, and a uniform global longitude grid. Column
mass integration uses `gravity = 9.80665` m/s².

### 1. Load the selected domain

Discard the spin-up records and select tropical latitudes and tropospheric
levels before loading the fields. Arrays have axes:

```text
(time, level, latitude, longitude)
```

The fields are temperature (`temp`) and four temperature tendencies:
`dt_tg_convection`, `dt_tg_condensation`, `tdt_rad`, and `dt_tg_diffusion`.
Data are read as 64-bit floating-point arrays.

### 2. Construct total diabatic heating

Sum convection, condensation, radiation, and diffusion tendencies. Convert
their stored units of K/s to K/day:

$$
Q = 86400\left(Q_{\mathrm{conv}} + Q_{\mathrm{cond}}
                  + Q_{\mathrm{rad}} + Q_{\mathrm{diff}}\right).
$$

`tdt_rad` already includes solar heating, so `tdt_solar` is not added separately.
Temperature remains in K.

### 3. Calculate anomalies

Subtract each field's time mean independently at every level and grid point:

$$
Q' = Q - \overline{Q}^{\,t}, \qquad
T' = T - \overline{T}^{\,t}.
$$

The means use the retained records after spin-up. The model has no seasonal
cycle, so no annual harmonics or seasonal-cycle fitting are applied.

### 4. Obtain a shared vertical EOF basis

Reshape heating anomalies into a matrix with levels as rows and all
time–latitude–longitude samples as columns:

$$
A_Q \in \mathbb{R}^{N_z\times(N_tN_yN_x)}, \qquad
A_Q = U\,\operatorname{diag}(s)\,V^\mathsf{T}.
$$

The first two columns of `U` are the leading heating EOFs. The SVD uses an
ordinary Euclidean inner product: there is no pressure or latitude weighting.
The native pressure levels are retained without interpolation. Inspect the
EOF structures when interpreting them as baroclinic modes.

### 5. Calculate heating and temperature PCs

For each mode $i=1,2$, the heating PC is

$$
q_i = s_i V_i^\mathsf{T} = U_i^\mathsf{T}A_Q.
$$

Including the singular value preserves the heating amplitude in K/day.
Project temperature anomalies onto the **same heating EOFs**:

$$
t_i = U_i^\mathsf{T}A_T.
$$

Reshape both PCs to `(time, latitude, longitude)`. Using the same vertical basis
makes their spectral covariance comparable. These PCs describe each mode's
amplitude; they are not full height-dependent reconstructed fields.

### 6. Segment, taper, and transform

Use complete 96-day segments starting at retained-record indices 0, 36, 72,
and so on. Discard any incomplete final segment. Apply the same five-day
half-Hann edge taper to heating and temperature; the segment interior has
weight one.

For each latitude, `_calc_fft2` applies:

1. A time FFT divided by the segment length.
2. A longitude IFFT, which includes the longitude-length normalization.

With increasing eastward longitude, positive wavenumber at positive frequency
denotes eastward propagation. Both fields use identical transforms and tapers,
so their common normalization cancels in the growth-rate ratio.

### 7. Calculate and average spectral products

For each segment and mode, calculate the complex cross-spectrum and the
temperature PC power spectrum, then take their arithmetic latitude means:

$$
C_{QT,i}(f,k) = \left\langle\widehat q_i\widehat t_i^*\right\rangle_y,
\qquad
P_{TT,i}(f,k) = \left\langle\widehat t_i\widehat t_i^*\right\rangle_y.
$$

The star denotes complex conjugation. Shift zero frequency and zero wavenumber
to the centers of their arrays. Apply one 1–2–1 smoothing pass along frequency
and one along wavenumber, using weights `[0.25, 0.5, 0.25]` and zero-padded
edges. Average the smoothed products across segments.

For the complete temperature anomaly field, also compute power at every
pressure level from a four-dimensional segment `(time, level, latitude,
longitude)`. Its taper has shape `(96, 1, 1, 1)`. After the latitude mean, this
spectrum has shape `(frequency, level, wavenumber)`. Shift and smooth only the
frequency and wavenumber axes, preserving the pressure-level order.

Divide all segment-averaged power and cross-spectra by `mean(window**2)` to
correct the common taper energy. This normalization cancels in the growth-rate
ratio. Spectra remain two-sided power or cross-power per discrete spectral bin;
there is no positive-frequency doubling or spectral-density conversion.

The script averages the **products before taking their ratio**. It does not
average segment-by-segment growth rates. Heating power and coherence squared
are not calculated.

### 8. Vertically integrate total temperature power

Use power from the **full temperature anomaly field**, including all vertical
modes, rather than summing only the first two temperature PC powers. Let
$P_T(f,p,k)$ denote its latitude/segment-averaged power. The requested column
mass integral is

$$
P_{T,\mathrm{column}}(f,k)
= \frac{1}{g}\int_{p_{\mathrm{top}}}^{p_{\mathrm{bottom}}} P_T(f,p,k)\,dp.
$$

Pressure is converted from hPa to Pa. Apply trapezoidal integration along the
level axis, using increasing pressure. Hydrostatic column mass is `dp/g`;
there is no extra multiplication by `p` and no division by the layer depth.
The units are **K²·kg/m²** per spectral bin. This is an integrated variance
diagnostic, not a normalized layer mean or an EAPE energy integral.

The actual limits are the shallowest and deepest retained native levels within
100–1000 hPa. The code does not extrapolate to the requested selection bounds.
The result for each segment is `cov_TT_intg`; the window-corrected segment mean
is `cov_TT_avg`.

### 9. Calculate the growth rate

For each mode, use the averaged spectra to calculate

$$
\sigma_i(f,k) = 2\frac{\operatorname{Re}\left[C_{QT,i}(f,k)\right]}
                         {\operatorname{Re}\left[P_{TT,i}(f,k)\right]}.
$$

Heating PCs are in K/day and temperature PCs are in K, giving growth rates in
day⁻¹. Positive values indicate EAPE generation by heating; negative values
indicate damping. The real part selects the in-phase or opposite-phase
contribution. This is a normalized heating contribution; interpreting it as
total wave growth also requires considering the other energy-budget terms.

The factor of two describes an **energy/variance growth rate**. For a scalar
amplitude satisfying `dT/dt = a T`, this diagnostic returns `2a`, while its
amplitude grows at rate `a`. Interpret the result as the heating contribution
to growth in that modal temperature metric. The calculation does not diagnose
the complete energy budget or independently estimate a wave's observed
exponential growth. No change to the factor of two is suggested.

Where temperature power is at most `1e-12` times that mode's maximum power,
the result is `NaN` to avoid division by numerical zero. This is a numerical
check, **not a significance test**. No coherence-based or other statistical
significance masking is applied.

### 10. Save the numerical diagnostics

Each experiment writes one compressed NPZ archive with double-precision diagnostics:

```text
Output/EAPE_budget/q{qflux}_rh{rh}_b{rad}_eape.npz
```

Set `EAPE_OUTPUT_DIR` to change its directory. The directory is created
automatically, and rerunning an experiment replaces its corresponding file.
The script prints the saved path and does not create figures.

| Variable | Dimensions | Units | Meaning |
|---|---|---|---|
| `frequency` | `frequency` | day⁻¹ | Cycles per day, in shifted FFT order |
| `wavenumber` | `wavenumber` | Dimensionless | Global zonal wavenumber, in shifted FFT order |
| `pfull` | `pfull` | hPa | Retained native levels, sorted by increasing pressure |
| `mode` | `mode` | EOF number | Values 1 and 2 |
| `temperature_power` | `frequency, wavenumber` | K²·kg/m² | Full-field, column mass-integrated temperature power per spectral bin; corresponds to `cov_TT_avg.real` |
| `sigma_1`, `sigma_2` | `frequency, wavenumber` | day⁻¹ | Heating contribution to each mode's energy/variance growth rate |
| `heating_eof` | `mode, pfull` | Dimensionless | Unit-norm vertical heating EOFs, sorted with pressure |
| `explained_variance` | `mode` | % | Fraction of total unweighted heating variance explained, expressed as a percentage |
| `modal_temperature_power` | `mode, frequency, wavenumber` | K² | Averaged temperature-PC power per spectral bin |
| `heating_temperature_cross_spectrum` | `mode, frequency, wavenumber` | K²/day | Complex averaged heating-temperature PC cross-spectrum |
| `metadata` | Scalar Unicode string | JSON | Experiment settings and per-array dimensions, units, and descriptions |

Both positive and negative frequencies, DC, and all zonal wavenumbers are
saved at the calculated bin resolution. Spectra are not clipped to a display
range or log-transformed. The existing numerical low-temperature-power cutoff
is preserved as missing (`NaN`) growth rates. Cross-spectra retain NumPy's
native complex dtype; use `.real` and `.imag` to access their components.

The `metadata` entry is a JSON string recording the input file, experiment parameters, actual
latitude/pressure bounds, time range and calendar, discarded spin-up, sample
and segment counts, taper settings, FFT convention, smoothing, EOF metric,
integration method, and numerical cutoff. All saved spectra include the
existing taper-energy correction. Numerical arrays are stored by name;
metadata uses a Unicode string rather than an object array, so the archive
can be loaded with `np.load(path, allow_pickle=False)`. Units and axis names
are available in `metadata["variables"][array_name]` after decoding the JSON.

### Visualization in a separate script

For a view matching the earlier spectral figures, select strictly positive
frequencies and zonal wavenumbers −15 to 15. At positive frequency, positive
wavenumber denotes eastward propagation. Frequency spacing is `1/96` cycles/day;
the even-length FFT places the Nyquist bin at −0.5 cycles/day. Spectra have
shape `(frequency, wavenumber)` and can be passed to `pcolormesh` without
transposing.

Use a sequential palette for `log10(temperature_power)` after masking invalid
or nonpositive values. Its units before transformation are K²·kg/m² per
spectral bin. Use a diverging palette centered on zero for `sigma_1` and
`sigma_2`, with day⁻¹ units: negative values indicate damping and positive
values indicate generation. The earlier shared ±10 day⁻¹ display range can
be reused for comparisons, with colorbar extensions for larger magnitudes.
Choose a common power scale when comparing experiments.

Plot `heating_eof[0]` and `heating_eof[1]` against `pfull`
with pressure increasing downward. Their signs are arbitrary, and
`explained_variance` describes heating variance rather than temperature
variance or EAPE. The earlier EOF amplitude display range was ±0.4.
The saved pressure coordinate and EOF profiles are already sorted together.

## Running the EAPE analysis

From `/home/b11209013/Isca_analysis`, use a Python environment with NumPy
and netCDF4. Run the default experiment with:

```bash
python Code/EAPE_budget.py
```

Override the experiment parameters at launch, for example:

```bash
EAPE_QFLUX=25 EAPE_RH=70 EAPE_RAD=_lock python Code/EAPE_budget.py
```

| Environment variable | Default | Purpose |
|---|---|---|
| `EAPE_QFLUX` | `0` | Ocean heat-uptake parameter for a direct Python run |
| `EAPE_RH` | `60` | Betts–Miller relative-humidity parameter for a direct Python run |
| `EAPE_RAD` | `original` | Radiation setting for a direct Python run |
| `EAPE_DATA_ROOT` | `/data92/b11209013/Isca_Gary/DATA` | Parent of the experiment directories |
| `EAPE_OUTPUT_DIR` | Project's `Output/EAPE_budget` | NPZ output directory |
| `EAPE_NUM_THREADS` | `1` | Numerical-library thread limit; must be positive |
| `PYTHON` | `python` | Python interpreter used by the shell runner |

To run all **18 combinations** of `qflux = 0, 25, 50`, `rh = 60, 70, 80`,
and `rad = original, _lock`:

```bash
bash Code/run_EAPE_budget.sh
```

The [batch runner](run_EAPE_budget.sh) sets `EAPE_QFLUX`, `EAPE_RH`, and
`EAPE_RAD` for each combination, runs one analysis at a time, and stops if any
analysis fails. Its experiment lists are defined in the shell script; edit
those loops to change the batch selection. Each run reads
`q{qflux}_rh{rh}_b{rad}/atmos_daily.nc` and writes its own NPZ archive
containing the spectra and vertical EOFs.

The default interpreter is `python` from the active environment. Override it
when needed:

```bash
PYTHON=/path/to/analysis-environment/bin/python bash Code/run_EAPE_budget.sh
```

Alternatively, execute the `# %%` cells in order in the IDE. Numerical libraries
use one CPU thread by default; set `EAPE_NUM_THREADS` before launching to change
that value.

## Vertical heating variance

[Heating_variance.py](Heating_variance.py) calculates separate variance profiles
for `dt_tg_convection` and `dt_tg_condensation`, following the sequential `# %%`
coding style used by `EAPE_budget.py`.

The default selection is 5°S–5°N and 100–1000 hPa after discarding the first
1800 daily records. Heating is converted from K/s to K/day before calculating
population variance (`ddof=0`) across all time, latitude, and longitude samples
at each pressure level, using `np.nanvar(data, axis=(0, 2, 3), ddof=0)`.
Each valid sample has equal weight; missing values are excluded. Both temporal
variability and spatial differences in time-mean heating contribute to this
diagnostic.

Each experiment writes one compressed NPZ archive with double-precision diagnostics:

```text
Output/Heating_variance/q{qflux}_rh{rh}_b{rad}_variance.npz
```

| Variable | Dimensions | Units | Meaning |
|---|---|---|---|
| `pfull` | `pfull` | hPa | Retained native pressure levels in increasing order |
| `convection_variance` | `pfull` | K²/day² | Population variance of `dt_tg_convection` |
| `condensation_variance` | `pfull` | K²/day² | Population variance of `dt_tg_condensation` |
| `metadata` | Scalar Unicode string | JSON | Experiment settings and per-array dimensions, units, and descriptions |

The profiles are sorted together with pressure. The JSON `metadata` entry records the
source file, experiment, selected domain, time range/calendar, retained sample
counts, spin-up, variance convention, and conversion from K/s to K/day.
Set `HEATING_OUTPUT_DIR` to change the output directory. It is created
automatically, and rerunning an experiment replaces its corresponding file.

For visualization, plot the profiles against pressure increasing downward.
The earlier figure used independent lower and upper variance axes for
convection and condensation because the component magnitudes differ.

### Running the heating-variance analysis

Use a Python environment with NumPy and netCDF4. From
`/home/b11209013/Isca_analysis`, run the default experiment with:

```bash
python Code/Heating_variance.py
```

Run all **18 combinations** of `qflux = 0, 25, 50`, `rh = 60, 70, 80`,
and `rad = original, _lock` with:

```bash
bash Code/run_heating_variance.sh
```

Inputs are defined near the top of the [shell runner](run_heating_variance.sh)
and can also be overridden at launch. For example, run only qflux 25 with RH 70
and the original radiation setting:

```bash
HEATING_QFLUX_VALUES="25" HEATING_RH_VALUES="70" HEATING_RAD_VALUES="original" \
    bash Code/run_heating_variance.sh
```

| Environment variable | Default | Purpose |
|---|---|---|
| `HEATING_QFLUX_VALUES` | `0 25 50` | Space-separated qflux inputs for the shell runner |
| `HEATING_RH_VALUES` | `60 70 80` | Space-separated RH inputs for the shell runner |
| `HEATING_RAD_VALUES` | `original _lock` | Space-separated radiation settings for the shell runner; takes precedence over `HEATING_RAD` |
| `HEATING_QFLUX` | `0` | qflux for a direct Python run |
| `HEATING_RH` | `60` | RH for a direct Python run |
| `HEATING_RAD` | `original` for a direct Python run | Radiation setting; also supplies the shell runner's radiation selection when `HEATING_RAD_VALUES` is unset; without either variable, the runner uses `original _lock` |
| `HEATING_SPINUP` | `1800` | Number of initial time records discarded; must be nonnegative |
| `HEATING_DATA_ROOT` | `/data92/b11209013/Isca_Gary/DATA` | Parent of the experiment directories |
| `HEATING_OUTPUT_DIR` | Project's `Output/Heating_variance` | NPZ output directory |
| `HEATING_NUM_THREADS` | `1` | Numerical-library thread limit; must be positive |
| `PYTHON` | `python` | Python interpreter used by the shell runner |

When running the Python script directly, use `HEATING_QFLUX`, `HEATING_RH`,
and `HEATING_RAD` for a single experiment. The runner sets those variables for
each selected combination, runs sequentially, and stops if any analysis fails.

## Reading the saved results

The calculation scripts require NumPy and netCDF4 to read the original model
input. Reading the saved `.npz` outputs requires only NumPy; the standard-library
`json` module decodes the optional metadata. Each archive holds multiple named
arrays, which can be listed with `saved.files`.

From the project directory, a separate visualization script can load them with:

```python
import json
import numpy as np

with np.load("Output/EAPE_budget/q0_rh60_boriginal_eape.npz", allow_pickle=False) as saved:
    frequency = saved["frequency"]
    wavenumber = saved["wavenumber"]
    display = np.ix_(frequency > 0, np.abs(wavenumber) <= 15)
    sigma_1 = saved["sigma_1"][display]
    sigma_2 = saved["sigma_2"][display]
    temperature_power = saved["temperature_power"][display]
    plot_frequency = frequency[frequency > 0]
    plot_wavenumber = wavenumber[np.abs(wavenumber) <= 15]
    heating_eof_1 = saved["heating_eof"][0]
    eof_pressure = saved["pfull"]
    explained_variance = saved["explained_variance"]
    eape_metadata = json.loads(saved["metadata"].item())

with np.load("Output/Heating_variance/q0_rh60_boriginal_variance.npz", allow_pickle=False) as saved:
    pressure = saved["pfull"]
    convection_variance = saved["convection_variance"]
    condensation_variance = saved["condensation_variance"]
    variance_metadata = json.loads(saved["metadata"].item())

sigma_units = eape_metadata["variables"]["sigma_1"]["units"]  # day-1
```

All numerical arrays remain usable after the archive is closed. Missing growth
rates remain ordinary NumPy NaNs. `heating_eof[0]` and `heating_eof[1]` correspond
to the saved `mode` values 1 and 2; NumPy indexing starts at zero. The spectral
array order is `(frequency, wavenumber)`, and the per-mode cross-spectrum and
temperature power use `(mode, frequency, wavenumber)`.

## Diagnostic review and proposed improvements

These scientific and robustness suggestions remain proposals; the current
implementation described above does not apply them. The shared heating EOF
basis, averaging of products before division, and separate frequency/wavenumber
smoothing are consistent with the modal spectral definition in equation (12)
of Chien and Kim (2023). The diagnostic's interpretation and factor of two are
explained in [step 9](#9-calculate-the-growth-rate).

### Highest-priority suggestions

1. **Quantify reliability before interpreting the full map.** The denominator
   cutoff of `1e-12 * max(P_TT)` is a numerical safeguard. It can retain bins
   with very little temperature variance and consequently large ratios.
   Calculate heating power with the same transform, taper, averaging, and
   smoothing, then examine
   `coherence_squared = abs(S_QT)**2 / (S_QQ * S_TT)`. Combine coherence with
   a documented temperature-power criterion and uncertainty in the rate;
   coherence alone does not establish a precise ratio or causality. A useful
   approach is a block bootstrap of paired heating/temperature records,
   recomputing both averaged spectra and their ratio in each replicate.
   Resampling overlapping 96-day segments as independent observations would
   underestimate uncertainty. Account for temporal overlap and correlation
   between latitudes. The paper's coherence cutoff of 0.03 depends on its
   sampling and degrees of freedom; it should not be copied automatically.

2. **Distinguish the modal proxy from a mass-integrated EAPE budget.** The
   present SVD is orthonormal in an unweighted level-index metric. That is a
   defined modal diagnostic, but unequal pressure spacing gives unequal
   representation per unit column mass. Inspect each EOF's vertical shape
   before calling EOF1/EOF2 first/second baroclinic modes. For an actual column
   energy growth rate, formulate a consistent energy metric for both
   numerator and denominator, including pressure mass and reference
   stratification/temperature weights. Simply replacing the denominator with
   `cov_TT_avg` would be dimensionally inconsistent. Weighting an EOF analysis
   is a different scientific choice and should be tested as a sensitivity,
   not treated as a mandatory correction to the present modal formula.

3. **Test sensitivity to segment preprocessing.** Removing the full-record
   time mean is reasonable for this nonseasonal aquaplanet, but individual
   segments may retain a mean or drift. Compare the current result with
   segment demeaning, and optionally linear detrending, before tapering.
   Examine whether low-frequency features and their growth-rate signs
   persist. Such preprocessing also removes some long-period variability, so
   record and justify the choice. Compare segment lengths and edge-taper
   widths relevant to the waves of interest. A 96-day segment has frequency
   spacing `1/96 = 0.01042` cycles/day; smoothing reduces independent
   resolution further. Features near the longest resolvable periods merit
   particular caution.

### Additional checks worth making

- **EOF stability across experiments.** These EOFs are estimated separately
  for each `qflux`/`rh` combination. Compare their shapes, explained variance,
  and separation from neighboring eigenvalues. An EOF number may not identify
  the same structure in every experiment. Consider a common reference basis
  if direct modal comparisons are central to the analysis. Simultaneous EOF
  sign reversal of both heating and temperature PCs leaves this ratio
  unchanged; the arbitrary SVD sign is not itself a problem.

- **Validity of input arrays.** NetCDF variables have explicit fill values.
  `np.asarray(masked_array)` discards the mask and can expose finite fill
  values, while `np.nanmean` only handles NaNs. Standard SVD and FFT do not
  support missing samples directly. Check masks and finite values before
  decomposition and choose an explicit missing-data policy if necessary.
  This is a robustness concern, not evidence that the present input contains
  missing values.

- **Heating completeness and provenance.** The inspected file reports K/s
  for all four tendencies, supporting the conversion to K/day. Temperature
  and those tendencies all declare `time: mean`, so there is no obvious
  instantaneous-versus-mean mismatch in these metadata. Confirm from the
  exact Isca configuration that their sum contains each intended diabatic
  process once, and that `tdt_rad` already includes solar heating. Variable
  names and units alone cannot establish that bookkeeping.

- **Daily averaging and aliasing.** The data are daily means. Their spectra
  represent daily-mean fields and may differ from spectra of subdaily fields,
  especially near the 0.5 cycles/day Nyquist limit. Compare with higher-time-
  resolution output before assigning physical significance to very rapid
  damping or generation. A rate larger than 0.5 day⁻¹ is not by itself an
  FFT-frequency error: growth rate and oscillation frequency are different
  quantities.

- **Smoothing boundaries.** Zero padding attenuates the outermost frequency
  and wavenumber bins. Common normalization cancels for proportional input
  fields, but boundary treatment can matter for general spectra. Examine
  treatment of the periodic FFT boundary and the Nyquist bin explicitly.
  Also, the current filter mixes bins across zero frequency and wavenumber;
  the lowest-frequency and smallest-wavenumber bins are therefore not
  isolated propagation-direction estimates.

- **Column-power integration.** The increasing-pressure `dp/g` integral of
  full-field temperature power is a column variance diagnostic in
  K²·kg/m² per discrete spectral bin. It is not an EAPE energy spectrum.
  The retained native full levels span approximately 116–973 hPa in the
  default file. If an exact 100–1000 hPa layer is required, define boundary
  interpolation or use appropriate layer-interface mass weights. The
  coordinate is labeled approximate full pressure; if pressure variations
  matter, assess using actual layer pressures/masses. These choices affect
  the column-power context panel, not the present modal growth-rate ratio.

- **Area weighting and symmetry.** The arithmetic latitude mean is a small
  approximation within 5°S–5°N; Gaussian-grid area weights would make its
  meaning explicit and matter more over a wider band. Averaging local
  latitude spectra does not isolate equatorially symmetric waves. If
  identifying Kelvin waves is the goal, consider a symmetric/antisymmetric
  decomposition and a justified wave-band selection. These would be
  additional analyses, not purely cosmetic overlays.

## Historical validation and observations

The earlier growth-rate review recorded the following checks and results for
its then-current revision and default experiment. Its statement that only the
final plotting cell and documentation changed refers to that earlier update.
These observations describe that earlier revision and are retained as
historical context, rather than a new validation of the current scripts.

- The calculation through assignment of `sigma_1, sigma_2` was compared
  byte-for-byte with a pre-edit copy and is unchanged.
- The default experiment contains 3600 retained daily samples, 28 retained
  pressure levels, four tropical latitudes, and 128 uniformly spaced global
  longitudes. Its segmentation produces 98 overlapping segments.
- The two heating EOFs explain approximately 84.2% and 5.9% of unweighted
  heating variance in that experiment. These percentages do not describe
  the fraction of temperature variance or EAPE explained.
- In the displayed positive-frequency, `|k| <= 15` domain, mode 1 ranges
  approximately from −39.81 to +2.19 day⁻¹, and mode 2 from −11.68 to
  −1.00 day⁻¹. These are raw estimates without statistical screening; they
  should not be read as established wave decay times.
- Analytic checks of the existing FFT/smoothing helpers verify Parseval
  normalization, the eastward sign convention, the `Q = a T` rate `2a`, and
  zero co-spectral contribution for quadrature heating.
- At the time of the review, the nine-experiment batch runner completed
  successfully in `atmos-env`. All nine figures were regenerated as PNG/PDF
  pairs, with representative figures visually inspected and PDF text
  extraction checked. The current runners default to 18 experiments, and
  the current scripts save numerical diagnostics to NPZ files for separate
  visualization.
