# Usage

## Choosing the right estimator

| Your data | Recommended function | Notes |
|---|---|---|
| Irregularly sampled, long (≥ 10 000 pts) | `compute_acf_gaussian_nufft` / `compute_acf_rectangle_nufft` | O(n log n), fastest |
| Irregularly sampled, any length | `compute_acf_gaussian_realspace` / `compute_acf_rectangle_realspace` | O(n) per lag, artifact-free |
| **Regularly sampled**, any method | `compute_acf_gaussian_fft` / `compute_acf_rectangle_fft` / `compute_acf_regular_fft` | 100–400× faster than NUFFT on regular grids |

All functions share the same signature:

```python
c, b = fn(lags, t, x)                     # no-kernel ("regular") variant
c, b = fn(lags, t, x, bin_width=0.5)      # gaussian / rectangle variants
```

- `lags` — array of lag values (same units as `t`)
- `t` — sample times, sorted ascending (float array, e.g. days since start)
- `x` — signal values, same length as `t`
- `bin_width` — kernel half-width (gaussian σ or rectangle half-width), same units as `t`
- Returns `(c, b)` — ACF estimate and effective pair count, both shape `(len(lags),)`

## Irregularly-sampled series (NUFFT path)

```python
import numpy as np
import pandas as pd
from nufftacf import compute_acf_gaussian_nufft, t_numeric_of

# A pandas Series with a DatetimeIndex works directly
idx = pd.date_range("2005-01-01", periods=3650, freq="D")
# simulate some irregular gaps
mask = np.random.default_rng(0).random(3650) > 0.3
sts = pd.Series(np.random.randn(3650)[mask], index=idx[mask])

lags = np.arange(1.0, 366.0)          # lags 1..365 days
t    = t_numeric_of(sts)               # timedelta -> float days
x    = sts.to_numpy()

c, b = compute_acf_gaussian_nufft(lags, t, x, bin_width=0.5)
```

## Regularly-sampled series (classic FFT path)

```python
from nufftacf import compute_acf_gaussian_fft, compute_acf_regular_fft

n = 3650
t = np.arange(n, dtype=float)         # integer day indices
x = np.random.default_rng(0).standard_normal(n)
lags = np.arange(0.0, 366.0)

# With Gaussian smoothing kernel
c_gauss, b = compute_acf_gaussian_fft(lags, t, x, bin_width=0.5)

# Without smoothing kernel (matches Pastas' bin_method="regular" exactly)
c_reg, b = compute_acf_regular_fft(lags, t, x)
```

!!! warning
    The `_fft` estimators raise `ValueError` if `t` is not regularly spaced.
    Use the `_nufft` or `_realspace` estimators for irregular data.

## Notebooks

Two Colab-ready notebooks are included in `notebook/`:

- [`pastas_vs_nufftact.ipynb`](https://github.com/jecampagne/nufftacf/blob/main/notebook/pastas_vs_nufftact.ipynb) —
  **irregularly-sampled** series (sine and AR(1)-like, with gaps), NUFFT path vs Pastas.
- [`pastas_vs_nufftacf_regular.ipynb`](https://github.com/jecampagne/nufftacf/blob/main/notebook/pastas_vs_nufftacf_regular.ipynb) —
  **regularly-sampled** series (sine, noisy sine, exponential, constant, square wave),
  classic-FFT path vs Pastas for all three bin methods.

Both install `nufftacf` directly from GitHub in their first cell — no local setup needed
to run them on Colab.

## Benchmark scripts

```bash
pip install -e ".[benchmark]"

# Irregular data: Pastas vs _nufft
python benchmark/benchmark_acf.py
python benchmark/fit_benchmark_acf.py benchmark_acf_results.csv

# Regular data: Pastas vs _fft vs _nufft, all three bin methods
python benchmark/benchmark_acf_regular.py
python benchmark/fit_benchmark_acf_regular.py --csv_path benchmark_acf_regular_results.csv
```

Example results (Apple Silicon Mac) are stored in `benchmark/`.
