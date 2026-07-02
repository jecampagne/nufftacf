# nufftacf

[![Tests](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml/badge.svg)](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml)
[![Lint](https://github.com/jecampagne/nufftacf/actions/workflows/lint.yml/badge.svg)](https://github.com/jecampagne/nufftacf/actions/workflows/lint.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/jecampagne/nufftacf/blob/main/LICENSE)

Fast autocorrelation function (ACF) estimation for **irregularly- and regularly-sampled** 1D
time series, scaling as O(n log n).

## Three estimator families

| Function(s) | Sampling | Method | Scaling |
|---|---|---|---|
| `compute_acf_gaussian_nufft` / `_rectangle_nufft` | irregular | NUFFT + Wiener-Khinchin | ~O(n log n) |
| `compute_acf_gaussian_realspace` / `_rectangle_realspace` | irregular or regular | direct weighted sum | O(n) per lag |
| `compute_acf_gaussian_fft` / `_rectangle_fft` / `_regular_fft` | **regular only** | classic FFT + filter | ~O(n log n) |

All functions share the same calling convention: `fn(lags, t, x, bin_width=0.5)`, returning
`(c, b)` — the ACF estimate and effective pair count.

## Quick start

```python
import numpy as np
from nufftacf import compute_acf_gaussian_nufft, t_numeric_of
import pandas as pd

sts = pd.Series(...)  # irregularly-sampled pandas Series with DatetimeIndex
lags = np.arange(1.0, 366.0)
c, b = compute_acf_gaussian_nufft(lags, t_numeric_of(sts), sts.to_numpy(), bin_width=0.5)
```

```python
# Regularly-sampled series: use the faster classic-FFT path
from nufftacf import compute_acf_gaussian_fft
t = np.arange(len(x), dtype=float)
c, b = compute_acf_gaussian_fft(lags, t, x, bin_width=0.5)
```

## Method

Built on two ingredients:

1. **[FINUFFT](https://github.com/flatironinstitute/finufft)** (Flatiron Institute)
   evaluates the power spectrum of the irregularly-sampled signal via a type-1 NUFFT,
   then inverts it at the requested lags (Wiener-Khinchin theorem) — giving ~O(n log n)
   scaling instead of the O(n²) real-space sum.
2. An analytical **"b" correction** for the effective pair count per lag, specific to
   the Gaussian and rectangular kernels, computed in O(n) via a two-pointer scan
   (`kernels.py`) — without this, the raw NUFFT spectrum would not be correctly normalised.

On a **regular** grid, `fft_acf.py` uses plain `scipy.signal.correlate` and recovers
the same `b` denominator by smoothing the triangular "raw pair-count" ramp with the
same discrete filter as the numerator — making any discretisation artefact cancel exactly
in the ratio.

## Notebooks

- [pastas_vs_nufftact.ipynb](https://github.com/jecampagne/nufftacf/blob/main/notebook/pastas_vs_nufftact.ipynb) —
  **irregularly-sampled** series, NUFFT path vs Pastas (Colab-ready).
- [pastas_vs_nufftacf_regular.ipynb](https://github.com/jecampagne/nufftacf/blob/main/notebook/pastas_vs_nufftacf_regular.ipynb) —
  **regularly-sampled** series, classic-FFT path vs Pastas (Colab-ready).

## Benchmark

`benchmark/` contains timing results on an Apple Silicon Mac (irregular case)
and a Linux sandbox (regular case), along with the scripts to reproduce them.
See [Usage](usage.md) for the commands.

## Citing

If you use `acf-nufft`, please also cite FINUFFT:

> A. H. Barnett, J. F. Magland, and L. af Klinteberg (2019).
> *A parallel non-uniform fast Fourier transform library based on an
> "exponential of semicircle" kernel.* SIAM J. Sci. Comput.
> [github.com/flatironinstitute/finufft](https://github.com/flatironinstitute/finufft)
