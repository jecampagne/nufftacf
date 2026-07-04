# nufftacf

[![Tests](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml/badge.svg)](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml)
[![Lint](https://github.com/jecampagne/nufftacf/actions/workflows/lint.yml/badge.svg)](https://github.com/jecampagne/nufftacf/actions/workflows/lint.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)


[![Python](https://img.shields.io/badge/python-3.11%20|%203.12%20|%203.13%20|%203.14-blue)](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml)
[![Linux](https://img.shields.io/badge/linux-ubuntu--latest-orange?logo=linux)](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml)
[![macOS](https://img.shields.io/badge/macOS-latest%20(arm64)%20%7C%20intel-orange?logo=apple)](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml)
[![macOS](https://img.shields.io/badge/macOS-arm64%20%7C%20x86__64-orange?logo=apple)](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml)
[![Windows](https://img.shields.io/badge/windows-latest-orange?logo=windows)](https://github.com/jecampagne/nufftacf/actions/workflows/tests.yml)



Fast autocorrelation function (ACF) estimation for **irregularly- and
regularly-sampled** time series, scaling as $O(n\log n)$ thanks noytably to the Nonuniform Fast Fourier Transform library developped by the Flatiron Institut ([FINUFFT](https://github.com/flatironinstitute/finufft)). 

With **`nufftacf`** three estimator families are provided:

| Function | Sampling | Method | Scaling | Notes |
|---|---|---|---|---|
| `compute_acf_gaussian_nufft` | irregular | NUFFT + Wiener-Khinchin | $\sim~O(n\log n)$ | fastest for long irregular series; ~1-3% residual amplitude bias on strongly periodic signals (see below) |
| `compute_acf_rectangle_nufft` | irregular | NUFFT + Wiener-Khinchin | $\sim~O(n\log n)$ | same caveat as above |
| `compute_acf_gaussian_realspace` | irregular or regular | direct real-space weighted sum | $O(n)$ per lag | artifact-free reference |
| `compute_acf_rectangle_realspace` | irregular or regular | direct real-space weighted sum | $O(n)$ per lag | artifact-free reference |
| `compute_acf_regular_fft` | **regular only** | classic FFT correlation, no kernel | $\sim~O(n)$ | matches Pastas `bin_method="regular"` to numerical precision |
| `compute_acf_rectangle_fft` | **regular only** | classic FFT correlation + box filter | $\sim~O(n\log n)$ | faster than `_nufft`/`_realspace` on regular data (no NUFFT/numba overhead) |
| `compute_acf_gaussian_fft` | **regular only** | classic FFT correlation + gaussian filter |$\sim~O(n\log n)$ | same |

All seven share the same calling convention: `fn(lags, t, x, bin_width=0.5)`
(`compute_acf_regular_fft` has no `bin_width`, since it applies no
smoothing kernel), and return `(c, b)` -- the ACF estimate and the
effective pair count, both shape `(len(lags),)`.

## Documentation

https://jecampagne.github.io/nufftacf/

## Installation

```bash
# from a local clone
git clone https://github.com/jecampagne/nufftacf.git
cd nufftacf
pip install -e .
# or with the optional benchmark (Pastas + matplotlib) :
pip install -e ".[benchmark]"

# from GitHub without any local clone
pip install "nufftacf @ git+https://github.com/jecampagne/nufftacf.git"
```

Requires Python >= 3.11. Core dependencies: numpy, pandas, numba, scipy, finufft.

## Quick start

```python
import numpy as np
import pandas as pd
from nufftacf import compute_acf_gaussian_nufft, t_numeric_of

# An irregularly-sampled series (any DatetimeIndex works)
idx = pd.date_range("2000-01-01", periods=5000, freq="D")[np.random.rand(5000) > 0.2]
x = pd.Series(np.random.randn(len(idx)), index=idx)

lags = np.arange(1.0, 366.0)         # 1 to 365 days
t = t_numeric_of(x)                   # elapsed days since first sample

c, b = compute_acf_gaussian_nufft(lags, t, x.to_numpy(), bin_width=0.5)
# c: ACF estimate per lag (c ~ 1 at lag -> 0)
# b: effective number of contributing pairs per lag (useful to flag
#    under-sampled lags, e.g. mask out lags where b is too small)
```



## Which estimator should I use?

- **Your data is regularly sampled** (a fixed time step, no gaps): use the
  `_fft` variants. They're faster than both `_nufft` (no NUFFT overhead) and
  `_realspace` (no numba two-pointer scan) on regular data, and
  `compute_acf_regular_fft` additionally gives you Pastas' "regular"
  bin_method (no smoothing kernel) at a fraction of its cost (see
  `benchmark/`).
- **Long, irregularly-sampled series** (tens of thousands of points or more)
  where Pastas' real-space approach becomes impractically slow: use the
  `_nufft` variants.
- **Strongly periodic, irregularly-sampled signals** (e.g. seasonal/annual
  cycles) where you need the most accurate possible ACF and series length is
  manageable: use the `_realspace` variants, or the `_nufft` variants with an
  increased `N1` (e.g. `N1>32*len(x)`), which reduces but does not fully
  eliminate the residual bias (see below).
- **Everything else, irregular case**: either `_nufft` or `_realspace` works;
  `_nufft` will generally be faster.

### A note on the NUFFT residual bias

The NUFFT-based estimators compute the power spectrum of the irregularly-sampled signal and invert it at the requested lags via the Wiener-Khinchin theorem. This implicitly relies on a finite-domain Fourier representation, which is mathematically equivalent to convolving the true spectrum with the "spectral window" induced by the irregular/gappy sampling pattern. A narrow spectral peak (a strongly periodic signal) is distorted much more visibly by this convolution than a broad, featureless spectrum (e.g. an AR(1)-type exponential decay), even though the absolute size of the distortion is similar in both cases.

In practice, with the default `N1 = 32 * len(x)` (the number of Fourier modes used internally by FINUFFT), this residual bias is on the order of 1–3% of the ACF amplitude for strongly periodic signals with irregular or gappy sampling, and negligible for smoothly-decaying, broadband signals. Reducing N1 speeds up the computation slightly at the cost of a larger bias; increasing it beyond `32 * len(x)` gives diminishing returns for most practical series.

The `_realspace` estimators do not have this limitation (no implicit
periodicity assumption), at the cost of O(n) scaling per lag rather than
O(n log n) -- for most practical series lengths both are fast; benchmark
on your own data if it matters (see `benchmark/`).

### Regularly-sampled data: the `_fft` estimators

When `t` is on a regular grid, `compute_acf_regular_fft` /
`compute_acf_rectangle_fft` / `compute_acf_gaussian_fft` (in `fft_acf.py`)
skip NUFFT entirely and use a plain `scipy.signal.correlate` (classic FFT
correlation) instead -- faster, and with no finufft/numba dependency in the
hot path. All three raise `ValueError` if `t` isn't regularly spaced (use
`_nufft`/`_realspace` for that).

- `compute_acf_regular_fft` reproduces Pastas' `bin_method="regular"`
  (a windowed Pearson correlation, no smoothing kernel) to numerical
  precision (`atol=1e-9` in `tests/test_fft_acf.py`), via an O(n) cumulative
  -moments computation (`E[X^2] - E[X]^2`) instead of one `np.corrcoef` call
  per lag.
- `compute_acf_rectangle_fft` / `compute_acf_gaussian_fft` match
  `compute_acf_rectangle_realspace` / `compute_acf_gaussian_realspace`
  almost exactly at the package's default `bin_width=0.5` (`atol=1e-9`).
  For other `bin_width` values, expect a small residual at `lag=0`
  specifically (a few %, decaying to <0.1% by lag~5) -- an inherent
  discretization artifact of approximating a continuous symmetric kernel
  window with a discrete digital filter, not a bug to chase further; see
  `tests/test_fft_acf.py::test_rectangle_fft_matches_realspace_various_bin_widths`
  for the exact numbers across `bin_width` values.

## Notebooks

- [`pastas_vs_nufftact.ipynb`](notebook/pastas_vs_nufftact.ipynb)
  compares **nufftacf** against **Pastas** on **irregularly**-sampled series
  (sine and AR(1)-like, with random gaps), using the `_nufft` estimators.
- [`pastas_vs_nufftacf_regular.ipynb`](notebook/pastas_vs_nufftacf_regular.ipynb)
  does the same on **regularly**-sampled series (sine, noisy sine,
  noisy exponential decay, square wave), using the `_fft` estimators,
  for all 3 of Pastas' bin methods (`regular`/`rectangle`/`gaussian`).
- [`zdcf_vs_nufftacf.ipynb`](`notebook/zdcf_vs_nufftacf.ipynb) compares **nufftacf** against **pyzdcf** on the same **irregularly**-sampled series used in the `pastas_vs_nufftact.ipynb`.

All are Colab-ready: the first cell installs **nufftacf** as well as **Pastas** or **pyzdcf** and third party libraries. Concerning **pyzdcf**, the repository was cloned and adapted to ensure compatibility with the pandas and other library versions used in this notebook, allowing it to run on Google Colab. These changes do not affect the quality of the computations.

## Method

`nufftacf` is built on two ingredients:

1. **[FINUFFT](https://github.com/flatironinstitute/finufft)** (Flatiron
   Institute) to evaluate the power spectrum of the irregularly-sampled
   signal via a type-1 non-uniform FFT, then invert it at the requested lags
   via a type-2 NUFFT (Wiener-Khinchin theorem) -- this is what gives the
   `_nufft` estimators their $\sim~O(n\ log\ n)$ scaling, instead of the $O(n^2)/O(n)$
   per-lag direct sum.
2. An analytical, kernel-specific correction for the number of
   contributing sample pairs per lag (the `b` denominator in `kernels.py`),
   for both the **Gaussian** and **rectangular/boxcar** smoothing kernels --
   computed with an $O(n)$ two-pointer scan (since `t` is sorted), rather than
   the naive $O(n^2)$ all-pairs count. This `b` is what turns the raw NUFFT
   power spectrum into a properly normalized correlation.

On a **regular** grid, `fft_acf.py` gets the same `b` correction for free,
without the two-pointer scan: smoothing the deterministic "raw pair count"
ramp (`n - lag`) with the *same* discrete filter (gaussian or box) used for
the correlation numerator reproduces `b` exactly -- see the "Regularly
-sampled data" section above for the validation numbers and the two bugs
this caught in the original prototype.

See `nufft_acf.py`, `kernels.py` and `fft_acf.py` docstrings for the full
derivation, and `notebook/` / `benchmark/` for empirical validation.

## Benchmark

- `benchmark/benchmark_acf.py` (+ `fit_benchmark_acf.py`): **Pastas** vs
  `_nufft`, on **irregularly**-sampled series of varying length, both
  kernels.
- `benchmark/benchmark_acf_regular.py` (+ `fit_benchmark_acf_regular.py`):
  Pastas vs `_fft` *and* `_nufft`, on **regularly**-sampled series of
  varying length, all 3 bin methods -- this is what lets you see, on
  regular data, how much the dedicated `_fft` path buys over just reusing
  the more general `_nufft` estimator.
- no real benchmarks are provided to compare **nuffacf** against **pyzdcf**, although in the plots obtained in `zdcf_vs_nufftacf.ipynb` one can appreciate that **nufftacf** is ~2 order of magnitude faster. 

```bash
pip install -e ".[benchmark]"
python benchmark/benchmark_acf.py            
# -> benchmark_acf_results.csv
python benchmark/fit_benchmark_acf.py

python benchmark/benchmark_acf_regular.py    
# -> benchmark_acf_regular_results.csv
python benchmark/fit_benchmark_acf_regular.py
```

Adjust `durations_years` / `n_points_list` and the **Pastas** cutoffs
(`pastas_max_years`, `pastas_max_n_regular`, `pastas_max_n_kernel` -- **Pastas**'
"gaussian"/"rectangle" bin methods are $O(n^2)$ on regular data too, just like
on irregular data, while "regular" is empirically $\sim O(n)$ and stays usable
much longer; both were measured directly before picking these defaults, not
assumed) at the top of each script as needed. Each measurement uses several
repeats and keeps the minimum, to reduce noise from shared/cloud
environments (Colab, background browser activity, etc.).

`benchmark/*_macosx.{csv,png}` give the results on MacBook Pro (2020) 2 GHz Intel Core i5 quatre cœurs (osx Tahoe 26.5.1)-- re-run on your own machine for comaparison. You can share your results on the `Discussions` section  of the repository.

## Citing

If you use `nufftacf`, please also cite FINUFFT, which it depends on:

> A. H. Barnett, J. F. Magland, and L. af Klinteberg (2019).
> *A parallel non-uniform fast Fourier transform library based on an
> "exponential of semicircle" kernel.* SIAM J. Sci. Comput.  41(5), C479-C504.
> https://github.com/flatironinstitute/finufft

> J.E Campagne (2026): *"Non Uniform FFT based Auto Correlation functions"*.  https://github.com/jecampagne/nufftacf

## License

[MIT](LICENSE)

## Development

```bash
pip install -e ".[dev]"
black .          # formatage
pytest tests/    # tests
```

## Tests

```bash
pip install -e ".[test]"
pytest tests/
```

`tests/test_nufftacf.py` (NUFFT vs realspace, irregular data) and
`tests/test_fft_acf.py` (fft vs realspace, and fft "regular" vs Pastas
itself, regular data) are correctness/sanity checks, not performance
benchmarks.
