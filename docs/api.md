# API Reference

All public functions are importable directly from `nufftacf`:

```python
from nufftacf import compute_acf_gaussian_nufft, compute_acf_regular_fft, ...
```

---

## NUFFT estimators (irregular and regular data)

::: nufftacf.nufft_acf.compute_acf_gaussian_nufft

::: nufftacf.nufft_acf.compute_acf_rectangle_nufft

---

## Real-space estimators (irregular and regular data)

::: nufftacf.realspace_acf.compute_acf_gaussian_realspace

::: nufftacf.realspace_acf.compute_acf_rectangle_realspace

---

## Classic-FFT estimators (regular data only)

::: nufftacf.fft_acf.compute_acf_regular_fft

::: nufftacf.fft_acf.compute_acf_rectangle_fft

::: nufftacf.fft_acf.compute_acf_gaussian_fft

---

## Kernel helpers

::: nufftacf.kernels.compute_b_gaussian

::: nufftacf.kernels.compute_b_rectangle

::: nufftacf.kernels.compute_c_gaussian

::: nufftacf.kernels.compute_c_rectangle

---

## Utilities

::: nufftacf.utils.t_numeric_of

::: nufftacf.utils.standardize
