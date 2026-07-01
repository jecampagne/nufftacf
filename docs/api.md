# API Reference

All public functions are importable directly from `acf_nufft`:

```python
from acf_nufft import compute_acf_gaussian_nufft, compute_acf_regular_fft, ...
```

---

## NUFFT estimators (irregular and regular data)

::: acf_nufft.nufft_acf.compute_acf_gaussian_nufft

::: acf_nufft.nufft_acf.compute_acf_rectangle_nufft

---

## Real-space estimators (irregular and regular data)

::: acf_nufft.realspace_acf.compute_acf_gaussian_realspace

::: acf_nufft.realspace_acf.compute_acf_rectangle_realspace

---

## Classic-FFT estimators (regular data only)

::: acf_nufft.fft_acf.compute_acf_regular_fft

::: acf_nufft.fft_acf.compute_acf_rectangle_fft

::: acf_nufft.fft_acf.compute_acf_gaussian_fft

---

## Kernel helpers

::: acf_nufft.kernels.compute_b_gaussian

::: acf_nufft.kernels.compute_b_rectangle

::: acf_nufft.kernels.compute_c_gaussian

::: acf_nufft.kernels.compute_c_rectangle

---

## Utilities

::: acf_nufft.utils.t_numeric_of

::: acf_nufft.utils.standardize
