# API Reference

All public functions are importable directly from `nufftacf`:

```python
from nufftacf import compute_acf_gaussian_nufft, compute_acf_regular_fft, ...
```

---

## NUFFT estimators

::: nufftacf.nufft_acf
    options:
      show_root_heading: true
      show_signature: true
      show_signature_annotations: true
      separate_signature: true
      members:
        - compute_acf_gaussian_nufft
        - compute_acf_rectangle_nufft


## Classic-FFT estimators (regular data only)

::: nufftacf.fft_acf
    options:
      show_root_heading: true
      show_signature: true
      show_signature_annotations: true
      separate_signature: true
      members:
        - compute_acf_regular_fft
        - compute_acf_rectangle_fft
        - compute_acf_gaussian_fft

## Real-space estimators

::: nufftacf.realspace_acf
    options:
      show_root_heading: true
      show_signature: true
      show_signature_annotations: true
      separate_signature: true
      members:
        - compute_acf_gaussian_realspace
        - compute_acf_rectangle_realspace
---

## Kernel helpers
::: nufftacf.kernels
    options:
      show_root_heading: true
      show_signature: true
      show_signature_annotations: true
      separate_signature: true
      members:
        - compute_b_gaussian
        - compute_b_rectangle
        - compute_c_gaussian
        - compute_c_rectangle

---

## Utilities
::: nufftacf.utils
    options:
      show_root_heading: true
      show_signature: true
      show_signature_annotations: true
      separate_signature: true
      members:
        - t_numeric_of
        - standardize
