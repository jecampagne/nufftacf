# nufftacf ➜ renamed to [nufftcf](https://github.com/jecampagne/nufftcf)

**This repository has been renamed and is no longer maintained here.**

`nufftacf` originally covered only autocorrelation (ACF) estimation for
irregularly-sampled 1D signals via the Non-Uniform FFT (hence the "a" in
the name). It has since been extended to also compute the
**cross-correlation function (CCF)** between two series, so the project
was renamed to reflect that broader scope.

## ➡️ Please go to [github.com/jecampagne/nufftcf](https://github.com/jecampagne/nufftcf)

- Full commit history has been preserved in the new repository.
- Latest code, documentation, and notebooks (ACF **and** CCF) now live
  there: <https://jecampagne.github.io/nufftcf/>
- The package is published on PyPI:
  ```bash
  pip install nufftcf
  ```
  (the `nufftacf` name on PyPI, if you happened to depend on it, is
  deprecated in favor of `nufftcf`)

If you have an existing local clone of this repository, update your
remote instead of re-cloning:

```bash
git remote set-url origin https://github.com/jecampagne/nufftcf.git
git fetch origin
git checkout main
```

and update any `import nufftacf` in your code to `import nufftcf`.

## License

[MIT](https://github.com/jecampagne/nufftacf/blob/main/LICENSE) -- unchanged in the new repository.
