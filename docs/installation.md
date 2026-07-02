# Installation

## Requirements

Python >= 3.11. Core dependencies (installed automatically): `numpy`, `pandas`,
`numba`, `scipy`, `finufft`.

## Recommended: install from a local clone

```bash
git clone https://github.com/jecampagne/nufftacf.git
cd nufftacf
python -m venv venv
source venv/bin/activate          # Windows: venv\Scripts\activate
pip install --upgrade pip setuptools wheel

# Force prebuilt wheels for the three compiled dependencies.
# This avoids source builds that can fail or produce mismatched OpenMP runtimes,
# particularly on macOS (see Troubleshooting below).
pip install --only-binary=:all: finufft numba llvmlite

pip install -e .
# with optional benchmark dependencies (Pastas + matplotlib):
pip install -e ".[benchmark,test]"
```

## Install directly from GitHub

```bash
pip install "nufftacf @ git+https://github.com/jecampagne/nufftacf.git"
# with benchmark extras:
pip install "nufftacf[benchmark] @ git+https://github.com/jecampagne/nufftacf.git"
```

!!! note
    For maximum reliability, run `pip install --only-binary=:all: finufft numba llvmlite`
    before the above when installing in a fresh environment on macOS.

## Troubleshooting (macOS)

If you hit a segfault or hang when running tests or importing the package,
see the **macOS troubleshooting** section in the [README](index.md).
The short version: make sure you have a clean Python environment (no conda
stacked on top of your venv), install `finufft`/`numba`/`llvmlite` as
prebuilt wheels (`--only-binary=:all:`), and clear any stale Numba
on-disk JIT cache:

```bash
rm -rf ~/.numba_cache
find . -type d -name "__pycache__" -exec rm -rf {} +
pip install -e ".[test]"
pytest tests/
```
