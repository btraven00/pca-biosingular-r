# pca-biocsingular-r

`BiocSingular`-backed PCA module for omnibenchmark scRNA pipelines.

`BiocSingular::runSVD` dispatches on a `BSPARAM` object, so a single code path covers several
SVD backends and the **solver becomes a benchmark parameter**:

| `--solver` | `BSPARAM` | algorithm |
|---|---|---|
| `exact` | `ExactParam()` | LAPACK full SVD, truncated to `k` |
| `random` | `RandomParam()` | Halko randomized SVD, via `rsvd` (`p=10`, `q=2`) |
| `irlba` | `IrlbaParam()` | implicitly-restarted Lanczos, via `irlba` (`work = k + 7`) |

## Setup

```sh
pixi install
pixi run check
```

`pixi run check` loads all runtime libraries and prints `OK`. Run it after install to confirm the
environment is healthy.

## Conda environment export

```sh
pixi run export-env
```

Exports the resolved environment to `envs/pca-biocsingular-r.yml`. The environment is named after
the repo root folder. The file is generated — change `pixi.toml` and re-export rather than editing
it by hand.

## Usage

```sh
pixi run Rscript pca.R \
  --output_dir <dir> \
  --name <name> \
  --normalized_selected_h5 <normalized_selected.h5> \
  --solver {exact|random|irlba} \
  --n_components <int> \
  --random_seed <int> \
  [--dense {true|false}] \
  [--deferred {true|false}] \
  [--blas_threads <int>]
```

Outputs, both written to `<output_dir>`:

- `<name>_pcas.tsv` — header `cell_id<TAB>PC1<TAB>...<TAB>PC{n_components}`, one row per cell
- `<name>_loadings.tsv` — header `gene<TAB>PC1<TAB>...`, one row per gene

## Axes

**`--dense`** materialises the matrix before PCA. Only the dense path puts real work through
BLAS-3, so it is the axis that interacts with `--blas_threads`.

**`--blas_threads`** calls `RhpcBLASctl::blas_set_num_threads` at runtime, which overrides the
`OMP_NUM_THREADS` that Snakemake exports (it defaults to 1, which would otherwise pin every job
single-threaded regardless of the linked BLAS). `0` leaves it alone. Note that changing the thread
count perturbs floating-point reduction order: results are stable *within* a thread count and
differ at ~1e-15 *across* thread counts, so compare with `all.equal`, not `identical`.

**`--deferred`** applies centering as an operator over the original matrix instead of
materialising the centered matrix. The two are mathematically identical (measured: `1.8e-14`
relative difference for `random`, exactly `0` for `exact`), so this isolates an implementation
choice from a method choice.

Measured effect on peak memory (BiocSingular 1.28.0, serial, `fold=Inf`):

| solver | `--deferred false` | `--deferred true` |
|---|---|---|
| `random` | 660 MB | **449 MB** |
| `exact` | 652 MB | 642 MB |
| `irlba` | 358 MB | 359 MB |

Only `random` is affected. `runExactSVD` calls `safe_svd(as.matrix(x))` and materialises either
way; `runIrlbaSVD`'s serial branch never reads `deferred`, because irlba's native `center=`
already defers. The inert cells are kept in the grid deliberately — that deferral is *not* a
per-method property is the claim under test, and those cells are the evidence.

**`--random_seed`** is consumed by `--solver random` only. `exact` is deterministic, and irlba
converges to `tol` regardless of its random start vector. BiocSingular itself contains no RNG;
randomness enters only through the delegated `rsvd` and `irlba` calls, off R's global stream.

## Citation

If you use this module in your research, please cite it using the information in `CITATION.cff`.
Please also cite `BiocSingular` (Lun), and the underlying solver you used: Baglama & Reichel
(2005) for `irlba`, or Halko, Martinsson & Tropp (2011) for `random`.
