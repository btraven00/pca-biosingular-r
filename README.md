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
  [--blas_threads <int>] \
  [--oversampling <int>] \
  [--power_iters <int>]
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

**`--oversampling` (`p`) and `--power_iters` (`q`)** size the randomized sketch, `--solver random`
only; they are forwarded to `rsvd::rsvd` through `RandomParam(...)`. Defaults restate rsvd's own
(`p=10`, `q=2`) and are echoed in the run log rather than left implicit. Passing either to `exact`
or `irlba` warns on stderr and is ignored — `ExactParam` has no sketch, and `IrlbaParam`'s analogue
is `extra.work`, not `p`/`q`.

`q` is the parameter that decides whether trailing PCs are resolved at all, which makes it the
axis rather than a default worth inheriting. **scikit-learn's `randomized_svd` picks `n_iter=7` at
these shapes** (`n_iter="auto"` → `7` when `n_components < 0.1 * min(dim)`; verified identical in
sklearn 1.2.0 and 1.8.0), i.e. 3.5× more than rsvd. Comparing the two libraries on their
respective defaults compares the defaults, not the libraries.

Measured on be1-fixture (2000 × 1715, `k=50`), mean per-PC `|cor|` against `ExactParam`:

| `random` config | PC1-20 | PC21-30 | PC31-40 | PC41-50 |
|---|---|---|---|---|
| `q=2, p=10` *(default)* | 0.98 | 0.619 | 0.275 | 0.115 |
| `q=7, p=10` | 1.00 | 0.994 | 0.866 | 0.222 |
| `q=7, p=20` | 1.00 | 0.998 | 0.954 | 0.365 |
| *scanpy randomized, for reference* | 1.00 | 0.996 | 0.710 | 0.302 |

At `q=7` this module lands on scanpy's randomized arm; the whole gap was the default. The decay
past PC30 that survives at every setting is the flat tail of the spectrum — consecutive singular
values differ by <1% there, so those components are not individually identifiable and no solver
setting recovers them.

Seeds are not a substitute for `q`: the error is bias, not zero-mean noise. Over 10 seeds at
`q=2` the median is 0.579 in PC21-30 — no better than a single run — while the spread runs 0.013
to 0.942. Raise `q`, do not average seeds.

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
