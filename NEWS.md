# mitodrift 1.0.2

- Initial public release.
- Add opt-in Metropolis-coupled MCMC (MC3) with persistent heated states,
  adjacent-temperature swaps, checkpoint/resume support, and swap diagnostics.
- Make native linking portable across Linux and macOS by using R's standard
  BLAS/LAPACK/Fortran linker variables.
- Qualify data.table helpers so the package works without attaching data.table.
- Accumulate ASDSF clade counts incrementally across MCMC batches instead of
  rescanning the full trace at every convergence check.
- Add configurable MCMC persistence cadence through `checkpoint_every` and
  `--tree_mcmc_checkpoint_every`; final traces are always persisted.
- Cache internal-edge ordinal lookup in the NNI proposal engine.
- Enforce `max_iter` as a safety cap when using convergence-threshold mode.
