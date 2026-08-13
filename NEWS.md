# mitodrift 1.0.2

- Initial public release.
- Add opt-in Metropolis-coupled MCMC (MC3) with persistent heated states,
  adjacent-temperature swaps, checkpoint/resume support, and swap diagnostics.
- Make native linking portable across Linux and macOS by using R's standard
  BLAS/LAPACK/Fortran linker variables.
- Qualify data.table helpers so the package works without attaching data.table.
