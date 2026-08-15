# Bidirectional (depth-independent) NNI cache

## Goal

`NNICache::compute_new_loglik_impl` evaluates an NNI proposal by walking from
the modified edge's parent up to the root, an O(depth) operation on every
single MCMC/MC3 iteration (accepted or not). On imbalanced, ladder-like trees
this is expensive: on `MM_BB094_t1_500` the mean proposal walk-depth is 85.9
edges (vs. ~9 for a balanced 500-tip tree), and measured ordinary-MCMC cost
was 5.64 ms/iteration there vs. 0.69 ms/iteration on the near-balanced
`pL1000_200_2` (200 tips, mean depth 6.9) — an 8x gap despite only 2x more
variants. The goal was to make proposal evaluation depth-independent by
maintaining a second "outside" cache (`G`) alongside the existing "inside"
cache (`F`).

## Attempt 1: eager `G` maintenance — correct, but a net loss

Kept `G` eagerly consistent on every accept by rebuilding whatever subtree(s)
went stale immediately, using `recompute_G_subtree`. Validated correct
(`compute_new_loglik_fast` agreed with the O(depth) baseline to floating-point
noise across thousands of proposals, including repeated-accept sequences),
but **not a net win**: the necessary eager rebuild touches the *combined
size* of the sibling subtrees hanging off the accepted path, not the path's
*length*. On anything but a deep ladder shape those subtrees are large
(measured: `pL1000_200_2`'s root has child subtrees of 81 and 317 nodes out
of 399), so wiring it in made `pL1000_200_2` **18x slower** (51.4s vs 2.8s
for 4000 iterations). Reverted; the `G` cache and its supporting math were
kept in place as a validated foundation (see the message-passing math
cross-checked against `em_helpers.cpp`'s already-production inside-outside
implementation for the EM E-step).

## Attempt 2: lazy invalidation for both `G` and `F` — correct, and a real win

The eager rebuild's cost comes from doing work at *commit* time that only
*some future evaluation* might actually need. The fix is to defer it:
mark a node's cache stale in O(1) at commit time, and only pay to repair it
the first time some later query actually reaches it.

### Design

Both `G` and `F` get a two-layer lazy scheme:

- **Content version** (`{g,f}_content_version[v]`): bumped whenever a node's
  cached value is *actually recomputed*. A node's cache is trustworthy only
  if the version of whatever it was last computed from (its parent for `G`,
  both children for `F`) still matches the snapshot it recorded at that time.
- **Per-epoch verified cache** (`{g,f}_verified_epoch[v]`): an O(1) fast path
  on top of the version check. Once a node is confirmed valid within the
  current epoch (unchanged since the last commit anywhere in the cache), a
  later query for it — or for a descendant/ancestor whose walk reaches it —
  can skip straight past without re-deriving anything. `mark_g_dirty` /
  writing new `F` at commit time bump the epoch, invalidating this cache
  cache-wide so the next query is forced through the real check again.

`G[v]` depends on `(G[parent(v)], F[sibling(v)])`, so `ensure_g_valid`
recurses through `parent_of`. `F[v]` depends on `(F[child0(v)], F[child1(v)])`,
so `ensure_f_valid` recurses through `children_of` — the mirror-image
direction. This asymmetry matters for what needs explicit marking at commit
time: `G`'s dependency (parent *and sibling*) doesn't align with the
parent-of recursion the same way `F`'s does, so `G` still needs explicit
per-level "mark this sibling dirty" calls walking from `p1` to the root (see
`commit_fast`) — cheap (O(depth) marks, no recompute), but not free the way
`F`'s cascade is. `F` needs no such explicit ancestor marking at all: once
`p1`/`p2` get their new content version, any ancestor's content-version
check against them fails automatically the first time something queries that
ancestor, with no separate bookkeeping.

`commit_fast` writes `F[p1]`/`F[p2]` directly, reusing values
`compute_new_loglik_fast` already computed during evaluation (one extra
cheap combine for `F[p1]`, since evaluation alone doesn't need it) — zero
extra compute at commit time, just a copy. Everything above `p1` is left
exactly as-is, to be repaired lazily only if and when something later
actually reads it.

### Why this needed a different entry point than attempt 1

Reusing `compute_new_loglik_fast` + `commit_staged_nni(sync_outside_cache=true)`
as before was tried first and made things *worse* even with lazy `G`: every
accepted proposal (which, empirically, is 72-83% of proposals on these test
cases — see below) paid for the fast evaluation *and then* a full redundant
O(depth) re-evaluation via `compute_new_loglik(..., true)` to get the F
values `commit_staged_nni` needs, since that function's own O(depth) walk has
no way to reuse work `compute_new_loglik_fast` already did. The fix: a
dedicated `commit_fast(new_total_loglik)` that applies what evaluation
already computed directly, and a separate `fast_staged_*` state so this path
never touches `compute_new_loglik`/`commit_staged_nni`'s own `staged_*`
fields. The two commit paths coexist untouched by each other.

### Correctness: two more real bugs found before this validated clean

Both found by exact trajectory comparison against the O(depth) baseline
(`tree_mcmc_cpp_cached_threadsafe` vs. `tree_mcmc_cpp_cached_threadsafe_fast`,
same seed, same RNG draws, comparing the full accepted/rejected tree sequence
— not just isolated evaluations), not by inspection:

1. **Snapshot-before-topology-update ordering bug.** `commit_fast` recorded
   each node's child-content-version snapshot using `children_of[p1]`/
   `children_of[p2]` *before* the topology bookkeeping updated them to the
   new post-swap children — so the recorded snapshot didn't match what the
   node's children actually were going forward, breaking the very next
   evaluation that touched it. Fixed by moving the version/snapshot writes
   to after the topology update.
2. **Missing epoch bump on direct writes.** `commit_fast` writes `F[p1]`/
   `F[p2]` directly (not via `ensure_f_valid`'s own repair path) and set
   their `f_verified_epoch` to the *current* `f_epoch` — but never advanced
   `f_epoch` itself. Every previously-verified ancestor kept short-circuiting
   on the O(1) `f_verified_epoch[v] == f_epoch` fast path forever after,
   because nothing ever made that comparison fail again — so the *correct*
   content-version check one line below it was never reached. `ensure_f_valid`
   itself was right the whole time; nothing had ever told it a query was
   overdue. Fixed with `++f_epoch` in `commit_fast`, mirroring what
   `mark_g_dirty` already does for `G`. Found via a dedicated debug tool
   (`debug_verify_f_consistency`, kept in the tree as
   `nni_cache_replay_and_verify_f_cpp`) that recomputes `F` from scratch,
   bottom-up, and diffs it against the lazily-maintained cache after a given
   replayed proposal sequence — isolated the bug to "wrong after commit 3"
   rather than needing to guess from a diverging downstream loglik.

After both fixes: exact trajectory identity confirmed across `small_test`,
`pL1000_200_2`, and `MM_BB094_t1_500`, multiple seeds each, up to 15,000
iterations, plus the original isolated-evaluation stress comparison
(thousands of proposals, repeated-accept sequences) from attempt 1.

### Performance: a real, reproducible win

Measured NNI acceptance rate on these two cases (uniform-random proposal,
matched to the actual MCMC loop) is high — 72% (`pL1000_200_2`) and 83%
(`MM_BB094_t1_500`) — which is *why* attempt 1's redundant-evaluation-on-accept
problem dominated before `commit_fast` existed: most proposals are accepted,
so "only rejects get cheaper" isn't enough on its own.

`tree_mcmc_cpp_cached_threadsafe_fast` vs. `tree_mcmc_cpp_cached_threadsafe`,
single chain, `RhpcBLASctl` pinned to 1 thread:

| case | iterations | baseline | fast | speedup |
|---|---:|---:|---:|---:|
| `pL1000_200_2` | 4,000 | 2.80s (0.700 ms/iter) | 2.69s (0.672 ms/iter) | 1.04x |
| `pL1000_200_2` | 15,000 | 10.71s (0.714 ms/iter) | 9.57s (0.638 ms/iter) | 1.12x |
| `MM_BB094_t1_500` | 4,000 | 22.57s (5.642 ms/iter) | 12.37s (3.093 ms/iter) | 1.82x |
| `MM_BB094_t1_500` | 15,000 | 69.96s (4.664 ms/iter) | 40.38s (2.692 ms/iter) | 1.73x |

Modest but real on the near-balanced case; substantial on the deep,
imbalanced case that motivated this whole investigation — consistent with
the original diagnosis (depth-dependent cost hurts imbalanced trees most).

## Current status: wired into production, validated at the integration point

Both real production kernels now use the fast path:

- **`SeededTreeChainWorker`** (ordinary MCMC, called via
  `tree_mcmc_parallel_seeded`) now calls `tree_mcmc_cpp_cached_threadsafe_fast`
  instead of `tree_mcmc_cpp_cached_threadsafe`. The latter is kept, unused by
  any driver, purely as a frozen O(depth) reference.
- **`MC3UpdateWorker`** (the flattened, actually-used-in-production MC3
  kernel, called via `tree_mc3_parallel_seeded`) now calls
  `compute_new_loglik_fast`/`commit_fast` instead of
  `compute_new_loglik`/`commit_staged_nni`.

`tree_mc3_cpp_cached_threadsafe` / `tree_mc3_parallel_seeded_serial` (the
serial-temperature validation reference `AGENTS.md` already documents as
existing specifically to stay fixed for comparison) is deliberately left
untouched. `commit_staged_nni`'s eager `sync_outside_cache` option (attempt 1)
is also untouched and still defaults to `false` -- unused now that both real
callers use `commit_fast` instead.

### Integration-point validation

The concern flagged when this was still unwired -- whether the lazy state
holds up under `MC3UpdateWorker`'s concurrent `(ensemble, temperature)` tasks
and the temperature-swap mechanism (`std::swap` on `unique_ptr<NNICache>`
between task slots) -- is now checked, not just reasoned about:

- **Ordinary MCMC, full pipeline**: `scripts/run_regression.sh` against the
  archived `small_test` and `pL1000_200_1` keys -- `mcmc_trace_shape_identical
  = TRUE`, clade-support `MAE = 0`, `max_delta = 0` on both, exact match to
  the archived reference.
- **MC3, raw kernel, before/after**: captured `tree_mc3_parallel_seeded`'s
  full output (traces, final states, swap attempts/accepts) on
  `pL1000_200_2` and `MM_BB094_t1_500` (4 ensembles, 3 temperatures, 3000
  iterations, swap interval 10) immediately before and after wiring
  `MC3UpdateWorker` to the fast path, same seeds. `identical()` on the
  complete result list: `TRUE` on both cases -- the swap mechanism and
  concurrent-task execution do not disturb the lazy state.
- **MC3, full pipeline**: `scripts/benchmark_mc3.R` (EM fit, tree build,
  MC3 sampling, clade annotation, comparison against the archived key) on
  `pL1000_200_2` runs cleanly end to end, ASDSF converging batch over batch
  as expected, support MAE and swap acceptance in the expected range.
- Full package test suite (`R CMD check`, including `test_mc3.R`'s
  checkpoint/resume and deferred-checkpoint tests) passes after each stage.

The performance numbers in the table above now apply directly in production,
since `MC3UpdateWorker`'s per-task inner loop is structurally the same
evaluate/commit pattern as the single-chain loop that was benchmarked.

## What's in the tree

- `G`, `compute_G_vectorized`, `recompute_G_subtree`: attempt 1's outside-cache
  machinery, unchanged, still used to build `G` once at construction.
- `g_content_version`, `g_parent_content_snapshot`, `g_verified_epoch`,
  `g_dirty`, `mark_g_dirty`: lazy `G` invalidation state and marking.
- `f_content_version`, `f_child_content_snapshot0/1`, `f_verified_epoch`,
  `f_dirty`, `mark_f_dirty` (currently unused — see design note above,
  `F`'s cascade needs no explicit marking beyond the direct `commit_fast`
  writes): lazy `F` invalidation state.
- `ensure_g_valid`, `ensure_f_valid`: the O(1)-fast-path, recursive-repair-
  on-miss lazy accessors.
- `compute_new_loglik_fast(edge_n, which)`: evaluates via `G[p1]`, also
  staging `F[p2]`/`F[p1]` into `fast_staged_*` for a cheap commit.
- `commit_fast(new_total_loglik)`: applies a proposal staged by
  `compute_new_loglik_fast`, separate from and non-interfering with
  `commit_staged_nni`.
- `tree_mcmc_cpp_cached_threadsafe_fast`: benchmark-only twin of the
  production single-chain loop, used for the measurements above.
- `nni_cache_compare_eval_cpp`: attempt 1's isolated-evaluation stress test,
  still valid and passing.
- `nni_cache_replay_and_verify_f_cpp` / `debug_verify_f_consistency`: the
  from-scratch consistency checker that found bug 2 above. Kept as a
  debugging tool for any future work on this cache.

All of the above is additive and inert to existing callers: nothing reads
`G`, `f_content_version`, etc. unless it explicitly calls into this new path.
