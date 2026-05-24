# TODO

Open work items for the v1 rewrite. Each entry names the failure mode and where it
shows up so the fix can be verified end to end.

## Examples status

- `examples/retrodiction_homodyne.jl`: file runs end-to-end, but the
  docs build entry in `docs/make.jl` is still commented out to match
  master. Origin: PR #266 shipped it commented; no recorded reason.
  Suspected cause is the 5 × 10⁴-step SDE solve being too heavy for the
  docs build. Decide whether to uncomment after a build-time check.
