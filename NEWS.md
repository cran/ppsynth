# ppsynth 0.2.1

* CRAN resubmission. The placebo-inference example in `?ppsynth` is
  `\donttest{}` rather than `\dontrun{}`. `ppsynth_selftest()` no longer
  calls `set.seed()`: its three toy designs are built from a fixed internal
  number sequence, so the function neither sets nor disturbs the user's
  random-number generator state.

# ppsynth 0.2.0

* Compiled inner solver (`engine = "cpp"`, the default), mirroring the numba
  engine behind the Stata command; the pure-R solver stays available as
  `engine = "r"`. Both are certified against the Python engines (see
  `inst/certification`).
* `cores` argument: the leave-one-out folds and the placebo units run in
  parallel; results do not depend on `cores`.
* Package renamed from `ppsynthR` to `ppsynth`.
* Germany placebo-inference fixture added to the certification set.

# ppsynthR 0.1.0

* First release: pure-R port of `ppsynth_core.py` (PPSYNTH v2.6.1), the
  Stata command's option set, display and graphs, certified against the
  Python engine on the Germany panel, the shipped test panel and the three
  self-test designs.
