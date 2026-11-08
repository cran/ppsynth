# Certification materials

`fixtures/*.json` are the Python reference engine's outputs
(ppsynth_core.py, PPSYNTH v2.6.1, numpy 2.5.3) on five designs, written by
`tools/make_fixtures.py` with shortest round-trip doubles (so the JSON
carries the reference doubles exactly; `Inf` is encoded as the string
"inf").

| file | design | fits / inference |
|---|---|---|
| germany_fits.json | Germany panel, trunit 8, trperiod 1990 | aggregate; cv; cv holdout 3; lambda 0, 0.3, 3, Inf |
| testpanel.json | shipped 9-unit test panel, trunit 1, trperiod 2010 | aggregate; cv; inference: aggregate (both conventions), cv, lambda 1 |
| selftest_toy.json | `_ppsynth_selftest` design (seed 7) | aggregate; monotone-grid fits |
| selftest_guardrail.json | `_ppsynth_selftest_guardrail` design (seed 11) | aggregate; cv |
| selftest_support.json | `_ppsynth_selftest_support` designs (seed 7) | main; boundary; single |

`results/` holds the gate output: one CSV per fixture set
(`gate_*.csv`: fixture, fit, field, max absolute deviation, note, R
seconds), the union `certification_results.csv`, the machine description
`environment.txt`, and `certification_report.pdf`.

The package tests (`tests/testthat/test-certification.R`) re-run the small
fixtures at every check; set `PPSYNTH_FULL_GATE=1` to also run the Germany
fits and the placebo-inference fixtures (slow).
