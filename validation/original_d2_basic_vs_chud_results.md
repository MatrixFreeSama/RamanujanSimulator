# Original class-number-2 evaluator vs Chudnovsky

This benchmark returns to the earliest low-dimensional plug-in family, before the later modular-polynomial state-transport machinery.

## Exact D2 member

The dimension convention is

`d = h(Delta)`

with the class-number-2 member

`Delta = -427`, `d = 2`.

Its frozen exact algebraic data are the two quadratic polynomials

`H(J) = J^2 + 15611455512523783919812608000 J + 155041756222618916546936832000000`

and

`L(U) = U^2 + 6049980860956530737897555251200 U + 9989238519497195119784714748139929600`.

The evaluator takes the tiny reciprocal roots `y = 1/J` and `v = 1/U` directly from the stable quadratic formula, then uses

`z = 1728 y`

and

`alpha = y / (427 v (1 - 1728 y))`.

The explicit series is

`1/pi = sqrt(427)/6 * sqrt(1-z) * sum c_n (6n + 1 - alpha) z^n`

with

`c_n = (1/6)_n (1/2)_n (5/6)_n / (n!)^3`.

No modular polynomial, no `(z,u,v)` state transport, and no repeated modular Newton solve is used here.

The D2 convergence scale is approximately

`-log10(|z|) = 24.9558996576542667...` decimal digits per series term.

For comparison, the D=163 class-number-1 seed is algebraically Chudnovsky and gives approximately

`14.1816474627254777...` decimal digits per term.

The asymptotic term-count ratio is therefore about

`14.1816474627 / 24.9558996577 = 0.568...`,

so D2 needs only about 56.8% as many direct-series terms.

## Benchmark protocol

GitHub-hosted Ubuntu 24.04 runner, GMP 6.3.0, GCC C17, `-O3 -march=native`.

Every timing is end-to-end for one run: initialize algebraic constants, compute pi, serialize the requested decimal digits, and close the file.

Three paths are compared:

1. Chudnovsky binary splitting.
2. The D=163/Chudnovsky formula evaluated by the same direct hypergeometric recurrence used by D2.
3. The original D2 quadratic plug-in evaluated by that direct recurrence.

Every D2 output was byte-identical to both Chudnovsky outputs.

| digits | Chud BS | Chud direct | D2 direct | D2 / Chud direct | D2 / Chud BS | D2 terms / Chud direct terms |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 0.000153 s | 0.000284 s | 0.000203 s | 0.714x | 1.33x | 50 / 86 |
| 10,000 | 0.001178 s | 0.040243 s | 0.023710 s | 0.589x | 20.14x | 410 / 720 |
| 30,000 | 0.004964 s | 0.522281 s | 0.299104 s | 0.573x | 60.25x | 1212 / 2131 |
| 100,000 | 0.026548 s | 8.469281 s | 4.835685 s | 0.571x | 182.15x | 4017 / 7067 |
| 200,000 | 0.070367 s | 43.400008 s | 24.699853 s | 0.569x | 351.01x | 8024 / 14118 |

## Result

The original D2 construction does show a real formula-level speed advantage when both formulas use the same direct recurrence. At 200,000 digits it takes about 56.9% of the Chudnovsky-direct time, corresponding to a speedup of about 1.76x. The measured runtime ratio approaches the predicted term-count ratio.

This means the low-dimensional plug-in itself is not the bottleneck. Its quadratic algebraic setup becomes negligible, and its higher per-term convergence translates almost one-for-one into wall-time improvement under matched arithmetic organization.

However, this direct D2 implementation does not beat optimized Chudnovsky binary splitting. The reason is arithmetic organization rather than convergence rate: Chudnovsky binary splitting replaces a long chain of full-precision floating-point recurrence updates with a balanced exact-integer multiplication tree. The gap therefore grows with requested precision in this benchmark.

The next meaningful comparison is not another modular transport benchmark. It is a binary-splitting or equivalent product-tree evaluator for the D2 algebraic series itself. Only that comparison can determine whether the roughly 1.76x convergence advantage survives after both sides use asymptotically efficient summation.
