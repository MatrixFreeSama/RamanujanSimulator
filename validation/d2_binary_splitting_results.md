# D2 binary splitting vs classical pi algorithms

This benchmark removes the bare direct-series comparison from the active timing path. The D2 formula is accelerated with a balanced exact product tree, directly analogous to Chudnovsky binary splitting.

## Algebraic basis

For the selected class-number-2 branch (`Delta = -427`), use the integral quadratic basis

`omega = (1 + sqrt(61))/2`, so `omega^2 = omega + 15`.

The two algebraic quantities required by the explicit formula reduce to

`z = (-59818419102592333 + 13579278976889262*omega) / 609541351191872000`

and

`beta = 1-alpha = (320004750671 - 56552805760*omega) / 766923569391`.

This basis is materially cheaper than the earlier large-coefficient `x = H2/J` basis. A quadratic-ring multiplication

`(a+b*omega)(c+d*omega)`

needs three large products:

`m0 = a*c`, `m1 = b*d`, `m2 = (a+b)(c+d)`,

followed by

`constant = m0 + 15*m1`, `omega_coeff = m2 - m0`.

The binary-splitting tree carries only `Q`, ratio product `P`, and the already-weighted partial sum `W`, so the structure is the quadratic-ring analogue of the usual Chudnovsky `P/Q/T` tree.

## Protocol

GitHub-hosted Ubuntu 24.04, GMP 6.3.0, GCC C17, `-O3 -march=native`.

The D2-BS and Chudnovsky-BS timings are end-to-end and include initialization, computation, decimal serialization, and file close. AGM and Borwein quartic are run immediately afterward on the same runner with the same precision budget and serializer. Their outputs are compared against the exact Chudnovsky-BS reference files produced in the first phase.

The final comparison does not execute the bare direct D2 or bare direct Chudnovsky series.

AGM and Borwein quartic use full target precision for every iteration; no adaptive precision-doubling schedule is used, so their timings are conservative rather than maximally optimized.

Every output at every tested precision was byte-identical to the Chudnovsky-BS reference.

| digits | Chud BS | D2 BS | D2 / Chud | AGM | Borwein-4 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 0.000127 s | 0.000139 s | 1.09x | 0.000165 s | 0.000146 s |
| 10,000 | 0.001029 s | 0.001718 s | 1.67x | 0.002711 s | 0.003610 s |
| 30,000 | 0.004262 s | 0.007249 s | 1.70x | 0.013150 s | 0.019744 s |
| 100,000 | 0.022866 s | 0.037463 s | 1.64x | 0.076566 s | 0.115266 s |
| 200,000 | 0.059652 s | 0.095458 s | 1.60x | 0.201342 s | 0.280549 s |
| 300,000 | 0.104132 s | 0.163471 s | 1.57x | 0.326274 s | 0.471700 s |
| 1,000,000 | 0.464147 s | 0.736330 s | 1.59x | 1.463959 s | 2.043881 s |

Term counts for the one-million-digit run were 70,530 for Chudnovsky and 40,082 for D2, preserving the approximately 0.568 term-count ratio implied by the formula-level convergence constants.

## Result

After receiving the same class of product-tree acceleration, D2 is no longer hundreds of times slower than Chudnovsky. At one million digits it is about 1.59x slower than the current Chudnovsky binary-splitting implementation while using only about 56.8% as many series terms.

On this same-run benchmark, D2-BS is faster than the full-precision Gauss-Legendre AGM and Borwein quartic implementations at every tested precision above the tiny startup regime. At one million digits:

- D2-BS / AGM = about 0.503, so D2-BS is about 1.99x faster.
- D2-BS / Borwein-4 = about 0.360, so D2-BS is about 2.78x faster.

These AGM/Borwein comparisons are implementation results, not a claim that D2-BS beats every optimized variant of those algorithms. Adaptive-precision implementations can reduce their cost.

The remaining gap to Chudnovsky is therefore no longer caused by a missing asymptotically efficient summation scheme. It is now primarily the constant-factor cost of performing the product tree in a quadratic coefficient ring instead of the ordinary integer ring. The current measured gap is approximately 1.6x at 100k-1M digits.
