# Low-order explicit modular pi vs Chudnovsky

This note records an end-to-end explicit-output benchmark for repeated low-order modular transport against Chudnovsky binary splitting.

## Scope

Each method starts from scratch, computes the requested precision, writes the decimal expansion to a file, closes the file, and is timed over the whole operation. The p=2 and p=3 outputs are byte-compared with the Chudnovsky output for every tested precision.

The modular paths use the frozen exact p=2 and p=3 z-polynomials already present in the repository, the common state transport, 14 Newton iterations per transform, and stop adding modular layers once `-log10(|z|)` exceeds the requested digits plus a 160-digit guard.

The benchmark source is `validation/benchmark_low_order_explicit_vs_chud.c` and can be run with `make loworderbench`.

## GitHub Actions measurement

Runner: GitHub-hosted Ubuntu 24.04, GMP 6.3.0, GCC with `-O3 -march=native -std=c17`.

| digits | Chud BS | p=2 layers | p=2 time | p=2 / Chud | p=3 layers | p=3 time | p=3 / Chud | byte equality |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| 1,000 | 0.000158 s | 7 | 0.009484 s | 59.9x | 4 | 0.009442 s | 59.6x | PASS |
| 10,000 | 0.001190 s | 10 | 0.310260 s | 260.6x | 6 | 0.329431 s | 276.7x | PASS |
| 30,000 | 0.004996 s | 11 | 1.539468 s | 308.1x | 7 | 1.748895 s | 350.0x | PASS |
| 100,000 | 0.026690 s | 13 | 9.035463 s | 338.5x | 8 | 9.783033 s | 366.5x | PASS |
| 300,000 | 0.114060 s | 15 | 35.688828 s | 312.9x | 9 | 37.510037 s | 328.9x | PASS |

The corresponding terminal depths at 300,000 digits were approximately 570,788.82 for p=2 and 342,858.70 for p=3.

## Result

No explicit-output crossover against Chudnovsky was observed through 300,000 digits. The low-order modular paths remain roughly two to three orders of magnitude slower in wall time even though they require only O(log N) modular layers.

The earlier implicit-depth sweet-spot result and this explicit-output result answer different questions:

- For fixed arithmetic precision, p=3 gave the highest sampled *implicit convergence-depth gain per second* among p=2,3,5.
- For actual explicit decimal output, p=2 is slightly cheaper than p=3 over the tested range, but neither is close to Chudnovsky binary splitting.

The limiting factor is therefore not the number of convergence layers. It is the cost of each full-precision modular state transport: repeated polynomial/derivative evaluation, Newton root refinement, divisions, square roots, and transport of `(u,v)` at the target precision. Chudnovsky binary splitting spends most of its work in GMP's highly optimized large-integer multiplication tree.

The p=2 and p=3 timings scale somewhat differently from Chudnovsky between individual points, but the present data do not support extrapolating a future crossover. Any such claim would require a substantially cheaper low-order transport kernel or a different arithmetic organization, followed by another end-to-end benchmark.
