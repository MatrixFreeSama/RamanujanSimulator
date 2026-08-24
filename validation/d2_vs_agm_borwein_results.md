# Original D2 vs AGM and Borwein quartic

This benchmark compares the original class-number-2 direct evaluator with two classical iterative pi algorithms under one explicit-output protocol.

## Methods

- Original D2 direct evaluator: `d = h(Delta) = 2`, `Delta = -427`, about 24.9558996577 decimal digits per hypergeometric term.
- Gauss-Legendre AGM: quadratic convergence.
- Borwein quartic algorithm: quartic convergence.
- Chudnovsky binary splitting: reference implementation and byte-equality oracle.

All methods use GMP `mpf`, the same decimal precision budget (`digits + 192` decimal guard converted to bits), the same decimal serializer, and end-to-end timing that includes initialization, computation, serialization, and file close. AGM and Borwein are deliberately run at full target precision for every iteration; no precision-doubling schedule is used, so their timings are conservative rather than aggressively optimized.

Every D2, AGM, and Borwein output was byte-identical to the Chudnovsky binary-splitting output.

## GitHub Actions measurement

Runner: GitHub-hosted Ubuntu 24.04, GMP 6.3.0, GCC C17, `-O3 -march=native`.

| digits | Chud BS | D2 direct | AGM | AGM iters | Borwein-4 | Borwein iters | D2 / AGM | D2 / Borwein-4 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 0.000103 s | 0.000124 s | 0.000075 s | 15 | 0.000121 s | 9 | 1.65x | 1.02x |
| 10,000 | 0.000674 s | 0.013888 s | 0.001785 s | 18 | 0.002511 s | 10 | 7.78x | 5.53x |
| 30,000 | 0.002906 s | 0.173707 s | 0.009015 s | 19 | 0.013672 s | 11 | 19.27x | 12.71x |
| 100,000 | 0.015599 s | 2.961271 s | 0.053192 s | 21 | 0.080243 s | 12 | 55.67x | 36.90x |
| 200,000 | 0.041018 s | 14.568708 s | 0.133629 s | 22 | 0.191821 s | 12 | 109.02x | 75.95x |

## Result

The original D2 construction beats the original Chudnovsky formula when both are evaluated as the same direct hypergeometric recurrence, as recorded separately in `original_d2_basic_vs_chud_results.md`. It does not beat Gauss-Legendre AGM or the Borwein quartic iteration as a complete explicit-output algorithm in this direct-recurrence implementation.

The reason is structural. D2 improves a linear-in-term-count direct series by increasing the digits gained per term from about 14.18 to about 24.96. AGM and Borwein instead multiply the number of correct digits per iteration, so the number of high-precision iterations grows only logarithmically with target precision. The D2 per-term convergence advantage therefore cannot compensate for thousands of full-precision recurrence steps at large `N`.

At 200,000 digits, D2 uses about 8,025 direct series terms, while AGM uses 22 full-precision iterations and Borwein quartic uses 12. Even without adaptive precision, the iterative methods are respectively about 109x and 76x faster than D2 direct on this runner.

This benchmark does not answer whether the D2 algebraic series can compete after it receives an asymptotically efficient product-tree or binary-splitting evaluator. That remains a separate arithmetic-organization question. The present result only establishes that the basic D2 direct recurrence is not the overall speed leader among classical high-order pi algorithms.
