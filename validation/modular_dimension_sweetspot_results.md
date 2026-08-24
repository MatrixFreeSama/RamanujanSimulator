# Modular dimension sweet-spot benchmark

This note records one C17/GMP measurement of the frozen modular kernels currently present in the repository. It is an engineering diagnostic, not a theorem about all modular orders or all machines.

## Environment

- CPU: AMD EPYC 9V74
- compiler: GCC 14.2.0
- GMP: 6.3.0
- flags: `-O3 -march=native -std=c17`
- working precision: 6000 decimal digits plus the repository guard conversion
- Newton budget: 14 iterations per transform
- timing: median of 5 repetitions per transform

The benchmark source is `validation/benchmark_modular_dimension_sweetspot.c`.

## Repeated p=2 chain

The purpose of this run is to isolate the relationship between per-step cost and convergence-depth gain while the arithmetic precision is held fixed.

| layer | depth before | depth after | gain | median step time | gain digits/s |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 14.181647 | 31.600839 | 17.419191 | 0.012936 s | 1,346.6 |
| 2 | 31.600839 | 66.439221 | 34.838382 | 0.013023 s | 2,675.1 |
| 3 | 66.439221 | 136.115986 | 69.676765 | 0.012511 s | 5,569.1 |
| 4 | 136.115986 | 275.469515 | 139.353530 | 0.012007 s | 11,606.5 |
| 5 | 275.469515 | 554.176575 | 278.707059 | 0.011976 s | 23,272.6 |
| 6 | 554.176575 | 1111.590693 | 557.414118 | 0.012266 s | 45,443.2 |
| 7 | 1111.590693 | 2226.418930 | 1114.828237 | 0.012313 s | 90,541.6 |
| 8 | 2226.418930 | 4456.075404 | 2229.656474 | 0.012218 s | 182,489.5 |

Across these eight layers the median transform time remains near 12.3 ms while the depth gain doubles exactly from layer to layer. The measured gain-per-second therefore rises by about 135.5x from layer 1 to layer 8. This is the expected signature of the `D -> 2D + constant` structure when the arithmetic precision and kernel cost are held approximately fixed.

This should not be interpreted as explicit decimal-output throughput. It measures implicit convergence depth per wall-clock second.

## Dimension sweep at the common seed

The frozen exact kernels available in the repository are p=2, p=3, and p=5. The first transform from the same seed gives the cleanest apples-to-apples comparison.

| p | polynomial terms | depth gain | median step time | gain digits/s | normalized gain rate `(gain/s)/D_in` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 11 | 17.419191 | 0.012228 s | 1,424.5 | 100.45 s^-1 |
| 3 | 17 | 34.838382 | 0.021606 s | **1,612.5** | **113.70 s^-1** |
| 5 | 38 | 69.676765 | 0.057515 s | 1,211.5 | 85.42 s^-1 |

Relative to p=2:

- p=3 provides 2x the depth gain for about 1.77x the step cost, so gain-per-second improves by about 13.2%.
- p=5 provides 4x the depth gain for about 4.70x the step cost, so gain-per-second falls to about 85.0% of p=2.

Under this implementation and this machine, the current sampled sweet spot is therefore around **p=3**, with p=2 close behind and p=5 already on the cost-dominated side.

Because only p=2, p=3, and p=5 frozen kernels are currently implemented, this result should be stated as a sampled sweet spot, not as proof that order 3 is globally optimal. A direct p=4 kernel or other orders could move the optimum.

## Interpretation

The present bottleneck is consistent with a high single-step cost competing against a high single-step convergence reward. The p=2 chain demonstrates that implicit depth-per-second can grow rapidly when step cost remains nearly flat. The dimension sweep shows that raising p is useful only while the extra depth multiplier grows faster than the modular-polynomial evaluation and transport cost.
