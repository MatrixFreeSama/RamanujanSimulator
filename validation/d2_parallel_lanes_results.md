# Exact two-lane D2 binary splitting

This experiment tests the original computational interpretation that a two-dimensional algebraic state should expose two concurrent one-dimensional coordinate lanes rather than forcing both coordinates through one serial ring kernel.

## Exact lane decomposition

Use the integral basis

`omega = (1 + sqrt(61))/2`, `omega^2 = omega + 15`.

For one quadratic-ring product

`(a + b*omega)(c + d*omega)`

the two output coordinates are

`A = ac + 15 bd`

and

`B = (a+b)(c+d) - ac`.

For large operands, the implementation assigns `A` and `B` to two OpenMP workers. The `ac` multiplication is deliberately duplicated, so the two coordinate jobs do not need an inter-lane barrier while the large products are running. For small operands the existing three-product serial kernel is retained because thread scheduling would cost more than it saves.

This remains exact `mpz_t` arithmetic and preserves the same D2 binary-splitting tree, term count, algebraic basis, guard precision, and final byte-level validation.

A mathematical boundary is important: the two basis coordinates are not two completely uncoupled scalar recurrences. Ring multiplication couples the input coordinates. The implementation therefore parallelizes the two exact *output-coordinate jobs* rather than pretending that the quadratic field is a direct product of two independent integer rings.

The alternative conjugate-embedding diagonalization is also not a symmetric replacement. The physical branch has

`-log10(|z_plus|) = 24.9558996576542667...`,

while the other conjugate has only

`-log10(|z_minus|) = 0.7594615452292167...`.

So evaluating two independent conjugate hypergeometric series would make the second lane vastly slower and would throw away the convergence advantage. The exact coordinate-lane schedule avoids that trap.

## Protocol

- GitHub-hosted Ubuntu 24.04 runners.
- GMP 6.3.0, GCC C17, `-O3 -march=native -fopenmp`.
- D2 is restricted to two OpenMP threads with `OMP_NUM_THREADS=2`, `OMP_DYNAMIC=false`, `OMP_PROC_BIND=true`, `OMP_PLACES=cores`.
- Chudnovsky-BS remains the existing one-thread one-dimensional reference.
- Output is written on tmpfs during the threshold sweep.
- Five samples per precision and threshold; table values are medians.
- Every individual D2 result is byte-identical to the Chudnovsky-BS reference.
- The complete threshold sweep was run twice on independent hosted runners.

## Threshold sweep

The best stable threshold is approximately 4096 operand bits. Parallelizing substantially smaller products adds scheduling overhead; waiting until tens or hundreds of thousands of bits leaves useful two-lane work serialized.

### Independent run 1

| threshold bits | 100k D2/Chud | 300k D2/Chud | 1M D2/Chud |
| ---: | ---: | ---: | ---: |
| 2,048 | 1.4442x | 1.4005x | 1.3281x |
| **4,096** | **1.4232x** | **1.3920x** | **1.3170x** |
| 8,192 | 1.4251x | 1.3973x | 1.3202x |
| 16,384 | 1.4434x | 1.4049x | 1.3280x |
| 32,768 | 1.4646x | 1.4182x | 1.3435x |
| 65,536 | 1.4903x | 1.4441x | 1.3631x |

At one million digits and threshold 4096:

- Chudnovsky-BS: `0.461914539 s`
- exact two-lane D2-BS: `0.608612299 s`
- ratio: `1.316987x`

### Independent run 2

| threshold bits | 100k D2/Chud | 300k D2/Chud | 1M D2/Chud |
| ---: | ---: | ---: | ---: |
| 2,048 | 1.3953x | 1.3686x | 1.2963x |
| **4,096** | **1.3884x** | **1.3664x** | **1.2924x** |
| 8,192 | 1.3987x | 1.3688x | 1.2954x |
| 16,384 | 1.4236x | 1.3794x | 1.3120x |
| 32,768 | 1.4429x | 1.4002x | 1.3254x |
| 65,536 | 1.4760x | 1.4234x | 1.3509x |

At one million digits and threshold 4096:

- Chudnovsky-BS: `0.481166124 s`
- exact two-lane D2-BS: `0.622991562 s`
- ratio: `1.292441x`

## Interpretation

The previous exact serial quadratic-ring D2 path clustered near a `1.54x` high-precision ratio against Chudnovsky-BS. The exact two-lane coordinate schedule reduces the one-million-digit ratio to approximately `1.29x` to `1.32x` on two independent runners.

So the dimensional parallelism is real and materially useful, but it does not halve the total wall time. The reason is structural: only the two algebraic coordinate products and coordinate accumulation can be split this way. The common denominator tree, leaf construction, final algebraic observation, decimal conversion, and the coupling inherent in quadratic multiplication remain shared work.

The important result is therefore not `2D = exactly 2x speedup`. It is:

`2D exposes two exact concurrent coordinate jobs, reducing the remaining D2/Chud constant-factor gap from about 1.54x to about 1.30x at one million digits.`

This comparison intentionally leaves Chudnovsky as the one-thread one-dimensional reference. Chudnovsky's binary-splitting tree can also be parallelized by range decomposition, so these measurements should be interpreted as a test of the requested dimension-to-thread execution model, not as a universal claim about the fastest possible multithreaded Chudnovsky implementation.
