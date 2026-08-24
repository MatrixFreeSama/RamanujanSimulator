# Precision-doubling Newton materialization

This note validates the arithmetic schedule used by `src/ramanujan_certificate_restore.c` after the precision-ladder rewrite.

The mathematical recovery path is unchanged:

\[
\text{exact seed}\to 32\text{ certified modular transforms}\to
R=uF(z)+v\theta F(z)\to \pi=1/R.
\]

The implementation change is purely arithmetic.  The old reference path initialized every transform at the requested final precision and executed 40 Newton corrections per modular layer.  The current path rebuilds the same 32-layer state sequence at geometrically increasing GMP precisions, reusing the previous stage's modular root as the next stage's starting value.  In the quadratic basin, one Newton correction approximately doubles the number of correct bits.

The final reciprocal is also evaluated with Newton-Schulz refinement,

\[
x_{k+1}=x_k(2-Rx_k),
\]

on a doubled-precision ladder.

## 1000-digit ladder trace

One local C17/GMP validation run produced:

```text
PATH stage bits:       512 -> 1024 -> 2048 -> 3905
PATH Newton totals:    118 ->   67 ->   68 ->   69

RECIP residual digits:
128 bits  ->   38.6 digits
256 bits  ->   77.0 digits
512 bits  ->  154.0 digits
1024 bits ->  307.9 digits
2048 bits ->  615.7 digits
3905 bits -> 1195.0 digits
```

The reciprocal trace is the expected quadratic pattern: the number of correct decimal digits approximately doubles at each precision step until it reaches the arithmetic precision ceiling.

## Local timing comparison

Environment used for this diagnostic: GCC 14.2, GMP 6.3, `-O3 -march=native`, single process.  Times are wall-clock measurements from one cloud CPU environment and are not portable performance guarantees.

| Explicit digits | Fixed-final-precision reference | Precision-doubling Newton | Speedup |
| ---: | ---: | ---: | ---: |
| 1,000 | 0.25 s | 0.03 s | 8.3x |
| 10,000 | 5.91 s | 1.22 s | 4.8x |
| 30,000 | 26.68 s | 4.23 s | 6.3x |
| 100,000 | not rerun | 25.16 s | n/a |

These measurements include the complete 32-transform certificate path.  They should not be compared with a truncated-transform direct generator, because that is a different materialization path.

## Output checks

The 1000-digit precision-ladder output has SHA-256

```text
e898fea26734a6d3af5396b9f4c60ae5dcc88fc40944d835911a9ee8a672ea1b
```

which is the repository's recorded reference hash.

The 30,000-digit output has SHA-256

```text
1f180ef04f63891fff05f0c1b8c9b5dcf6fd0bc97faa5a544c52b71606368bbb
```

and matched the fixed-final-precision path byte-for-byte in the local validation run.

## Interpretation

The speedup does **not** mean that decimal output has become logarithmic in the requested digit count.  Writing \(N\) explicit digits still has an unavoidable \(\Omega(N)\) output cost.  The optimization removes repeated high-precision work inside the recovery calculation by aligning Newton's quadratic convergence with a geometric precision schedule.
