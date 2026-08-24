# D2 binary-splitting noise-floor study

This validation pass separates small implementation overhead from cloud-runner timing noise for the accelerated D2 path. It does not execute the bare direct D2 or bare direct Chudnovsky series.

## Change under test

The D2 quadratic-ring product tree now allocates its leaf and ring-multiplication scratch `mpz_t` objects once per complete solve and reuses them across the recursive traversal. This removes repeated GMP `init/clear` traffic without changing the arithmetic, term count, algebraic basis, binary-splitting topology, or requested precision.

## Timing protocol

- GitHub-hosted Ubuntu 24.04 runner.
- GCC C17 with `-O3 -march=native`, GMP 6.3.0.
- Baseline source is fetched from `main` and compiled beside the optimized source on the same runner.
- Output files are written on `/dev/shm`, retaining decimal conversion and memory-copy cost while removing persistent-filesystem latency.
- One warm-up solve is run before sampling.
- Nine samples are collected for each precision and each variant.
- Baseline/optimized execution order alternates on consecutive samples.
- Results are summarized by the median and median absolute deviation (MAD).
- Every individual solve retains the existing byte-equality check against its Chudnovsky-BS reference.
- The complete benchmark job was run twice on independent hosted runners.

## Independent run 1

| digits | baseline D2 median | optimized D2 median | optimized reduction | optimized D2 / Chud median | optimized D2 MAD |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100,000 | 0.039278269 s | 0.038039684 s | 3.15% | 1.556638x | 0.000279665 s |
| 300,000 | 0.168915272 s | 0.166069984 s | 1.68% | 1.609413x | 0.000348091 s |
| 1,000,000 | 0.754058123 s | 0.743612289 s | 1.39% | 1.542578x | 0.001044270 s |

## Independent run 2

| digits | baseline D2 median | optimized D2 median | optimized reduction | optimized D2 / Chud median | optimized D2 MAD |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100,000 | 0.039191723 s | 0.038138390 s | 2.69% | 1.565030x | 0.000212908 s |
| 300,000 | 0.169033289 s | 0.166160345 s | 1.70% | 1.610780x | 0.000296354 s |
| 1,000,000 | 0.754204750 s | 0.746591806 s | 1.01% | 1.544237x | 0.000851154 s |

## Combined interpretation

Averaging the two independent median measurements gives an implementation-overhead reduction of about 2.92% at 100k digits, 1.69% at 300k digits, and 1.20% at one million digits.

The corresponding optimized D2/Chud median ratios from the two runs cluster tightly around:

- about 1.56x at 100k digits,
- about 1.61x at 300k digits,
- about 1.54x at one million digits.

The one-million-digit D2 median MAD is about 0.14% in both runs. The paired D2/Chud ratio MAD is below about 0.4% there. Therefore the remaining roughly 1.54x high-precision gap is much larger than the measured timing noise and should be treated as a real arithmetic constant-factor gap, not a startup or filesystem artifact.

The reusable-scratch change is retained because it gives a repeatable low-single-digit speedup and reduces avoidable allocation traffic. Further progress should target the number and size of quadratic-ring multiplications rather than additional benchmark warm-up or file-I/O tuning.
