# Convergence-tapered direct recurrence

This experiment deliberately stays on the original scalar direct-recurrence path for the class-number-2 member `Delta = -427`. It does not replace the series by binary splitting and it does not carry an exact quadratic-field product tree through the summation.

## Accelerator

For the direct hypergeometric term recurrence, the magnitude of term `n` falls approximately as

`|t_n| ~ 10^(-q n)`.

For D2,

`q = 24.9558996576542667...` decimal digits per term.

For the D=163 / Chudnovsky direct recurrence,

`q = 14.1816474627254777...` decimal digits per term.

The Convergence-Tapered Direct Recurrence (CTDR) therefore reduces the relative working precision of the recurrence as the terms become smaller. In the implementation, the active bit precision is approximately

`P_bits - q*n*log2(10) + safety_bits`,

with a conservative safety allowance and a minimum precision floor. Precision changes are made only at coarse term-block boundaries to reduce precision-management overhead. The full-precision `F` and `T = sum(n*t_n)` accumulators are retained.

The formula, term count, term recurrence, and final observation are otherwise unchanged. The same accelerator is applied to D2 and to the one-dimensional D=163 / Chudnovsky direct recurrence.

## Rejected variant

An earlier magnitude-window prototype also reduced the precision of local block accumulators. It was faster, but byte validation failed. The reason is structural: although the magnitude of a later term is small, its mantissa still contributes information across the remaining requested digits. Limiting a local accumulator to only the inter-term magnitude span discards deep digits. That variant is rejected and is not the reported CTDR method.

## Validation

GitHub-hosted Ubuntu 24.04 runners, GMP, GCC C17, `-O3 -march=native`.

Every reported D2 CTDR and Chud CTDR output is byte-identical to the Chudnovsky binary-splitting reference at the requested decimal precision.

The high-precision validation run used GitHub Actions run `32731870165`, job `97445661238`.

| digits | taper block | Chud BS | D2 direct | D2 CTDR | Chud direct | Chud CTDR | D2 CTDR / D2 direct | D2 CTDR / Chud CTDR |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 100,000 | 16 | 0.024451971 s | 4.376951456 s | 1.952798843 s | 7.641611338 s | 3.396393776 s | 0.446155x | 0.574962x |
| 100,000 | 64 | 0.024074316 s | 4.366041899 s | 1.973523617 s | 7.648643017 s | 3.419784546 s | 0.452017x | 0.577090x |
| 200,000 | 16 | 0.062802076 s | 22.369212151 s | 9.841870785 s | 39.258158445 s | 17.253968477 s | 0.439974x | 0.570412x |
| 200,000 | 64 | 0.062806129 s | 22.368253231 s | 9.981220245 s | 39.240109205 s | 17.294273138 s | 0.446223x | 0.577140x |

At 200,000 digits with a 16-term taper block:

- CTDR reduces D2 direct time from `22.369212151 s` to `9.841870785 s`, about `2.27x` faster.
- CTDR reduces Chudnovsky direct time from `39.258158445 s` to `17.253968477 s`.
- Under the same CTDR accelerator, `D2 / Chud = 0.570412x`, so D2 is about `1.75x` faster than the equally accelerated Chudnovsky direct recurrence.

The measured `0.5704x` ratio is close to the formula-level asymptotic term-count ratio

`14.1816474627 / 24.9558996577 ~= 0.568`.

This shows that the original D2 formula-level convergence advantage continues to translate almost one-for-one into wall-clock advantage when both direct recurrences receive the same convergence-aware precision accelerator.

## Boundary

CTDR is an accelerator for the direct recurrence state. It is not competitive with optimized Chudnovsky binary splitting at these precisions. At 200,000 digits in the same run, Chudnovsky BS took about `0.0628 s`, while D2 CTDR took about `9.84 s`.

Therefore the result supports two narrower claims:

1. Precision tapering is a substantial valid acceleration of the original direct D2 recurrence.
2. D2 retains approximately its full formula-level advantage over Chudnovsky when the same direct-state accelerator is applied to both.

It does not establish superiority over binary-splitting Chudnovsky.
