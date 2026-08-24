# Class-number-10 single-term convergence certificate

This validation isolates the formula-level convergence question. It does not benchmark binary splitting and it does not claim a globally fastest finite `1/pi` series.

## D10 choice

For fundamental imaginary quadratic discriminants with class number `h = 10`, the largest absolute discriminant is

`D = -13843`.

The C17 validator independently enumerates reduced primitive positive binary quadratic forms and reproduces

- `h(-13843) = 10`,
- 87 fundamental discriminants of class number 10,
- `13843` as the largest fundamental absolute discriminant in that class.

The same scan gives the following maxima through class number 10:

| class number h | largest fundamental |D| |
| ---: | ---: |
| 1 | 163 |
| 2 | 427 |
| 3 | 907 |
| 4 | 1555 |
| 5 | 2683 |
| 6 | 3763 |
| 7 | 5923 |
| 8 | 6307 |
| 9 | 10627 |
| 10 | 13843 |

## Exact algebraic convergence carrier

Let `C = j^1/3` for the principal CM embedding. For `D = -13843`, the selected degree-10 factor is

`P10(C) = 0`,

with

`P10(x) = x^10`

`+ 322959704590412332725686274638226069771593878764442080 x^9`

`+ 53975830800870124552573137895158808579232875400442163200 x^8`

`+ 17531495856827170677348882327567868131596697367321055232000 x^7`

`+ 1234312959090914569645822713213259492862969867888586915840000 x^6`

`+ 33483032339905697583284254695534390058767103248033329971200000 x^5`

`- 211565635591557038546781996408514792373254059354608369664000000 x^4`

`+ 421486683565974884866530269709770250139307349495722803200000000 x^3`

`- 63536734918606054543228298377769047547857440465492954316800000000 x^2`

`+ 148269606635177168317423665367562949221038833187015360512000000000 x`

`- 718150912949642523148693959180027408505991659454736629760000000000`.

The validator stores an integer Hilbert class polynomial `H_{-13843}(j)` and checks by exact GMP polynomial division that

`P10(x)` divides `H_{-13843}(x^3)`.

The principal real root has

`|C| = 3.22959704590412297e53`,

so `|j| = |C|^3` and the factorial kernel contributes an asymptotic factor `1728`. The resulting asymptotic decimal gain per term is therefore

`log10(|j| / 1728) = 157.289901279323630 digits/term`.

This corresponds to about

- 636 terms for 100,000 decimal digits,
- 1,272 terms for 200,000 decimal digits,
- 6,358 terms for 1,000,000 decimal digits,

before guard terms and implementation details.

For comparison, the class-number-1 Chudnovsky member is about `14.181647` digits/term, the class-number-2 `D=-427` member is about `24.955900` digits/term, and the classical class-number-4 `D=-1555` example is about `50.56` digits/term. Thus this D10 member has substantially higher single-term convergence than those lower-class examples.

## Reproduction

Run

```text
make d10convergencebench
```

The GitHub Actions validation run `32738231741`, job `97466158984`, completed successfully and reported:

```text
D10 convergence certificate
discriminant=-13843 class_number=10 reduced_forms=10
exact_factor_H_of_x3_by_C_polynomial=PASS
principal_abs_C=3.22959704590412297e+53
asymptotic_digits_per_term=157.289901279323630
terms_for_100k=636
terms_for_200k=1272
terms_for_1M=6358
```

## Claim boundary

This result supports the narrower statement:

**Among fundamental imaginary quadratic fields of class number 10, `D=-13843` is the maximum-discriminant member and its associated Ramanujan/Borwein-type series has about 157.29 asymptotic decimal digits of convergence per term.**

It does **not** establish a global world record for single-term convergence. The general Ramanujan/Borwein construction admits arbitrarily higher class number/discriminant and therefore increasingly rapid convergence at the price of increasingly complicated algebraic constants. There is no finite unrestricted champion under the metric `digits per term` alone.

This certificate also does not yet supply the full non-circular degree-10 algebraic constants `a(t)` and `b(t)` needed for a standalone D10 explicit pi evaluator. The present artifact certifies the class number, the exact CM `j^(1/3)` carrier, and the resulting single-term convergence rate.
