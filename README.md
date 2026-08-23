# RamanujanSimulator

**RamanujanSimulator** is an experimental C17 research implementation for studying extremely deep **implicit precision** in Ramanujan/Chudnovsky-family hypergeometric constructions.

The project separates two costs that are usually coupled:

1. propagating a compact mathematical state into an extremely deep convergence regime; and
2. materializing and storing an enormous explicit representation of \(\pi\).

The active implementation propagates

\[
S_k=(z_k,u_k,v_k)
\]

through a fixed nested sequence of modular transformations. The certifier does **not** write a giant decimal file. Instead, it emits a **constructive certificate** containing the seed, modular-polynomial data, branch rule, state-transport equations, transform sequence, and recovery contract. A separate materializer can replay the same mathematical path at a requested precision and convert the resulting state into an explicit representation.

This repository is research software. It does **not** present an implicit certificate as a conventional stored-decimal world record.

> **Research question.** How far can a precision state be propagated and checked while the full decimal expansion remains unmaterialized, and can the resulting certificate still reconstruct explicit digits without changing the mathematical solve path?

---

## Abstract

High-precision computations of \(\pi\) normally entangle convergence, arbitrary-precision arithmetic, representation, and storage. RamanujanSimulator investigates a different organization of the same mathematical object. It starts from an exact CM/Chudnovsky seed and repeatedly applies modular transformations of orders 2 and 5. Each transformation drives the hypergeometric argument toward zero while transporting the coefficients that represent \(1/\pi\).

The current `10K^4` experiment uses

\[
(2,2,2,2,5,5,5,5)^4,
\]

which contains 32 elementary transforms and has effective modular power

\[
10^{16}.
\]

Under the current validation model, the C17 certifier emits the research-build implicit-depth field

\[
139,999,999,999,999,979
\]

decimal digits while deliberately materializing no full decimal expansion of \(\pi\). The certificate is constructive rather than hash-only. A standalone C reference restorer can read it, replay the same state transport, evaluate the final hypergeometric observer, and emit explicit digits. A 1000-digit recovery smoke test is reproduced byte-for-byte.

The project therefore studies **computation/representation decoupling**. It does not remove the information-theoretic cost of writing \(N\) digits. Instead, it asks whether convergence depth can be represented and checked by a much smaller state before any \(\Theta(N)\) output object is created.

---

## 1. Separation of concerns

| Layer | Object | Responsibility |
| --- | --- | --- |
| Implicit solve | compact state \((z,u,v)\) | propagate the mathematical state |
| Certification | constructive certificate | retain the recipe, validation metadata, and claimed depth |
| Materialization | explicit representation of \(\pi\) | re-evaluate the certificate at requested precision and serialize it |

The intended contract is

\[
\text{same seed}
\rightarrow
\text{same modular transforms}
\rightarrow
\text{same transported state}
\rightarrow
\text{optional observer/storage layer}.
\]

A third party may implement a different storage backend, distributed representation, checkpoint format, decimal serializer, binary limb format, or RNS layout. None of those changes are intended to modify the mathematical solver.

---

## 2. Hypergeometric state

Define

\[
F(z)={}_3F_2\!\left(\frac16,\frac12,\frac56;1,1;z\right)
=\sum_{n=0}^{\infty}c_nz^n,
\]

with

\[
c_n=\frac{(1/6)_n(1/2)_n(5/6)_n}{(n!)^3},
\qquad
\theta=z\frac{d}{dz}.
\]

The working representation is

\[
\frac1\pi=uF(z)+v\,\theta F(z).
\]

The implementation evaluates the coefficients with the exact recurrence

\[
c_{n+1}=c_n\frac{(6n+1)(2n+1)(6n+5)}{72(n+1)^3},
\]

so explicit recovery does not depend on Gamma-function evaluation.

---

## 3. Exact seed

The active branch starts from

\[
z_0=-\frac1{53360^3}=-\frac1{151931373056000},
\]

and

\[
\alpha_0=\frac{77265280}{90856689}.
\]

With

\[
C_0=\frac{\sqrt{163}}6\sqrt{1-z_0},
\]

set

\[
u_0=C_0(1-\alpha_0),
\qquad
v_0=6C_0.
\]

The constructive certificate stores the exact rational/algebraic recipe for this seed rather than a decimal approximation.

---

## 4. Modular state transport

For a modular step of order \(p\), let \(y\) be related to the current argument \(x\) by

\[
P_p(x,y)=0.
\]

The implementation works in the \(z\)-variable after the substitution

\[
J=\frac{1728}{z}
\]

into the corresponding modular polynomial. The active `10K^4` certificate carries the complete integer coefficient tables for \(P_2\) and \(P_5\).

On the selected small branch,

\[
y\sim\frac{x^p}{1728^{p-1}},
\qquad x\rightarrow0.
\]

Implicit differentiation gives

\[
r=\frac{dy}{dx}=-\frac{P_x}{P_y},
\]

\[
s=\frac{d^2y}{dx^2}
=-\frac{P_{xx}+2P_{xy}r+P_{yy}r^2}{P_y}.
\]

The transport factors are

\[
M_p=\frac1p\frac{x}{y}r\sqrt{\frac{1-x}{1-y}},
\]

\[
L_p=1-\frac{xr}{y}+\frac{xs}{r}-\frac{x}{2(1-x)}+\frac{xr}{2(1-y)}.
\]

The state update is

\[
u'=\frac{u-vL_p}{M_p},
\qquad
v'=v\frac{xr}{yM_p}.
\]

Only the current and next states are required. The implementation does not retain an expanding numerical history.

---

## 5. The `10K^4` construction

One macro block is

\[
(2,2,2,2,5,5,5,5),
\]

whose product is

\[
2^4 5^4=10^4.
\]

Repeating the block four times gives

\[
(2,2,2,2,5,5,5,5)^4
\]

with 32 elementary transforms and effective modular power

\[
(10^4)^4=10^{16}.
\]

For a sufficiently small modular parameter the selected branch behaves schematically as

\[
z_{k+1}=O(z_k^{p_k}).
\]

Thus the logarithmic convergence depth scales multiplicatively with the product of the modular orders. The current certifier records the conservative lower bound

\[
-\log_{10}|z_{32}|>140,000,000,000,000,000,
\]

then applies a 21-digit prefactor safety margin and writes

\[
N_{\mathrm{claim}}=139,999,999,999,999,979.
\]

This is a **research-build certificate field** under the current validation model. It is not a claim that this many decimal digits have already been explicitly generated and stored.

---

## 6. Constructive certificate

A useful certificate cannot be only

```text
PASS
precision = ...
hash = ...
```

because a one-way digest cannot reconstruct the mathematical state. The generated JSON contains the information needed for an independent materializer:

- exact seed construction;
- complete 32-step transform sequence;
- full integer coefficient tables for \(P_2\) and \(P_5\);
- small-branch selection rule;
- state-transport equations;
- hypergeometric coefficient recurrence;
- observer relation \(\pi=1/(uF+v\theta F)\);
- current validation metadata and implicit-depth field.

The intended semantics are

\[
\boxed{
\text{certificate}
\rightarrow
\text{replay the same mathematical path}
\rightarrow
\text{explicit observer}
\rightarrow
\text{chosen representation}
}
\]

The official implicit path may discard transient high-precision numerical state after verification while retaining the construction required to reproduce it.

---

## 7. Explicit recovery

`src/ramanujan_certificate_restore.c` is deliberately separate from the certifier. Given a constructive certificate and requested precision, it:

1. reconstructs \(z_0,u_0,v_0\);
2. reads the embedded \(P_2\) and \(P_5\) tables;
3. replays the transform sequence;
4. evaluates \(F(z)\) and \(\theta F(z)\);
5. forms \(R=uF+v\theta F\);
6. evaluates \(\pi=1/R\);
7. serializes the requested explicit digits.

The reference restorer is not intended to dictate how extremely large output must be stored. A large-scale implementation may replace only the final conversion/storage layer with streaming decimal blocks, binary limbs, RNS reconstruction, distributed files, or another representation while preserving the solve path.

If \(N\) decimal digits are actually written, output still has an unavoidable \(\Omega(N)\) information cost. RamanujanSimulator does not claim otherwise.

---

## 8. C17 implementation

The active runtime is C17. GMP provides exact integer/rational primitives and arbitrary-precision floating evaluation. `json-c` is used by the certificate writer and reader.

No Python, SymPy, or mpmath is required by the C runtime.

| File | Purpose |
| --- | --- |
| `src/ramanujan_c_common.c` | shared implementation unit |
| `src/ramanujan_c_common.h` | shared structures and declarations |
| `src/ramanujan_c_common_core.inc` | exact arithmetic, modular-polynomial and state-transport core |
| `src/ramanujan_c_common_qseries.inc` | exact q-series machinery |
| `src/ramanujan_c_common_tail.inc` | remaining shared high-precision helpers |
| `src/ramanujan_nested_core.c` | modular-generation and nested-convergence diagnostics |
| `src/ramanujan_10k4_certifier.c` | `10K^4` certifier and constructive JSON writer |
| `src/ramanujan_certificate_restore.c` | certificate-to-explicit-\(\pi\) reference materializer |

The current mainline is CPU-oriented. The repository does not currently contain CUDA or distributed-HPC code.

---

## 9. Validation model

The repository deliberately distinguishes **finite exact checks**, **numerical diagnostics**, **constructive recovery**, and **theorem-level proof**.

### 9.1 Modular-polynomial checks

The C implementation generates relevant modular polynomials from exact q-series arithmetic using GMP rationals. During the C port, generated \(\Phi_2\), \(\Phi_3\), and \(\Phi_5\) were compared term-for-term with the frozen exact reference tables.

The `10K^4` certifier also evaluates exact rational residuals of

\[
\Phi_p(j(q),j(q^p))
\]

through \(q^{120}\) for \(p=2\) and \(p=5\). The current validation snapshot reports zero residual in both tested finite ranges.

A finite q-series check is implementation evidence. It is **not**, by itself, an infinite-order symbolic proof.

### 9.2 Small-branch and depth checks

The certifier constructs rational envelope bounds in a small-\(q\) region and checks the contraction assumptions used by the conservative depth estimate. Stored decimal digits of \(\pi\) are not used as the source of the claimed convergence depth.

### 9.3 Numerical diagnostics

A separate finite-precision path propagates the 32 states and records selected values of

\[
-\log_{10}|z_k|.
\]

These values are diagnostics rather than the sole basis of the conservative depth field.

### 9.4 Constructive recovery smoke test

The C certificate was used to reconstruct 1000 digits after the decimal point. The explicit output SHA-256 is

```text
e898fea26734a6d3af5396b9f4c60ae5dcc88fc40944d835911a9ee8a672ea1b
```

and matched the earlier Python reference output byte-for-byte.

This demonstrates that the certificate is constructive and that explicit recovery is operational at the tested precision. It does **not** independently prove the entire \(10^{17}\)-scale implicit-depth field.

### 9.5 Formal-proof boundary

The repository is not currently a Lean, Coq, Isabelle, or other proof-assistant formalization. Independent mathematical review remains appropriate for the global transform identities, branch semantics, prefactor/error propagation, and interpretation of the final implicit-depth field.

---

## 10. Reproducibility

### Debian/Ubuntu dependencies

```bash
sudo apt-get update
sudo apt-get install -y build-essential pkg-config libgmp-dev libjson-c-dev
```

### Build

```bash
make
```

The default build enables `-O3 -march=native`. For a portable build:

```bash
make clean
make NATIVE=0
```

### Run nested-core diagnostics

```bash
make core
```

### Generate a fresh certificate

```bash
make certify
```

### Restore 1000 explicit digits

```bash
make restore1000
```

### End-to-end smoke test

```bash
make smoke
```

The smoke test generates a fresh certificate, restores 1000 digits, hashes the output, and compares it with the recorded reference hash. The repository also contains a GitHub Actions workflow for the portable C17 build and smoke test.

---

## 11. Representative C-port measurements

Measurements from one cloud environment during the C17 port are shown below. They are diagnostics, not portable performance guarantees.

| Workload | Python prototype | C17 port |
| --- | ---: | ---: |
| Nested modular core | 2.06 s, ~148948 KB RSS | 0.09 s, ~3608 KB RSS |
| `10K^4` certifier | 2.35 s, ~148664 KB RSS | 0.14 s, ~3904 KB RSS |
| 1000-digit certificate restore | 0.82 s, ~114664 KB RSS | 0.27 s, ~2504 KB RSS |

These timings must not be interpreted as the cost of materializing \(10^{17}\) digits. Huge explicit materialization enters a different regime dominated by multiprecision arithmetic, base conversion, memory traffic, and storage I/O.

---

## 12. Bottleneck inversion

The project is motivated by a representation question. If a computation insists on maintaining an explicit \(N\)-digit object throughout the solve, representational cost becomes entangled with convergence. The nested-state approach instead preserves the mathematical object through a compact generative state and postpones large-scale materialization until observation is requested.

The resulting engineering transition is

\[
\text{convergence bottleneck}
\longrightarrow
\text{materialization and storage bottleneck}.
\]

This is not a claim that information disappears. It is a claim about **when** the large representation is created and **which layer** is responsible for it.

---

## 13. Claim boundary

### This repository does claim

- a concrete C17 implementation of the stated nested modular state transport;
- a 32-step `10K^4` transform sequence with effective modular power \(10^{16}\);
- exact rational finite q-series checks used by the current certifier;
- a constructive certificate format containing the recovery recipe rather than only a hash;
- a separate C materializer that can recover explicit digits from the certificate;
- a reproduced 1000-digit explicit-recovery smoke test;
- the current research-build implicit-depth field exactly as emitted by the certifier.

### This repository does not claim

- that the current implicit-depth field is already a formally verified theorem;
- that 139,999,999,999,999,979 decimal digits have been explicitly generated and stored;
- that an implicit certificate is equivalent to a conventional \(\pi\)-digit world record;
- that explicit materialization at hundreds of trillions of digits takes the same time as the compact certifier;
- that the current branch has undergone independent peer review.

This separation is intentional. The purpose is to make the experiment inspectable without inflating a convergence-scale number into a claim about an output artifact that was never produced.

---

## 14. Independent verification targets

An external reviewer can attack the construction without trusting the executable by independently regenerating \(\Phi_2\) and \(\Phi_5\), checking the conversion to the \(z\)-space polynomials, deriving the \(M_p\) and \(L_p\) transport equations, verifying the small branch over all 32 transforms, replacing finite q-series checks with symbolic arguments, deriving a rigorous global prefactor/error enclosure, and implementing an independent certificate materializer.

The certificate format is deliberately intended to expose enough information for such work.

---

## 15. Design principle

The software rule behind the project is

\[
\boxed{
\text{the numerical representation is an observation of the state, not the state authority itself}
}
\]

The official implicit path therefore ends at a constructive certificate. Explicit decimal output is optional and replaceable. This is the sense in which the project uses an implicit, matrix-free-style representation: a very large representational object is not materialized merely because it can eventually be reconstructed.

---

## 16. Research status

**Status: experimental.**

The repository is intended for reproduction, code review, numerical experimentation, and independent mathematical checking. The active C17 mainline contains the components required for the `10K^4` construction, certificate generation, and explicit recovery path.

---

## 17. License

This repository is distributed under the **MatrixFreeSama Permissive License 2.0 (MFSPL 2.0)**. See [`LICENSE`](LICENSE) for the complete terms.

SPDX custom reference:

```text
LicenseRef-MatrixFreeSama-Permissive-2.0
```

The license text is MFSPL 2.0.

---

## 18. Citation

Machine-readable research-software citation metadata is provided in [`CITATION.cff`](CITATION.cff).

---

## Minimal end-to-end experiment

```bash
make NATIVE=0
make certify
make restore1000
make smoke
```

The experiment demonstrates the repository's central architecture:

\[
\boxed{
\text{compact implicit solve}
\rightarrow
\text{constructive certificate}
\rightarrow
\text{optional explicit recovery}
}
\]

without changing the mathematical solve path between the implicit and explicit modes.
