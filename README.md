# RamanujanSimulator

RamanujanSimulator is an experimental research implementation for extremely deep implicit precision in Ramanujan/Chudnovsky-family hypergeometric constructions.

The project separates two tasks that are usually coupled:

1. **Implicit solve and certification**: propagate a compact mathematical state and certify an extremely small final convergence parameter without materializing the full decimal expansion of \(\pi\).
2. **Optional explicit materialization**: a third party may replay the same constructive certificate at arbitrary precision and convert the result into a conventional decimal representation. The mathematical solve path does not need to be changed; only the representation/storage layer is added.

The current runtime has been ported from the original Python research prototype to **C17**.

## Core idea

The active branch starts from the exact CM/Chudnovsky seed

```text
z0 = -1 / 53360^3
```

and transports the state

```text
(z, u, v)
```

through modular transformations of orders 2 and 5 while preserving the hypergeometric representation

```text
1/pi = u F(z) + v theta(F)(z)
```

with

```text
F(z) = 3F2(1/6, 1/2, 5/6; 1, 1; z).
```

The current 10K^4 experiment uses the 32-step sequence

```text
(2, 2, 2, 2, 5, 5, 5, 5)^4
```

whose effective modular power is

```text
10^16
```

The certifier does **not** write a giant decimal file. It records the exact seed, transformation sequence, modular-polynomial data, branch rule, state-transport equations, and recovery contract in a constructive certificate.

## Constructive certificate

The certificate is intended to be more than a hash or a precision claim. It contains enough mathematical information for an independent implementation to reconstruct the same state trajectory and then materialize an explicit value if desired.

Conceptually:

```text
seed
  -> nested modular state transport
  -> implicit certified state
  -> constructive certificate
  -> optional third-party materializer
  -> explicit decimal pi
```

The official implicit path stops after certification and may discard transient high-precision numerical state. Decimal serialization, distributed storage, checkpointing, and large-scale I/O are deliberately outside the solver core.

## C17 implementation

The C port is organized around the following components:

| Component | Purpose |
| --- | --- |
| `ramanujan_nested_core.c` | Runtime generation of modular relations and nested state transport |
| `ramanujan_10k4_certifier.c` | 10K^4 implicit-depth checks and constructive certificate generation |
| `ramanujan_certificate_restore.c` | Reference explicit recovery from a constructive certificate |
| `ramanujan_d100_evaluator.c` | Legacy D100/class-number experiment |
| `ramanujan_c_common.c/.h` | Shared GMP exact arithmetic, modular-polynomial generation, state transport, hypergeometric evaluation |

Dependencies:

```text
C17 compiler
GMP
json-c     # certifier/restorer only
```

No Python, SymPy, or mpmath is required by the C runtime.

## Build

```bash
make
```

For a more portable build without `-march=native`:

```bash
make clean
make NATIVE=0
```

Typical commands:

```bash
./ramanujan_core
./ramanujan_certifier
./ramanujan_restore Ramanujan_10K4_Implicit_Certificate_C.json --digits 1000 --guard 140 --output pi1000.txt
```

## Current validation snapshot

In the current research build:

```text
10K^4 elementary transforms: 32
effective modular power:      10^16
decimal pi materialized by certifier: NO
```

The C port reproduced the existing 1000-digit explicit recovery test byte-for-byte against the Python reference implementation.

Representative measurements from one cloud environment showed a substantial reduction in runtime and memory after the C17 port. These measurements are implementation diagnostics, not portable performance guarantees.

## What the precision claim means

An implicit precision claim is **not** the same thing as a conventional record consisting of a stored decimal expansion of \(\pi\).

The project studies how far the convergence state can be propagated and certified while leaving the enormous decimal representation unmaterialized. If an explicit expansion is requested, the output itself still has unavoidable large-scale arithmetic, conversion, storage, and I/O costs.

Accordingly, this repository should not describe an implicit certificate as a conventional decimal-digit world record unless the corresponding explicit result has actually been generated and independently accepted under the relevant record rules.

## Research status

This is experimental numerical/mathematical research software. The current certifier combines exact rational checks, generated modular relations, constructive recovery data, and numerical diagnostics. Independent mathematical review is welcome, especially for the modular-transform identities, branch selection, error bounds, and certificate semantics.

## License

A license file should be added before treating the repository as a finalized public release.
