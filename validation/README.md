# Validation

This directory records reproducibility and implementation-validation artifacts for the active C17 mainline.

## Reproduce locally

```bash
make NATIVE=0
make NATIVE=0 smoke
```

The smoke test performs the following chain:

```text
C17 certifier
  -> freshly generated constructive certificate
  -> standalone C17 restorer
  -> 1000 explicit digits after the decimal point
  -> SHA-256 comparison
```

Expected SHA-256 for the generated `pi1000.txt`:

```text
e898fea26734a6d3af5396b9f4c60ae5dcc88fc40944d835911a9ee8a672ea1b
```

`PORT_VALIDATION.txt` records the C-port parity checks, representative cloud timings, and the boundary between current implementation evidence and theorem-level proof.

The validation artifacts do not redefine the mathematical algorithm. They are intended to make the implementation independently inspectable and repeatable.
