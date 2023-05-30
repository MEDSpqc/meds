# MEDS
## Matrix Equivalence Digital Signature

This repository procides the reference implentation of the PQC signature scheme MEDS
as submitted to the NIST PQC Signature standardization process.

## Compile

The referemce implementation in directory `/ref` can be compiled
using the provided Makefile.

We provide three programs:
`test` to run at test of key generation, signing, and verification,
`bench` for benchmarling the implentation using several rund, and
`KAT_test` for computing known answert tests.

The test can be compiled and run by

`make RUN`

the benchmark using

`make BENCH`

and the KAT test using

`make KAT`.

The default parameter set is the `toy` parameter set. Another parameter set can be selected using `PARAM`, e.g.:

`make RUN PARAM=MEDS9923`

A list of available parameter sets can be obtained by:

`./params.py -l`

To run all targets, add `_ALL` to `RUN`, `BENCH`, and `KAT`, e.g.:

`make RUN_ALL`

When the code is compiled with `DEBUG` defined, exhaustive step-by-step debugging is produced and written to `stderr`.


## Africacrypt 2023

The source code accompanying the Africacrypt 2023 paper:

- Tung Chou, Ruben Niederhagen, Edoardo Persichetti,
  Tovohery Hajatiana Randrianarisoa, Krijn Reijnders, Simona Samardjiska,
  and Monika Trimoska:
  *"Take your MEDS: Digital Signatures from Matrix Code Equivalence"*.
  Progress in Cryptology - AfricaCrypt 2023.
  Lecture Notes in Computer Science, Springer, 2023.

is in branch `Africacrypt`.

