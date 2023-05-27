# MEDS

This is the MEDS implementation accompanying the Africacrypt paper:

  Tung Chou, Ruben Niederhagen, Edoardo Persichetti,
  Tovohery Hajatiana Randrianarisoa, Krijn Reijnders, Simona Samardjiska,
  and Monika Trimoska:
  *"Take your MEDS: Digital Signatures from Matrix Code Equivalence"*.
  Progress in Cryptology - AfricaCrypt 2023.
  Lecture Notes in Computer Science, Springer, 2023.

## How to compile the code

Call

  `make RUN`

to compile and run a toy-parameter set example.

You can get a list of all parameter sets from `params.py`:

  `./params.py -l`

Then compile and run the code by, e.g.:

  `make RUN PARAM=MEDS2826st`

Run all parameter sets by:

  `make RUN_ALL`

The same works for the benchmarking, e.g.:

  `make BENCH`

  `make BENCH PARAM=MEDS2826st`

  `make BENCH_ALL`

