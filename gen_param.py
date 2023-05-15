#!/usr/bin/python

import sys

from math import log
from math import ceil

import params

print("#ifndef PARAMS_H")
print("#define PARAMS_H")

for param in params.params:
  print(f"#ifdef {param.name}")

  print(f'  #define MEDS_name "{param.name}"')

  print()

  print(f"  #define MEDS_sec_bytes {param.sec}")
  print(f"  #define MEDS_salt_sec_bytes 32")
  print(f"  #define MEDS_p {param.q}")
  print(f"  #define GFq_t uint{ceil(log(param.q, 256))*8}_t")
  print(f"  #define GFq_bits {ceil(log(param.q, 2))}")

  print()

  print(f"  #define MEDS_m {param.m}")
  print(f"  #define MEDS_n {param.n}")
  print(f"  #define MEDS_k {param.k}")

  print()

  print(f"  #define MEDS_s {param.s}")
  print(f"  #define MEDS_t {param.t}")
  print(f"  #define MEDS_w {param.w}")

  if param.have_seed_tree:
    print()
    print(f"  #define HAVE_SEED_TREE")

  print()
   
  print(f"  #define MEDS_PK_BYTES {param.pk_size}")
  print(f"  #define MEDS_SIG_BYTES {param.sig_size}")

  print(f"#endif")
  print()

print("#endif")
print()

