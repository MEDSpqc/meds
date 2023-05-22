#!/usr/bin/python

import sys
import argparse

from math import log
from math import ceil

sys.path.append("../sage-ref")
import params


parser = argparse.ArgumentParser()

parser.add_argument('parset', nargs='?', default=None, help = "parameter set")

args = parser.parse_args()


print("#ifndef PARAMS_H")
print("#define PARAMS_H")
print()

if not args.parset:
  par_list = params.params
else:
  par_list = [params.par_list[args.parset]]

for param in par_list:
  if not args.parset:
    print(f"#ifdef {param.name}")
    ind = "  "
  else:
   ind = ""

  print(f'{ind}#define MEDS_name "{param.name}"')

  print()

  print(f"{ind}#define MEDS_digest_bytes {param.digest_bytes}")
  print(f"{ind}#define MEDS_pub_seed_bytes {param.pub_seed_bytes}")
  print(f"{ind}#define MEDS_sec_seed_bytes {param.sec_seed_bytes}")
  print(f"{ind}#define MEDS_st_seed_bytes {param.st_seed_bytes}")

  print(f"{ind}#define MEDS_st_salt_bytes {param.st_salt_bytes}")

  print()

  print(f"{ind}#define MEDS_p {param.q}")
  print(f"{ind}#define GFq_t uint{ceil(log(param.q, 256))*8}_t")
  print(f"{ind}#define GFq_bits {ceil(log(param.q, 2))}")
  print(f"{ind}#define GFq_bytes {ceil(ceil(log(param.q, 2))/8)}")

  print()

  print(f"{ind}#define MEDS_m {param.m}")
  print(f"{ind}#define MEDS_n {param.n}")
  print(f"{ind}#define MEDS_k {param.k}")

  print()

  print(f"{ind}#define MEDS_s {param.s}")
  print(f"{ind}#define MEDS_t {param.t}")
  print(f"{ind}#define MEDS_w {param.w}")

  print()

  print(f"{ind}#define MEDS_seed_tree_height {ceil(log(param.t, 2))}")
  print(f"{ind}#define SEED_TREE_size {((1 << (ceil(log(param.t, 2)) + 1)) - 1)}")

  print(f"{ind}#define MEDS_max_path_len {param.seed_max_tree_len}")

  print()

  print(f"{ind}#define MEDS_t_mask 0x{2**ceil(log(param.t, 2)) - 1:08X}")
  print(f"{ind}#define MEDS_t_bytes {ceil(log(param.t-1, 2)/8)}")
  print()
  print(f"{ind}#define MEDS_s_mask 0x{2**ceil(log(param.s, 2)) - 1:08X}")

  print()
   
  print(f"{ind}#define MEDS_PK_BYTES {param.pk_size}")
  print(f"{ind}#define MEDS_SK_BYTES {param.sk_size}")
  print(f"{ind}#define MEDS_SIG_BYTES {param.sig_size}")

  if not args.parset:
    print(f"#endif")
  print()

print("#endif")
print()

