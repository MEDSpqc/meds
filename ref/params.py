#!/usr/bin/env python3

from dataclasses import dataclass
from math import log
from math import log2
from math import ceil
from math import comb
import re

@dataclass
class Param:
  digest_bytes: int
  pub_seed_bytes: int
  sec_seed_bytes: int
  st_seed_bytes: int
  st_salt_bytes: int
  q: int
  n: int
  m: int
  s: int
  t: int
  w: int
  k: int
  seed_tree: bool = True

  _name: str=None

  def __str__(self):
    if self._name:
      return self._name
    else:
      return f"MEDS-{self.pk_size}"

  def __repr__(self):
    return self.name

  @property
  def name(self):
    if self._name:
      return re.sub(r'\W+', '', self._name)
    else:
      return f"MEDS{self.pk_size}"

  @property
  def pk_size(self):
    pub_seed_bytes = self.pub_seed_bytes
    s = self.s
    k = self.k
    m = self.m
    n = self.n
    q = self.q

    q_bits = ceil(log2(q))

    return (s - 1) * ceil(((k-2) * (m*n - k) * ceil(log2(q))) / 8) + pub_seed_bytes

  @property
  def sk_size(self):
    sec_seed_bytes = self.sec_seed_bytes
    s = self.s
    k = self.k
    m = self.m
    n = self.n
    q = self.q

    q_bits = ceil(log2(q))

    def mat_bytes(i, j):
       return ceil(i*j*q_bits/8)

    #return (s-1)*(mat_bytes(m,m) + mat_bytes(n,n)) + sec_seed_bytes + self.pub_seed_bytes
    return (s-1)*(mat_bytes(m,m) + mat_bytes(n,n) + mat_bytes(k,k)) + sec_seed_bytes + self.pub_seed_bytes

  @property
  def fiat_shamir(self):
    return log2(comb(self.t, self.w) * (self.s - 1)**self.w)

  @property
  def seed_tree_cost(self):
    return self.seed_max_tree_len * self.st_seed_bytes

  @property
  def seed_max_tree_len(self):
    st_seed_bytes = self.st_seed_bytes
    t = self.t
    w = self.w

    if self.seed_tree:
      return (2**ceil(log2(w)) + w * (ceil(log2(t)) - ceil(log2(w)) - 1))
    else:
      return t - w


  @property
  def sig_size(self):
    digest_bytes = self.digest_bytes
    m = self.m
    n = self.n
    k = self.k
    q = self.q
    t = self.t
    w = self.w

    q_bits = ceil(log2(q))

    sig_size = digest_bytes + w*ceil(2*k*q_bits/8) + self.seed_tree_cost
    #sig_size = digest_bytes + w*(ceil(m*m*q_bits/8) + ceil(n*n*q_bits/8)) + self.seed_tree_cost

    sig_size += self.st_salt_bytes

    return sig_size

  def dump(self):
    print()

    print(f"digest_bytes   = {self.digest_bytes}")
    print(f"pub_seed_bytes = {self.pub_seed_bytes}")
    print(f"sec_seed_bytes = {self.sec_seed_bytes}")
    print(f"st_seed_bytes  = {self.st_seed_bytes}")
    print(f"st_salt_bytes  = {self.st_salt_bytes}")

    print()

    print(f"q = {self.q:>3}")
    print(f"n = {self.n:>3}")
    print(f"m = {self.m:>3}")
    print(f"s = {self.s:>3}")
    print(f"t = {self.t:>3}")
    print(f"w = {self.w:>3}")
    print(f"k = {self.k:>3}")

    print()

    print(f"pk size:      {self.pk_size:7}")
    print(f"max sig size: {self.sig_size:7}")
    print(f"sum:          {self.pk_size + self.sig_size:7}")

    print()

    print(f"FS security:   {self.fiat_shamir:12.5f}")

    print()

params = [
# Test Level 1
  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 128 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 4093,
        m = 14,
        n = 15,
        k = 14,
        s = 2,
        t = 144,
        w = 48,
        seed_tree = False),

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 128 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 4093,
        m = 14,
        n = 15,
        k = 14,
        s = 4,
        t = 128,
        w = 26),
#        _name = "MEDS-2084"),

# # Level I
#   Param(digest_bytes = 256 >> 3,
#         pub_seed_bytes = 256 >> 3,
#         sec_seed_bytes = 256 >> 3,
#         st_seed_bytes = 128 >> 3,
#         st_salt_bytes = 256 >> 3,
#         q = 4093,
#         m = 14,
#         n = 15,
#         k = 14,
#         s = 4,
#         t = 1152,
#         w = 14),
# 
#   Param(digest_bytes = 256 >> 3,
#         pub_seed_bytes = 256 >> 3,
#         sec_seed_bytes = 256 >> 3,
#         st_seed_bytes = 128 >> 3,
#         st_salt_bytes = 256 >> 3,
#         q = 4093,
#         m = 14,
#         n = 15,
#         k = 14,
#         s = 5,
#         t = 192,
#         w = 20),

# Level III

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 192 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 4093,
        m = 22,
        n = 23,
        k = 22,
        s = 2,
        t = 224,
        w = 66,
        seed_tree = False),

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 192 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 4093,
        m = 22,
        n = 23,
        k = 22,
        s = 4,
        t = 192,
        w = 38,
        seed_tree = False),

#   Param(digest_bytes = 256 >> 3,
#         pub_seed_bytes = 256 >> 3,
#         sec_seed_bytes = 256 >> 3,
#         st_seed_bytes = 192 >> 3,
#         st_salt_bytes = 256 >> 3,
#         q = 4093,
#         m = 22,
#         n = 23,
#         k = 22,
#         s = 4,
#         t = 608,
#         w = 26),
# 
#   Param(digest_bytes = 256 >> 3,
#         pub_seed_bytes = 256 >> 3,
#         sec_seed_bytes = 256 >> 3,
#         st_seed_bytes = 192 >> 3,
#         st_salt_bytes = 256 >> 3,
#         q = 4093,
#         m = 22,
#         n = 23,
#         k = 22,
#         s = 5,
#         t = 160,
#         w = 36),

# Level V

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 256 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 2039,
        m = 30,
        n = 31,
        k = 30,
        s = 2,
        t = 288,
        w = 94,
        seed_tree = False),

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 256 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 2039,
        m = 30,
        n = 31,
        k = 30,
        s = 4,
        t = 256,
        w = 50),

#   Param(digest_bytes = 256 >> 3,
#         pub_seed_bytes = 256 >> 3,
#         sec_seed_bytes = 256 >> 3,
#         st_seed_bytes = 256 >> 3,
#         st_salt_bytes = 256 >> 3,
#         q = 2039,
#         m = 30,
#         n = 31,
#         k = 30,
#         s = 5,
#         t = 192,
#         w = 52),
# 
#   Param(digest_bytes = 256 >> 3,
#         pub_seed_bytes = 256 >> 3,
#         sec_seed_bytes = 256 >> 3,
#         st_seed_bytes = 256 >> 3,
#         st_salt_bytes = 256 >> 3,
#         q = 2039,
#         m = 30,
#         n = 31,
#         k = 30,
#         s = 6,
#         t = 112,
#         w = 66),

# toy

  Param(digest_bytes = 128 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 128 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 13, #4093,
        m = 4, #10,
        n = 5, #11,
        k = 4, #10,
        s = 2, #5,
        t = 4, #16,
        w = 2,
        _name = "toy"),
]

for i, p in enumerate(params):
  exec(f"{p.name} = params[i]")

par_list = {p.name : p for p in params}


def print_table():
  from tabulate import tabulate

  tab = []

  for param in params:
    fs = -param.fiat_shamir

    pk =  param.pk_size
    sk =  param.sk_size
    sig = param.sig_size

    tab.append([param, param.q, param.n, param.m, param.k, param.s, param.t, param.w, pk, sk, sig, fs])

  headers = ["set", "q", "n", "m", "k", "s", "t", "w", "pk", "sk", "sig", "fs"]

  tab = tabulate(tab, headers, tablefmt="latex_booktabs")

  print(tab)

def print_list():
  for param in params:
    print(param.name)

def interactive():
  param = params[0]

  def inject_params(param):
    global digest_bytes, pub_seed_bytes, sec_seed_bytes, st_seed_bytes, st_salt_bytes, q, n, m, s, t, w, k

    digest_bytes = param.digest_bytes
    pub_seed_bytes = param.pub_seed_bytes
    sec_seed_bytes = param.sec_seed_bytes
    st_seed_bytes = param.st_seed_bytes
    st_salt_bytes = param.st_salt_bytes
    q = param.q
    n = param.n
    m = param.m
    s = param.s
    t = param.t
    w = param.w
    k = param.k

  inject_params(param)

  def opt_tw(s, range_t=range(1, 1024)):
    @dataclass
    class TMP:
      sig_size: int

    tmp = TMP(100000000000000)

    for t in range_t:
      for w in range(1, min(t, 100)):
        loc = Param(digest_bytes, pub_seed_bytes, sec_seed_bytes, st_seed_bytes, st_salt_bytes, q, n, m, s, t, w, k)

        if loc.fiat_shamir > st_seed_bytes*8:
          if loc.sig_size < tmp.sig_size:
            tmp = loc

    print(f"sig size: {tmp.sig_size}  pk size: {tmp.pk_size} -> t = {tmp.t}, w = {tmp.w}")

  def dump():
    Param(digest_bytes, pub_seed_bytes, sec_seed_bytes, st_seed_bytes, st_salt_bytes, q, n, m, s, t, w, k).dump()

  import code
  import readline
  import rlcompleter

  vars = globals()
  vars.update(locals())

  readline.set_completer(rlcompleter.Completer(vars).complete)
  readline.parse_and_bind("tab: complete")
  code.InteractiveConsole(vars).interact(banner="""

Explore paramter space interactively.

Local variables q, n, m, ... can be dumped via 'dump()'.

(Exit python prompt with Ctrl-D)
""")


def gen_api(par_set):
  if not par_set:
    par_set = "toy"

  par_set = par_list[par_set]

  print(f"""#ifndef API_H
#define API_H

#define CRYPTO_SECRETKEYBYTES {par_set.sk_size}
#define CRYPTO_PUBLICKEYBYTES {par_set.pk_size}
#define CRYPTO_BYTES {par_set.sig_size}

#define CRYPTO_ALGNAME "{par_set.name}"

int crypto_sign_keypair(
    unsigned char *pk,
    unsigned char *sk
  );

int crypto_sign(
    unsigned char *sm, unsigned long long *smlen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *sk
  );

int crypto_sign_open(
    unsigned char *m, unsigned long long *mlen,
    const unsigned char *sm, unsigned long long smlen,
    const unsigned char *pk
  );

#endif
""")

def gen_param(parset):
  print("#ifndef PARAMS_H")
  print("#define PARAMS_H")
  print()

  if not parset:
    plist = params
  else:
    plist = [par_list[parset]]

  for param in plist:
    if not parset:
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
    print(f"{ind}#define GFq_fmt \"%{ceil(log(param.q, 10))}u\"")

    print()

    print(f"{ind}#define MEDS_m {param.m}")
    print(f"{ind}#define MEDS_n {param.n}")
    print(f"{ind}#define MEDS_k {param.k}")

    print()

    print(f"{ind}#define MEDS_s {param.s}")
    print(f"{ind}#define MEDS_t {param.t}")
    print(f"{ind}#define MEDS_w {param.w}")

    print()

    if param.seed_tree:
      print(f"{ind}#define MEDS_seed_tree_height {ceil(log(param.t, 2))}")
      print(f"{ind}#define SEED_TREE_size {((1 << (ceil(log(param.t, 2)) + 1)) - 1)}")
      print(f"{ind}#define MEDS_max_path_len {param.seed_max_tree_len}")
    else:
      print(f"{ind}#define MEDS_no_seed_tree")
      print(f"{ind}#define MEDS_seed_tree_height 0")
      print(f"{ind}#define SEED_TREE_size {param.t}")
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

    if not parset:
      print(f"#endif")
    print()

  print("#endif")
  print()



if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument("-t", "--table", action='store_true', help = "print Latex table")
  parser.add_argument("-l", "--list", action='store_true', help = "list param set names")
  parser.add_argument("-i", "--interactive", action='store_true', help = "interactive python console")
  parser.add_argument("-a", "--api", action='store_true', help = "generate api.h")
  parser.add_argument("-p", "--param", action='store_true', help = "generate param.h")
  parser.add_argument('parset', nargs='?', help = "parameter set", default=None)

  args = parser.parse_args()

  if args.table:
    print_table()
  elif args.list:
    print_list()
  elif args.interactive:
    interactive()
  elif args.api:
    gen_api(args.parset)
  elif args.param:
    gen_param(args.parset)
  else:
    print_table()

