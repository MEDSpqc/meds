#!/usr/bin/python

from dataclasses import dataclass
from math import log2
from math import ceil
from math import comb

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
      return self._name
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

    return (s - 1) * ceil(((k * (m*n - k) - (m*n - k + (m-1)*n-k)) * ceil(log2(q))) / 8) + pub_seed_bytes

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

    return (s-1)*(mat_bytes(m,m) + mat_bytes(n,n)) + sec_seed_bytes + self.pub_seed_bytes

    #return ceil((((s-1)*m*m + (s-1)*n*n + k*(m*n-k))*ceil(log2(q)))/8) + sec_seed_bytes

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
  
    return (2**ceil(log2(w)) + w * (ceil(log2(t)) - ceil(log2(w)) - 1))
  
  @property
  def sig_size(self):
    digest_bytes = self.digest_bytes
    m = self.m
    n = self.n
    q = self.q
    t = self.t
    w = self.w

    q_bits = ceil(log2(q))

    sig_size = digest_bytes + w*(ceil(m*m*q_bits/8) + ceil(n*n*q_bits/8)) + self.seed_tree_cost

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
# Level I
  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 128 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 4093,
        m = 14,
        n = 14,
        k = 14,
        s = 4,
        t = 1152,
        w = 14),

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 128 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 4093,
        m = 14,
        n = 14,
        k = 14,
        s = 5,
        t = 192,
        w = 20),

# Level III

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 192 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 4093,
        m = 22,
        n = 22,
        k = 22,
        s = 4,
        t = 608,
        w = 26),

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 192 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 4093,
        m = 22,
        n = 22,
        k = 22,
        s = 5,
        t = 160,
        w = 36),

# Level V

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 256 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 16381,
        m = 25,
        n = 25,
        k = 25,
        s = 4,
        t = 768,
        w = 35),

  Param(digest_bytes = 256 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 256 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 16381,
        m = 25,
        n = 25,
        k = 25,
        s = 5,
        t = 304,
        w = 42),

# toy

  Param(digest_bytes = 128 >> 3,
        pub_seed_bytes = 256 >> 3,
        sec_seed_bytes = 256 >> 3,
        st_seed_bytes = 128 >> 3,
        st_salt_bytes = 256 >> 3,
        q = 8191,
        m = 10,
        n = 10,
        k = 10,
        s = 4,
        t = 16,
        w = 6,
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
    sig = param.sig_size

    tab.append([param, param.q, param.n, param.m, param.k, param.s, param.t, param.w, pk, sig, fs])

  headers = ["set", "q", "n", "m", "k", "s", "t", "w", "pk", "sig", "fs"]

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

if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument("-t", "--table", action='store_true', help = "print Latex table")
  parser.add_argument("-l", "--list", action='store_true', help = "list param set names")
  parser.add_argument("-i", "--interactive", action='store_true', help = "interactive python console")

  args = parser.parse_args()

  if args.table:
    print_table()
  elif args.list:
    print_list()
  elif args.interactive:
    interactive()
  else:
    print_table()

