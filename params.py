#!/usr/bin/python

from dataclasses import dataclass
from math import log2
from math import ceil
from math import comb

@dataclass
class Param:
  sec: int
  q: int
  n: int
  m: int
  s: int
  t: int
  w: int
  k: int
  have_seed_tree: bool

  _name: str=None

  def __str__(self):
    if self._name:
      return self._name
    else:
      return f"MEDS-{self.pk_size}{'-st' if self.have_seed_tree else ''}"

  def __repr__(self):
    return self.name

  @property
  def salt_sec(self):
    return 32

  @property
  def name(self):
    if self._name:
      return self._name
    else:
      return f"MEDS{self.pk_size}{'st' if self.have_seed_tree else ''}"

  @property
  def pk_size(self):
    sec = self.sec
    s = self.s
    k = self.k
    m = self.m
    n = self.n
    q = self.q
  
    return ceil(((s - 1) * (k * (m*n - k) - (m*n - k + (m-1)*n-k)) * ceil(log2(q))) / 8) + sec

  @property
  def sk_size(self):
    sec = self.sec
    s = self.s
    k = self.k
    m = self.m
    n = self.n
    q = self.q
  
    return ((s-1)*m*m + (s-1)*n*n + k*(m*n-k))*ceil(log2(q)/8) + sec

  
  @property
  def seed_tree_cost(self):
    sec = self.sec
    t = self.t
    w = self.w
  
    if self.have_seed_tree:
      return (2**ceil(log2(w)) + w * (ceil(log2(t)) - ceil(log2(w)) - 1)) * sec
    else:
      return (t-w) * sec
  
  @property
  def sig_size(self):
    sec = self.sec
    m = self.m
    n = self.n
    q = self.q
    t = self.t
    w = self.w

    sig_size = sec + w*ceil((m**2 + n**2)*ceil(log2(q))/8) + self.seed_tree_cost

    #if self.have_seed_tree:
    sig_size += self.salt_sec

    return sig_size


params = [
  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 2,
        t = 256,
        w = 30,
        have_seed_tree = True),
                
  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 4,
        t = 160,
        w = 23,
        have_seed_tree = True,
        _name = "MEDS8445stf"),

  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 4,
        t = 464,
        w = 17,
        have_seed_tree = True),

#  Param(sec = 128 >> 3,
#        q = 8191,
#        m = 13,
#        n = 13,
#        k = 13,
#        s = 4,
#        t = 1600,
#        w = 14,
#        have_seed_tree = True,
#        _name = "MEDS8445sts"),

  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 4,
        t = 1760,
        w = 13,
        have_seed_tree = True,
        _name = "MEDS8445sts"),

#  Param(sec = 128 >> 3,
#        q = 8191,
#        m = 13,
#        n = 13,
#        k = 13,
#        s = 5,
#        t = 448,
#        w = 16,
#        have_seed_tree = True),
 
  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 5,
        t = 224,
        w = 19,
        have_seed_tree = True),

  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 5,
        t = 224,
        w = 19,
        have_seed_tree = False),
 
  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 16,
        t = 128,
        w = 16,
        have_seed_tree = True),
 
  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 128,
        t = 80,
        w = 12,
        have_seed_tree = True),

  Param(sec = 128 >> 3,
        q = 8191,
        m = 13,
        n = 13,
        k = 13,
        s = 256,
        t = 64,
        w = 11,
        have_seed_tree = True),

  Param(sec = 128 >> 3,
        q = 8191,
        m = 10,
        n = 10,
        k = 10,
        s = 4,
        t = 16,
        w = 6,
        have_seed_tree = True,
        _name = "toy"),
]

for i, p in enumerate(params):
  exec(f"{p.name} = params[i]")

 

def print_table():
  from tabulate import tabulate
  
  tab = []

  for param in params:
    q = param.q
    n = param.n
    m = param.m
    s = param.s
    t = param.t
    w = param.w
    k = param.k
      
    fs = -log2(comb(t, w) * (s - 1)**w)
    #ac = cost.attack_cost(m, n, k, q)
    #wk = cost.weak_key(k, n, q)

    pk =  param.pk_size
    sig = param.sig_size

    #tab.append([param.name, param.q, param.n, param.k, param.s, param.t, param.w, 'X' if param.have_seed_tree else '--', pk, sig, ac, fs, wk])
    tab.append([param, param.q, param.n, param.m, param.k, param.s, param.t, param.w, 'checkmark' if param.have_seed_tree else '--', pk, sig, fs])

  #headers = ["set", "q", "n", "k", "s", "t", "w", "st", "pk", "sig", "ac", "fs", "wk"]
  headers = ["set", "q", "n", "m", "k", "s", "t", "w", "st", "pk", "sig", "fs"]

  tab = tabulate(tab, headers, tablefmt="latex_booktabs").replace('checkmark', '\\checkmark')

  print(tab)

def print_list():
  for param in params:
    print(param.name)


if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument("-t", "--table", action='store_true', help = "print Latex table")
  parser.add_argument("-l", "--list", action='store_true', help = "list param set names")

  args = parser.parse_args()

  if args.table:
    print_table()
  elif args.list:
    print_list()
  else:
    print_table()

