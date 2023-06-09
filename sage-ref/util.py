from sage.all_cmdline import *   # import sage libraryS

from Crypto.Hash import SHAKE256
from Crypto.Hash import SHA256

import logging

from SeedTree import SeedTree

def compress(m, bs):
  GFq = m[0,0].base_ring()
  GF_BITS = ceil(log(GFq.base().order(),2))
  if GFq.is_prime_field():
    for v in [j for v in m.rows() for j in v]:
      bs.write(int(v), GF_BITS)
  else:
    raise "Not implemented..."

cache = {}

def rnd_GFq(seed, GFq, min=0):
  global cache

  if GFq in cache:
    numbytes, mask, q = cache[GFq]
  else:
    numbytes = ceil(log(GFq.base().order(),2)/8)
    mask = (1 << ceil(log(GFq.base().order(),2))) - 1
    q = GFq.base().order()

    cache[GFq] = numbytes, mask, q 

  if type(seed) == SHAKE256.SHAKE256_XOF:
    shake = seed
  else:
    shake = SHAKE256.new()
    shake.update(seed)

  while True:
    val = 0

    for i in range(numbytes):
      val += ord(shake.read(1)) << (i*8)

    val = val & mask

    if val >= min and val < q:
      return GFq(val)

def rnd_matrix(seed, GFq, m, n=None):
  if not n:
    n = m

  shake = SHAKE256.new()
  shake.update(seed)

  return matrix(GFq, m, n, [rnd_GFq(shake, GFq) for i in range(m*n)])

def rnd_inv_matrix(seed, GFq, n):
  while True:
    M = rnd_matrix(seed, GFq, n, n)

    if M.is_invertible():
      return M

def rnd_sys_matrix(seed, GFq, k, m, n):
  I_k = matrix.identity(ring=GFq, n=k)

  return I_k.augment(rnd_matrix(seed, GFq, k, m*n-k))


def pi(A, B, G):
  m = A.nrows()
  n = B.nrows()
  k = G.nrows()

  GFq = A[0,0].base_ring()

  G = [matrix(GFq, m, n, row) for row in G.rows()]

  AGB = [A*v*B for v in G]

  return matrix(GFq, k, m*n, [AGB[i][j,g] for i in range(k) for j in range(m) for g in range(n)])


def XOF(seed, length):
  shake = SHAKE256.new()
  shake.update(seed)

  if not isinstance(length, list):
    return shake.read(length)

  data = shake.read(sum(length))

  return tuple(data[sum(length[:i]):sum(length[:i+1])] for i in range(len(length)))


def H(params):
  def hash(m):
    shake = SHAKE256.new()
    shake.update(m)
  
    return shake.read(params.digest_bytes)

  return hash

def G(params):
  def hash_pair(value):
    shake = SHAKE256.new()
    shake.update(value)

    return shake.read(params.st_seed_bytes), shake.read(params.st_seed_bytes)

  return hash_pair


def seeds_from_path(h, path, salt, param):
  # prepare tree structure
  seeds = SeedTree(param.t, salt = salt)

  for i, v in enumerate(h):
    if v > 0:
      seeds[i] = object

  parsed = []

  for i, v in enumerate(seeds.path()):
    if v == object:
      parsed.append(None)
    else:
      parsed.append(path[:param.st_seed_bytes])

      path = path[param.st_seed_bytes:]

  # apply patch once tree structure is set up
  seeds.patch(parsed, G(param))

  return seeds

def parse_hash(digest, params):
  t = params.t
  s = params.s
  w = params.w

  logging.debug(f"digest: %s", [int(v) for v in digest])
  logging.debug("digest len: %i", len(digest))

  shake = SHAKE256.new()

  shake.update(digest)

  h = [0] * t

  num = 0

  while num < w:
    pos = 0

    for i in range(ceil(log(t,2)/8)):
      pos += ord(shake.read(1)) << (i*8)

    pos &= (1 << ceil(log(t,2))) - 1


    if pos >= t:
      continue

    if h[pos] > 0:
      continue  

    logging.debug(f"pos: {pos}")

    val = 0

    while val < 1 or val > s-1:
      val = 0

      for i in range(ceil(log(s,2)/8)):
        val += ord(shake.read(1)) << (i*8)
    
      val &= (1 << ceil(log(s,2))) - 1

    h[pos] = val

    logging.debug(f"p: {pos}  v: {val}")

    num += 1

  return h


def solve_opt(P0prime, Amm):
  m = P0prime[0].nrows()
  n = P0prime[0].ncols()

  GFq = Amm.base_ring()

  N = -P0prime[0].transpose()

  #logging.debug(f"N:\n%s", N)


  M = matrix(GFq, n, m + m + 2)

  for i in range(m):
    for j in range(n):
      M[j, i] = -P0prime[1][i, j]

  for i in range(m):
    for j in range(n):
      M[j, i+n] = P0prime[0][i, j]

  for j in range(n):
    M[j, m+n] = - P0prime[0][m-1,j]*Amm

  for j in range(n):
    M[j, m+n+1] = P0prime[1][m-1, j]*Amm

  #logging.debug(f"M:\n%s", M)

  Ms = M.matrix_from_rows_and_columns(range(M.nrows()-1), range(M.ncols()))

  Ms = Ms.rref()
  M = Ms.stack(M.rows()[-1])

  if not all([Ms[i,i] == 1 for i in range(n-1)]):
    return None, None

  #logging.debug(f"M part:\n%s", M)

  for i in range(n-1):
    M.add_multiple_of_row(n-1, i, -M[n-1,i])

  if M[n-1,n-1] == 0:
    return None, None

  M.set_row_to_multiple_of_row(n-1, n-1, 1/M[n-1,n-1])

  M[-1, -1] = GFq(0)

  #logging.debug(f"M red:\n%s", M) 

  for i in range(n-1):
    M.add_multiple_of_row(i, n-1, -M[i,n-1])

  #logging.debug(f"M done:\n%s", M) 


  sol = [0] * (n*n + m*m)

  sol[-1] = Amm

  for i in range(n-1):
    sol[n*n + m*m - n + i] = M[i, -1]

  for i in range(n):
    sol[n*n + m*m - 2*n + i] = M[i, -2]

  for i in range(n):
    sol[n*n - n + i] = P0prime[0][m-1,i]*Amm

  #logging.debug(f"initial sol:\n%s", sol)

  # incomplete blocks:

  for i in range(n):
    for j in range(n-1):
      sol[n*n + m*m - 2*n + i] += - M[i, n+n-2-j] * sol[n*n + m*m - 2 - j]
      sol[n*n - n + i] += - N[i, m-2-j] * sol[n*n + m*m - 2 - j]

  #logging.debug(f"incomplete blocks:\n%s", sol)


  # complete blocks:

  for block in range(3,n+1):
    for i in range(n):
      for j in range(n):
        sol[n*n + m*m - block*n + i] += - M[i, n+n-1-j] * sol[n*n + m*m - 1 - (block-2)*n - j]

  for block in range(2,n+1):
    for i in range(n):
      for j in range(n):
        sol[n*n - block*n + i] += - N[i, m-1-j] * sol[n*n + m*m - 1 - (block-1)*n - j]

  #logging.debug(f"complete blocks:\n%s", sol)


  A = matrix(GFq, m,m, sol[n*n:])
  B_inv = matrix(GFq, n,n, sol[:n*n])

  #logging.debug(f"A:\n%s", A)
  #logging.debug(f"B_inv:\n%s", B_inv)

  return A, B_inv

def solve_symb(P0prime, Amm):
  m = P0prime[0].nrows()
  n = P0prime[0].ncols()

  GFq = Amm.base_ring()

  Pj = [None] * 2

  Pj[0] = matrix(GFq, m, n, [[GFq(1) if i==j else GFq(0) for i in range(n)] for j in range(m)])
  Pj[1] = matrix(GFq, m, n, [[GFq(1) if i==j else GFq(0) for i in range(n)] for j in range(1,m)] + [[GFq(0)]*n])

  R = PolynomialRing(GFq, m*m + n*n,
     names = ','.join([f"b{i}_{j}" for i in range(n) for j in range(n)]) + "," \
           + ','.join([f"a{i}_{j}" for i in range(m) for j in range(m)]))

  A     = matrix(R, m, var(','.join([f"a{i}_{j}" for i in range(m) for j in range(m)])))
  B_inv = matrix(R, n, var(','.join([f"b{i}_{j}" for i in range(n) for j in range(n)])))

  A[m-1,m-1] = Amm

  eqs1 = Pj[0] * B_inv - A*P0prime[0]
  eqs2 = Pj[1] * B_inv - A*P0prime[1]

  eqs = eqs1.coefficients() + eqs2.coefficients()[:-1]

  rsys = matrix(GFq, [[eq.coefficient(v) for v in R.gens()[:-1]] + [-eq.constant_coefficient()] for eq in eqs])

  rsys_rref = rsys.rref()

  if not all([rsys_rref[i][i] == 1 for i in range(rsys_rref.nrows())]):
    #logging.debug("no sol")
    return None, None

  sol = rsys_rref.columns()[-1].list()

  A = matrix(GFq, m, sol[n*n:] + [Amm])
  B_inv = matrix(GFq, m, sol[:n*n])

  #logging.debug(f"A:\n%s", A)
  #logging.debug(f"B_inv:\n%s", B_inv)

  return A, B_inv

def solve_comp(P0prime, Amm):
  m = P0prime[0].nrows()
  n = P0prime[0].ncols()

  GFq = Amm.base_ring()

  Pj = [None] * 2

  Pj[0] = matrix(GFq, m, n, [[GFq(1) if i==j else GFq(0) for i in range(n)] for j in range(m)])
  Pj[1] = matrix(GFq, m, n, [[GFq(1) if i==j else GFq(0) for i in range(n)] for j in range(1,m)] + [[GFq(0)]*n])

  rsys = []

  for l in range(n):
    for j in range(m):
      tmp = [GFq(0)] * (m^2 + n^2)

      for i in range(m):
        tmp[m*m + l*m + i] = - P0prime[0][i,j]

      for i in range(n):
        tmp[i*n + j] = Pj[0][l,i]

      if l == n-1:
        tmp[m*m + n*n - 1] = Amm * P0prime[0][n-1,j]

      rsys.append(tmp)

  for l in range(n):
    for j in range(m if l < n-1 else m-1):
      tmp = [GFq(0)] * (m^2 + n^2)

      for i in range(m):
        tmp[m*m + l*m + i] = - P0prime[1][i,j]

      for i in range(n):
        tmp[i*n + j] = Pj[1][l,i]

      if l == n-1:
        tmp[m*m + n*n - 1] = Amm * P0prime[1][n-1,j]

      rsys.append(tmp)

  rsys = matrix(GFq, ncols=m^2 + n^2, entries=rsys)

  rsys_rref = rsys.rref()

  if not all([rsys_rref[i][i] == 1 for i in range(rsys_rref.nrows())]):
    #logging.debug("no sol")
    return None, None

  sol = rsys_rref.columns()[-1].list()

  A = matrix(GFq, m, sol[n*n:] + [Amm])
  B_inv = matrix(GFq, m, sol[:n*n])

  #logging.debug(f"A:\n%s", A)
  #logging.debug(f"B_inv:\n%s", B_inv)

  return A, B_inv


