from sage.all_cmdline import *   # import sage librarys

from Crypto.Hash import SHAKE256
from Crypto.Hash import SHA256

import logging

from seedtree import SeedTree
from bitstream import BitStream


GFq_cache = {}

def get_cached(GFq):
  global GFq_cache

  if GFq not in GFq_cache:
    numbytes = ceil(log(GFq.base().order(),2)/8)
    GF_BITS = ceil(log(GFq.base().order(),2))
    mask = (1 << ceil(log(GFq.base().order(),2))) - 1
    q = GFq.base().order()

    GFq_cache[GFq] = numbytes, mask, q, GF_BITS

  return GFq_cache[GFq]



def Compress(M):
  GFq = M[0,0].base_ring()
  _, _, _, GF_BITS = get_cached(GFq)

  bs = BitStream()

  for v in [j for v in M.rows() for j in v]:
    bs.write(int(v), GF_BITS)

  return bs.data

def Decompress(b, GFq, r, c):
  _, _, _, GF_BITS = get_cached(GFq)

  bs = BitStream(b)

  data = [GFq(bs.read(GF_BITS)) for _ in range(r*c)]

  return matrix(GFq, r, c, data)


def CompressG(M, k, m, n):
  GFq = M[0,0].base_ring()
  _, _, _, GF_BITS = get_cached(GFq)

  bs = BitStream()

  for i in range(n):
    bs.write(int(M[1,m*n-n+i]), GF_BITS)

  for i in range(2, k):
    for j in range(k, m*n):
      bs.write(int(M[i, j]), GF_BITS)

  return bs.data

def DecompressG(b, GFq, k, m, n):
  _, _, _, GF_BITS = get_cached(GFq)

  bs = BitStream(b)

  I_k = matrix.identity(ring=GFq, n=k)

  G = I_k.augment(matrix(GFq, k, m*n-k))

  for i in range(1, m):
    G[0, i*(n+1)] = 1

  for i in range(1, m-1):
    G[1, i*(n+1)+1] = 1

  for i in range(n):
    G[1,m*n-n+i] = GFq(bs.read(GF_BITS))

  for i in range(2, k):
    for j in range(k, m*n):
      G[i, j] = GFq(bs.read(GF_BITS))

  return G


def ExpandFqs(seed, num, GFq):
  numbytes, mask, q, _ = get_cached(GFq)

  if type(seed) == SHAKE256.SHAKE256_XOF:
    shake = seed
  else:
    shake = SHAKE256.new()
    shake.update(seed)

  ret = []

  for _ in range(num):
    while True:
      val = 0

      for i in range(numbytes):
        val += ord(shake.read(1)) << (i*8)

      val = val & mask

      if val < q:
        ret.append(GFq(val))
        break

  return ret

def ExpandInvMat(seed, GFq, n):
  shake = SHAKE256.new()
  shake.update(seed)

  while True:
    M = matrix(GFq, n, n, ExpandFqs(shake, n*n, GFq))

    if M.is_invertible():
      return M

def ExpandSystMat(seed, GFq, k, m, n):
  I_k = matrix.identity(ring=GFq, n=k)

  return I_k.augment(matrix(GFq, k, m*n-k, ExpandFqs(seed, k*(m*n-k), GFq)))


def pi(A, B, G):
  m = A.nrows()
  n = B.nrows()
  k = G.nrows()

  GFq = A[0,0].base_ring()

  G = [matrix(GFq, m, n, row) for row in G.rows()]

  AGB = [A*v*B for v in G]

  return matrix(GFq, k, m*n, [AGB[i][j,g] for i in range(k) for j in range(m) for g in range(n)])

def SF(M):
  M = M.echelon_form()

  # check if we got systematic form
  if sum([M[j,j] for j in range(M.nrows())]) != M.nrows():
    return None

  return M


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

def G(params, salt):
  def hash_pair(value, addr):
    shake = SHAKE256.new()
    shake.update(salt + value + addr.to_bytes(4, "little"))

    return shake.read(params.st_seed_bytes), shake.read(params.st_seed_bytes)

  return hash_pair

def SeedTreeToPath(h, root, salt, param):
  seeds = SeedTree(param.t, root, G(param, salt))

  for i, h_i in enumerate(h):
    if h_i > 0:
      seeds.delete(i)

  ret = b"".join(seeds.path())

  return ret + bytes([0]) * (param.seed_tree_cost - len(ret))

def PathToSeedTree(h, path, salt, param):
  # prepare tree structure
  seeds = SeedTree(param.t)

  for i, v in enumerate(h):
    if v > 0:
      seeds.delete(i)

  parsed = []

  for _ in range(len(seeds.path())):
    parsed.append(path[:param.st_seed_bytes])
    path = path[param.st_seed_bytes:]

  # apply patch once tree structure is set up
  seeds.patch(parsed, G(param, salt))

  return seeds

def PaseHash(digest, params):
  t = params.t
  s = params.s
  w = params.w

  #logging.debug(f"digest:\n%s", [int(v) for v in digest])
  #logging.debug("digest_len:\n%i", len(digest))

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

    #logging.debug(f"pos: {pos}")

    val = 0

    while val < 1 or val > s-1:
      val = 0

      for i in range(ceil(log(s,2)/8)):
        val += ord(shake.read(1)) << (i*8)

      val &= (1 << ceil(log(s,2))) - 1

    h[pos] = val

    #logging.debug(f"p: {pos}  v: {val}")

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

def solve_symb(P0prime):
  m = P0prime[0].nrows()
  n = P0prime[0].ncols()

  GFq = P0prime[0][0,0].base_ring()

  Pj = [None] * 2

  Pj[0] = matrix(GFq, m, n, [[GFq(1) if i==j   else GFq(0) for i in range(n)] for j in range(m)])
  Pj[1] = matrix(GFq, m, n, [[GFq(1) if i==j+1 else GFq(0) for i in range(n)] for j in range(m)])

  R = PolynomialRing(GFq, m*m + n*n,
     names = ','.join([f"a{i}_{j}" for i in range(m) for j in range(m)]) + "," \
           + ','.join([f"b{i}_{j}" for i in range(n) for j in range(n)]))

  A     = matrix(R, m, var(','.join([f"a{i}_{j}" for i in range(m) for j in range(m)])))
  B_inv = matrix(R, n, var(','.join([f"b{i}_{j}" for i in range(n) for j in range(n)])))

  eqs1 = Pj[0] * B_inv - A*P0prime[0]
  eqs2 = Pj[1] * B_inv - A*P0prime[1]

  eqs = eqs1.coefficients() + eqs2.coefficients()

  rsys = matrix(GFq, [[eq.coefficient(v) for v in R.gens()] for eq in eqs])

  logging.debug(f"rsys:\n%s", rsys)

  rsys_rref = rsys.rref()

  logging.debug(f"rsys:\n%s", rsys_rref)

  if not all([rsys_rref[i][i] == 1 for i in range(rsys_rref.nrows())]):
    #logging.debug("no sol")
    return None, None

  sol = rsys_rref.columns()[-1].list()

  logging.debug("sol:\n%s", sol)

  A = matrix(GFq, m, sol[:m*m])
  B_inv = matrix(GFq, n, sol[m*m:] + [GFq(-1)])

  logging.debug(f"A:\n%s", A)
  logging.debug(f"B_inv:\n%s", B_inv)

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


