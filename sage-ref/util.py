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

  for i in range(1, m):
    G[1, i*(n+1)+1] = 1

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

def SF(G):
  M = G.submatrix(0,0,G.nrows(),G.nrows())

  if M.is_invertible():
    return M.inverse() * G

  return None


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

#def solve_opt(P0prime, Amm):
#  m = P0prime[0].nrows()
#  n = P0prime[0].ncols()
#
#  GFq = Amm.base_ring()
#
#  N = -P0prime[0].transpose()
#
#  #logging.debug(f"N:\n%s", N)
#
#
#  M = matrix(GFq, n, m + m + 2)
#
#  for i in range(m):
#    for j in range(n):
#      M[j, i] = -P0prime[1][i, j]
#
#  for i in range(m):
#    for j in range(n):
#      M[j, i+n] = P0prime[0][i, j]
#
#  for j in range(n):
#    M[j, m+n] = - P0prime[0][m-1,j]*Amm
#
#  for j in range(n):
#    M[j, m+n+1] = P0prime[1][m-1, j]*Amm
#
#  #logging.debug(f"M:\n%s", M)
#
#  Ms = M.matrix_from_rows_and_columns(range(M.nrows()-1), range(M.ncols()))
#
#  Ms = Ms.rref()
#  M = Ms.stack(M.rows()[-1])
#
#  if not all([Ms[i,i] == 1 for i in range(n-1)]):
#    return None, None
#
#  #logging.debug(f"M part:\n%s", M)
#
#  for i in range(n-1):
#    M.add_multiple_of_row(n-1, i, -M[n-1,i])
#
#  if M[n-1,n-1] == 0:
#    return None, None
#
#  M.set_row_to_multiple_of_row(n-1, n-1, 1/M[n-1,n-1])
#
#  M[-1, -1] = GFq(0)
#
#  #logging.debug(f"M red:\n%s", M) 
#
#  for i in range(n-1):
#    M.add_multiple_of_row(i, n-1, -M[i,n-1])
#
#  #logging.debug(f"M done:\n%s", M) 
#
#
#  sol = [0] * (n*n + m*m)
#
#  sol[-1] = Amm
#
#  for i in range(n-1):
#    sol[n*n + m*m - n + i] = M[i, -1]
#
#  for i in range(n):
#    sol[n*n + m*m - 2*n + i] = M[i, -2]
#
#  for i in range(n):
#    sol[n*n - n + i] = P0prime[0][m-1,i]*Amm
#
#  #logging.debug(f"initial sol:\n%s", sol)
#
#  # incomplete blocks:
#
#  for i in range(n):
#    for j in range(n-1):
#      sol[n*n + m*m - 2*n + i] += - M[i, n+n-2-j] * sol[n*n + m*m - 2 - j]
#      sol[n*n - n + i] += - N[i, m-2-j] * sol[n*n + m*m - 2 - j]
#
#  #logging.debug(f"incomplete blocks:\n%s", sol)
#
#
#  # complete blocks:
#
#  for block in range(3,n+1):
#    for i in range(n):
#      for j in range(n):
#        sol[n*n + m*m - block*n + i] += - M[i, n+n-1-j] * sol[n*n + m*m - 1 - (block-2)*n - j]
#
#  for block in range(2,n+1):
#    for i in range(n):
#      for j in range(n):
#        sol[n*n - block*n + i] += - N[i, m-1-j] * sol[n*n + m*m - 1 - (block-1)*n - j]
#
#  #logging.debug(f"complete blocks:\n%s", sol)
#
#
#  A = matrix(GFq, m,m, sol[n*n:])
#  B_inv = matrix(GFq, n,n, sol[:n*n])
#
#  #logging.debug(f"A:\n%s", A)
#  #logging.debug(f"B_inv:\n%s", B_inv)
#
#  return A, B_inv

def solve(data):
  return solve_symb(data)
  ##return solve_comp(data)
  #return solve_opt(data)

def solve_symb(C):
  m = C[0].nrows()
  n = C[0].ncols()

  GFq = C[0][0,0].base_ring()


  D = [None] * 2

  D[0] = matrix(GFq, m, n, [[GFq(1) if i==j   else GFq(0) for i in range(n)] for j in range(m)])
  D[1] = matrix(GFq, m, n, [[GFq(1) if i==j+1 else GFq(0) for i in range(n)] for j in range(m)])

  R = PolynomialRing(GFq, m*m + n*n,
     names = ','.join([f"b{i}_{j}" for i in range(n) for j in range(n)]) + "," \
           + ','.join([f"a{i}_{j}" for i in range(m) for j in range(m)]))

  A_tilde     = matrix(R, m, var(','.join([f"a{i}_{j}" for i in range(m) for j in range(m)])))
  B_tilde_inv = matrix(R, n, var(','.join([f"b{i}_{j}" for i in range(n) for j in range(n)])))

  eqs1 = D[0] * B_tilde_inv - A_tilde * C[0]
  eqs2 = D[1] * B_tilde_inv - A_tilde * C[1]

  eqs1 = eqs1.coefficients()
  eqs2 = eqs2.coefficients()

  eqs = eqs1 + eqs2

  rsys = matrix(GFq, [[eq.coefficient(v) for v in R.gens()] for eq in eqs])

#  logging.debug(f"rsys:\n%s", rsys)

  rsys_rref = rsys.rref()

#  logging.debug(f"rsys:\n%s", rsys_rref)

  pivot_row = -1

  off = 0
 
  for i in range(rsys_rref.nrows()):
    if rsys_rref[i,i+off] != 1:
      if pivot_row == -1:
        pivot_row = i
        off = 1
        continue
      if pivot_row >= 0:
        logging.debug(f"no sol")
        return None, None

  if pivot_row < 0:
    pivot_row = rsys_rref.ncols() - 1

  sol = rsys_rref.columns()[pivot_row].list()[:pivot_row] + [GFq(-1)] + (m*m + n*n - pivot_row - 1) * [0]

  logging.debug("sol:\n%s", sol)

  A_tilde     = matrix(GFq, m, sol[n*n:])
  B_tilde_inv = matrix(GFq, n, sol[:n*n])

  logging.debug(f"A_tilde:\n%s", A_tilde)
  logging.debug(f"B_tilde_inv:\n%s", B_tilde_inv)

  return A_tilde, B_tilde_inv


def solve_opt(P0prime):
  m = P0prime[0].nrows()
  n = P0prime[0].ncols()

  GFq = P0prime[0][0,0].base_ring()

  ###################
  N = -P0prime[1].transpose()
  N = N.augment(P0prime[0].transpose())

  N0 = N.submatrix(0,0,N.nrows(), N.ncols())
  N1 = N.submatrix(N.nrows()-1,0,1, N.ncols())

  logging.debug(f"N:\n%s", N0)

  # compute systematic form of left half of the matrix
  _, _, U = N.LU()
  
  M = U[range(m),range(2*m)]
  
  N = M.rref().stack(U[[m],range(2*m)])

  for i in range(m):
    if N[i,i] != 1:
      #if (i == N.nrows()-1) and partial:
      #  continue
      #else:
        logging.debug(f"no sol")
        return None, None

  logging.debug(f"N:\n%s", N)


  tmp = N.submatrix(0, m, n, m)


  for i in range(m-1):
    N = N.stack(N1)


  N1 = N.submatrix(m, m, 1, N.ncols()-m)
  N1 = N1.stack(N.submatrix(0, m, m-2, N.ncols()-m))
  for i in range(1, m-1):
    for j in range(N1.ncols()):
      N1[i,j] = 0

  logging.debug(f"N1:\n%s", N1)

  N = N.submatrix(0, m, m, m)

  logging.debug(f"N:\n%s", N)

  for row in range(1, m-1):
    for i in range(0, m):
      for j in range(0, m):
        tmp0 = N1[row-1, i]
        tmp1 = N[i, j]

        prod = tmp0 * tmp1

        tmp2 = N1[row, j]

        diff = tmp2 - prod

        N1[row, j] = diff

#  for block in range(1, m):
#    for row in range(block):
#      for i in range(m):
#        for j in range(m):
#          N1[row, j+m] = N1[row, j+m] - N1[row, i] * N[i, j+m]
#
#      for i in range(m):
#        N1[row, i] = N1[row, i+m]
#        N1[row, i+m] = 0

  logging.debug(f"N1:\n%s", N1)

  N1 = N1.submatrix(0,0,m-1, m)

  N1 = N1.rref()

  logging.debug(f"N1:\n%s", N1)

  
#  logging.debug(f"tmp:\n%s", tmp)
#
#  for i in range(n):
#    for j in range(n):
#      if i != j:
#        tmp.add_multiple_of_row(i, j, GFq.random_element())
#
#  logging.debug(f"tmp:\n%s", tmp)
#
##  tmp.add_multiple_of_row(0,n-1,1)
##  tmp.add_multiple_of_row(0,n-2,1)
#  tmp1 = tmp.submatrix(0,0,n-2, m)
#
#  tmp1 = tmp1.rref()
#
#  logging.debug(f"tmp1:\n%s", tmp1)
#
#
#  for i in range(n):
#    for j in range(n):
#      if i != j:
#        tmp.add_multiple_of_row(i, j, GFq.random_element())
#
#
#  tmp2 = tmp.submatrix(0,0,n-2, m)
#
#  tmp2 = tmp2.rref()
#
#  logging.debug(f"tmp2:\n%s", tmp2)


#  tmp = P0prime[0].transpose()
#
#  for i in range(n-2):
#    for j in range(n-1):
#      if i != j:
#        tmp.add_multiple_of_row(i, j, GFq.random_element())
#
##  tmp.add_multiple_of_row(0,n-1,1)
##  tmp.add_multiple_of_row(0,n-2,1)
#  tmp = tmp.submatrix(0,0,n-2, m)
#
#  tmp = tmp.rref()
#
#  logging.debug(f"tmp:\n%s", tmp)
#


#  pivot_row = -1

  for pivot_row in range(m):
    if pivot_row == m-1:
      break 
    if N1[pivot_row, pivot_row] != 1:
#      pivot_row = i
      break

#  print(pivot_row, m-1)

  sol = [0] * (m*m + n*n)

#  if pivot_row > 0:
#    #if partial:
#      for i in range(pivot_row):
#        sol[2*m*n - (m-1) + i] = N1[i, pivot_row]
#
#      sol[2*m*n - (m-1) + pivot_row] = GFq(-1)
#    #else:
#    #  logging.debug(f"no sol")
#    #  return None, None
#  else:
#    for i in range(m-1):
#      sol[2*m*n - (m-1) + i] = N1[i, m-1]
#
#    sol[-1] = GFq(-1)
#
#    pivot_row = m-1

  for i in range(pivot_row):
    sol[2*m*n - (m-1) + i] = N1[i, pivot_row]

  sol[2*m*n - (m-1) + pivot_row] = GFq(-1)

  logging.debug(f"sol:\n%s", sol)

  logging.debug(f"N:\n%s", N)

  for i in range(m):
    sol[2*m*n - m - (m-1) + i] = N[i, pivot_row]

  logging.debug(f"sol:\n%s", sol)


  for c in reversed(range(pivot_row)):
    for r in range(m):
      sol[2*m*n - m - (m-1) + r] -= N[r, c] * sol[2*m*n - (m-1) + c]

  logging.debug(f"sol:\n%s", sol)


  P = -P0prime[1].transpose()

  logging.debug(f"P01nt:\n%s", P)

  for i in range(n):
    sol[m*n + i] = P[i, pivot_row]

  logging.debug(f"sol:\n%s", sol)

  for c in reversed(range(pivot_row)):
    for r in range(n):
      sol[m*n + r] -= P[r, c] * N1[c, pivot_row]

  logging.debug(f"sol:\n%s", sol)


  P = -P0prime[0].transpose()

  logging.debug(f"P00nt:\n%s", P)

  for i in range(n):
    sol[(m-1)*n + i] = P[i, pivot_row]

  logging.debug(f"sol:\n%s", sol)

  for c in reversed(range(pivot_row)):
    for r in range(n):
      sol[(m-1)*n +r] -= P[r, c] * N1[c, pivot_row]

  logging.debug(f"sol:\n%s", sol)

  for b in reversed(range(m-2)):
    for c in reversed(range(m)):
      for r in range(m):
        sol[(m+1)*n + b*m + r] -= N[r, c] * sol[(m+1)*n + b*m + m + c]

  logging.debug(f"sol:\n%s", sol)


  P = -P0prime[0].transpose()

  for b in reversed(range(m-1)):
    for c in reversed(range(m)):
      for r in range(n):
        sol[b*n + r] -=  P[r, c] * sol[2*m*n - (m-1) - (m-1-b)*m + c]

#  logging.debug(f"sol:\n%s", sol)
#
#  sol[-1] = GFq(-1)

  ###################


  logging.debug("sol:\n%s", sol)
  A     = matrix(GFq, m, sol[n*n:])
  B_inv = matrix(GFq, n, sol[:n*n])

  logging.debug(f"A:\n%s", A)
  logging.debug(f"B_inv:\n%s", B_inv)

  return A, B_inv

# def solve_comp(P0prime, Amm):
#   m = P0prime[0].nrows()
#   n = P0prime[0].ncols()
# 
#   GFq = Amm.base_ring()
# 
#   Pj = [None] * 2
# 
#   Pj[0] = matrix(GFq, m, n, [[GFq(1) if i==j else GFq(0) for i in range(n)] for j in range(m)])
#   Pj[1] = matrix(GFq, m, n, [[GFq(1) if i==j else GFq(0) for i in range(n)] for j in range(1,m)] + [[GFq(0)]*n])
# 
#   rsys = []
# 
#   for l in range(n):
#     for j in range(m):
#       tmp = [GFq(0)] * (m^2 + n^2)
# 
#       for i in range(m):
#         tmp[m*m + l*m + i] = - P0prime[0][i,j]
# 
#       for i in range(n):
#         tmp[i*n + j] = Pj[0][l,i]
# 
#       if l == n-1:
#         tmp[m*m + n*n - 1] = Amm * P0prime[0][n-1,j]
# 
#       rsys.append(tmp)
# 
#   for l in range(n):
#     for j in range(m if l < n-1 else m-1):
#       tmp = [GFq(0)] * (m^2 + n^2)
# 
#       for i in range(m):
#         tmp[m*m + l*m + i] = - P0prime[1][i,j]
# 
#       for i in range(n):
#         tmp[i*n + j] = Pj[1][l,i]
# 
#       if l == n-1:
#         tmp[m*m + n*n - 1] = Amm * P0prime[1][n-1,j]
# 
#       rsys.append(tmp)
# 
#   rsys = matrix(GFq, ncols=m^2 + n^2, entries=rsys)
# 
#   rsys_rref = rsys.rref()
# 
#   if not all([rsys_rref[i][i] == 1 for i in range(rsys_rref.nrows())]):
#     #logging.debug("no sol")
#     return None, None
# 
#   sol = rsys_rref.columns()[-1].list()
# 
#   A = matrix(GFq, m, sol[n*n:] + [Amm])
#   B_inv = matrix(GFq, m, sol[:n*n])
# 
#   #logging.debug(f"A:\n%s", A)
#   #logging.debug(f"B_inv:\n%s", B_inv)
# 
#   return A, B_inv


