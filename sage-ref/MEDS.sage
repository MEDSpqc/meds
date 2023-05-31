from Crypto.Hash import SHAKE256
from Crypto.Hash import SHA256

import secrets, time, sys
import logging

import os.path as path

import binascii

import bitstream
from bitstream import BitStream


from SeedTree import SeedTree

import params


class BadSignatureError(Exception):
  pass

class redo(Exception):
  pass
 

class MEDSbase:
  ENDIANESS = 'little'
  
  @staticmethod
  def int_to_bytes(i, length):
      return i.to_bytes(length, MEDSbase.ENDIANESS)
      
  @staticmethod
  def bytes_to_int(bs):
      return int.from_bytes(bs, MEDSbase.ENDIANESS)
  
  @staticmethod
  def matrix_to_bytes(m):
    GFq = m[0,0].base_ring()
    GF_BYTES = ceil(log(GFq.base().order(),2)/8)
    if GFq.is_prime_field():
      return b"".join([MEDSbase.int_to_bytes(int(j), GF_BYTES) for v in m.rows() for j in v])
    else:
      return bytes([j.integer_representation() for v in m.rows() for j in v])
  
  @staticmethod
  def matrix_to_bits(m, bs):
    GFq = m[0,0].base_ring()
    GF_BITS = ceil(log(GFq.base().order(),2))
    if GFq.is_prime_field():
      for v in [j for v in m.rows() for j in v]:
        bs.write(int(v), GF_BITS)
    else:
      raise "Not implemented..."
 
  @staticmethod
  def bytes_to_matrix(GFq, bs, m=0, n=0):
    GF_BYTES = ceil(log(GFq.base().order(),2)/8)
  
    if m == 0:
      m = sqrt(len(bs)/GF_BYTES)
      n = m

    if GFq.is_prime_field():
      data = [MEDSbase.bytes_to_int(bs[i:i + GF_BYTES]) for i in range(0, len(bs), GF_BYTES)]
    else:
      data = [GFq.fetch_int(b) for b in bs]
  
    return matrix(GFq, m, n, data)

  cache = {}

  @staticmethod
  def rnd_GFq(shake, GFq, min=0):
    if GFq in MEDSbase.cache:
      numbytes, mask, q = MEDSbase.cache[GFq]
    else:
      numbytes = ceil(log(GFq.base().order(),2)/8)
      mask = (1 << ceil(log(GFq.base().order(),2))) - 1
      q = GFq.base().order()

      MEDSbase.cache[GFq] = numbytes, mask, q 

    while True:
      val = MEDSbase.bytes_to_int(shake.read(numbytes))

      val = val & mask

      if val >= min and val < q:
        return val

  @staticmethod
  def rnd_GFqs(seed, GFq, length, min=0):
    if type(seed) == SHAKE256.SHAKE256_XOF:
      shake = seed
    else:
      shake = SHAKE256.new()
      shake.update(seed)

    if not isinstance(length, list):
      return [MEDSbase.rnd_GFq(shake, GFq, min) for _ in range(length)]

    data = [MEDSbase.rnd_GFq(shake, GFq, min) for _ in range(sum(length))]

    return tuple(data[sum(length[:i]):sum(length[:i+1])] for i in range(len(length)))

 
  @staticmethod
  def rnd_inv_matrix_pair(seed, GFq, n):
    if type(seed) == SHAKE256.SHAKE256_XOF:
      shake = seed
    else:
      shake = SHAKE256.new()
      shake.update(seed)
  
    while True:
      M = MEDSbase.rnd_matrix(shake, GFq, n, n)
    
      if M.is_invertible():
        return M, M.inverse()

  @staticmethod
  def rnd_inv_matrix(seed, GFq, n):
    if type(seed) == SHAKE256.SHAKE256_XOF:
      shake = seed
    else:
      shake = SHAKE256.new()
      shake.update(seed)
  
    while True:
      M = MEDSbase.rnd_matrix(shake, GFq, n, n)
    
      if M.is_invertible():
        return M
  
    ## L = matrix(GFq, n, n)
  
    ## for i in range(n):
    ##   L[i,i] = 1
    ##   for j in range(i+1, n):
    ##     L[j,i] = MEDSbase.rnd_GFq(shake, GFq)
  
    ## logging.debug(f"L:\n%s", L)
  
    ## U = matrix(GFq, n, n)
  
    ## for i in range(n):
    ##   # get only non-zero values for the diagonal
    ##   U[i, i] = MEDSbase.rnd_GFq(shake, GFq, 1)
  
    ##   for j in range(i+1, n):
    ##     U[i,j] = MEDSbase.rnd_GFq(shake, GFq)
  
    ## logging.debug(f"U:\n%s", U)
  
    ## # There is a small bias - since we do not apply a permutation!
    ## return L*U
  
  @staticmethod
  def rnd_matrix(shake, GFq, m, n=None):
    if not n:
      n = m
  
    return matrix(GFq, m, n, [MEDSbase.rnd_GFq(shake, GFq) for i in range(m*n)])
  
  @staticmethod
  def rnd_sys_matrix(seed, GFq, k, m, n):
    if type(seed) == SHAKE256.SHAKE256_XOF:
      shake = seed
    else:
      shake = SHAKE256.new()
      shake.update(seed)
   
    I_k = matrix.identity(ring=GFq, n=k)
  
    return I_k.augment(MEDSbase.rnd_matrix(shake, GFq, k, m*n-k))


  def __init__(self, params, rng=None):
    self.params = params

    q = params.q
    w = params.w
    t = params.t
    
    self.GFq = GF(q)
    
    self.GF_BYTES = ceil(log(q, 2) / 8)

    self.rng = rng
    
    if self.params.have_seed_tree:
      self.MEDSbase_MAX_PATH_LEN = ((1 << ceil(log(w, 2))) + w * (ceil(log(t, 2)) - ceil(log(w, 2)) - 1))


  def hash(self, m):
    shake = SHAKE256.new()
    shake.update(m)

    return shake.read(self.params.digest_bytes)
  
  def hash_pair(self, value, addr):
    shake = SHAKE256.new()
    shake.update(self.st_salt + value + self.int_to_bytes(addr, 4))
  
    return shake.read(self.params.st_seed_bytes), shake.read(self.params.st_seed_bytes)


  def parse_hash(self, digest):
    t = self.params.t
    s = self.params.s
    w = self.params.w

    logging.debug(f"digest: {[int(v) for v in digest]}")
    logging.debug("digest len: %i", len(digest))
  
    shake = SHAKE256.new()
  
    shake.update(digest)
  
    h = [0] * t

    num = 0

    while num < w:
      pos = MEDSbase.bytes_to_int(shake.read(ceil(log(t,2)/8))) & ((1 << ceil(log(t,2))) - 1)

      if pos >= t:
        continue
  
      if h[pos] > 0:
        continue  
  
      logging.debug(f"pos: {pos}")
  
      val = 0
  
      while val < 1 or val > s-1:
        val = MEDSbase.bytes_to_int(shake.read(ceil(log(s,2)/8))) & ((1 << ceil(log(s,2))) - 1)
  
      h[pos] = val

      logging.debug(f"p: {pos}  v: {val}")

      num += 1

    return h

   
  def seeds_from_path(self, h, path):
    q = self.params.q
    m = self.params.m
    n = self.params.n
    k = self.params.k
    s = self.params.s
    t = self.params.t
    w = self.params.w

    param = self.params

    GFq = self.GFq
    GF_BYTES = ceil(log(q, 2) / 8)

    if self.params.have_seed_tree:
      # prepare tree structure
      seeds = SeedTree(t, param.st_seed_bytes)
  
      for i, v in enumerate(h):
        if v > 0:
          #seeds[i] = (n*n + m*m)*GF_BYTES
          seeds[i] = ceil((n*n + m*m)*log(q, 2)/8)
    
      parsed = []
    
      for i, v in enumerate(seeds.path()):
        parsed.append(path[:v])
    
        path = path[v:]
  
      # apply patch once tree structure is set up
      seeds.patch(parsed, self.hash_pair)
  
    else:
      seeds = []
      for h_i in h:
        if h_i == 0:
          seeds.append(path[:self.params.st_seed_bytes])
          path = path[self.params.st_seed_bytes:]
        else:
          seeds.append(path[:ceil((m*m+n*n)*log(q, 2)/8)])
          path = path[ceil((m*m+n*n)*log(q, 2)/8):]
  
    return seeds


  def pi(self, A, B, G):
    m = self.params.m
    n = self.params.n
    k = self.params.k
  
    GFq = A[0,0].base_ring()
  
    G = [matrix(GFq, m, n, row) for row in G.rows()]
  
    AGB = [A*v*B for v in G]
  
    return matrix(GFq, k, m*n, [AGB[i][j,g] for i in range(k) for j in range(m) for g in range(n)])

  def randombytes(self, nbytes):
    if self.rng:
      return self.rng.randombytes(nbytes)
    else:
      return secrets.token_bytes(nbytes)

  def XOF(self, seed, length):
    shake = SHAKE256.new()
    shake.update(seed)

    if not isinstance(length, list):
      return shake.read(length)

    data = shake.read(sum(length))

    return tuple(data[sum(length[:i]):sum(length[:i+1])] for i in range(len(length)))

  def solve(self, G0, Tj, Amm):
    q = self.params.q
    m = self.params.m
    n = self.params.n
    k = self.params.k
    s = self.params.s
    t = self.params.t
    w = self.params.w

    GFq = self.GFq

    G0prime = Tj * G0
    
    logging.debug(f"G0prime:\n%s", G0prime)

    P0prime = [matrix(GFq, m, n, row) for row in G0prime.rows()]

    N = -P0prime[0].transpose()

    logging.debug(f"N:\n%s", N)


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

    logging.debug(f"M:\n%s", M)

    _, _, M = M.LU()

    for i in range(n):
      M.set_row_to_multiple_of_row(i, i, 1/M[i,i])
    logging.debug(f"M:\n%s", M) 

    M[-1, -1] = GFq(0)
    logging.debug(f"M:\n%s", M) 

    M = M.rref()
    logging.debug(f"M:\n%s", M) 

    sol = [0] * (n*n + m*m)

    sol[-1] = Amm

    for i in range(n-1):
      sol[n*n + m*m - n + i] = M[i, -1]

    for i in range(n):
      sol[n*n + m*m - 2*n + i] = M[i, -2]

    for i in range(n):
      sol[n*n - n + i] = P0prime[0][m-1,i]*Amm

    logging.debug(f"initial sol:\n%s", sol)

    # incomplete blocks:

    for i in range(n):
      for j in range(n-1):
        sol[n*n + m*m - 2*n + i] += - M[i, n+n-2-j] * sol[n*n + m*m - 2 - j]
        sol[n*n - n + i] += - N[i, m-2-j] * sol[n*n + m*m - 2 - j]

    logging.debug(f"incomplete blocks:\n%s", sol)


    # complete blocks:

    for block in range(3,n+1):
      for i in range(n):
        for j in range(n):
          sol[n*n + m*m - block*n + i] += - M[i, n+n-1-j] * sol[n*n + m*m - 1 - (block-2)*n - j]

    for block in range(2,n+1):
      for i in range(n):
        for j in range(n):
          sol[n*n - block*n + i] += - N[i, m-1-j] * sol[n*n + m*m - 1 - (block-1)*n - j]

    logging.debug(f"complete blocks:\n%s", sol)


    A = matrix(GFq, m,m, sol[n*n:])
    B_inv = matrix(GFq, n,n, sol[:n*n])

    logging.debug(f"A:\n%s", A)
    logging.debug(f"B_inv:\n%s", B_inv)

    return A, B_inv


  def solve_med(self, G0, Tj, Amm):
    q = self.params.q
    m = self.params.m
    n = self.params.n
    k = self.params.k
    s = self.params.s
    t = self.params.t
    w = self.params.w

    GFq = self.GFq

    G0prime = Tj * G0
    
    logging.debug(f"G0prime:\n%s", G0prime)

    P0prime = [matrix(GFq, m, n, row) for row in G0prime.rows()]

    AP0 = -P0prime[0].transpose()
    AP1 = -P0prime[1].transpose()
    
    logging.debug(f"AP0:\n%s", AP0)
    logging.debug(f"AP1:\n%s", AP1)

    rsys_top = matrix(GFq, n*m, m*n)
    rsys_bot = matrix(GFq, n*m - 1, m*n)

 
    for block in range(m):
      for row in range(m):
        for col in range(n if block < m-1 else n-1):
          rsys_top[block*n+row, block*m+col] = AP0[row, col]
    
    for row in range(m):
      rsys_top[(m-1)*n+row, (m-1)*m+n-1] = -Amm * AP0[row, n-1]
    
    
    for block in range(m):
      for row in range(m if block < m-1 else m-1):
        for col in range(n if block < m-1 else n-1):
          rsys_bot[block*n+row, block*m+col] = AP1[row, col]
    
    for row in range(m-1):
      rsys_bot[(m-1)*n+row, (m-1)*m+n-1] = -Amm * AP1[row, n-1]
    
    
    for block in range(m-1):
      for row in range(m):
        for col in range(n if block < m-2 else n-1):
          rsys_bot[block*n+row, (block+1)*m+col] = -AP0[row, col]
    
    for row in range(m):
          rsys_bot[(m-2)*n+row, (m-1)*m+n-1] = Amm * AP0[row, n-1]
    
    logging.debug(f"sys top:\n%s", rsys_top)
    logging.debug(f"sys bottom:\n%s", rsys_bot)
    
    rsys_bot = rsys_bot.rref()
    
    logging.debug(f"sys bottom rref:\n%s", rsys_bot.rref())

    for col in reversed(range(m*n-1)):
      for row in range(n*m):
        rsys_top[row, m*n-1] -= rsys_top[row, col] * rsys_bot[col, m*n - 1]
    
    logging.debug(f"sys top elim:\n%s", rsys_top)

    sol = list(rsys_top.column(m*n-1)) + list(rsys_bot.column(m*n-1))
    
    logging.debug(f"sol:\n%s", sol)

    A = matrix(GFq, n, sol[m*m:] + [Amm])
    B_inv = matrix(GFq, m, sol[:m*m])

    return A, B_inv

  def solve_slow(self, G0, Tj, Amm):
    q = self.params.q
    m = self.params.m
    n = self.params.n
    k = self.params.k
    s = self.params.s
    t = self.params.t
    w = self.params.w

    GFq = self.GFq

    G0prime = Tj * G0
    
    logging.debug(f"G0prime:\n%s", G0prime)

    P0prime = [matrix(GFq, m, n, row) for row in G0prime.rows()]

    Pj = [None] * 2

    Pj[0] = matrix(GFq, m, n, [[GFq(1) if i==j else GFq(0) for i in range(n)] for j in range(m)])
    logging.debug(f"Pj[0]:\n%s", Pj[0].str(rep_mapping=lambda x : f'{int(x):4}'))

    Pj[1] = matrix(GFq, m, n, [[GFq(1) if i==j else GFq(0) for i in range(n)] for j in range(1,m)] + [[GFq(0)]*n]) #[sec_GFqs[0]])
    logging.debug(f"Pj[1]:\n%s", Pj[1].str(rep_mapping=lambda x : f'{int(x):4}'))

    rsys = []

    for l in range(n):
      for j in range(m):
        tmp = [GFq(0)] * (m^2 + n^2)
        
        for ii in range(m):
          tmp[l*m + ii] = P0prime[0][ii,j]
        
        for ii in range(n if j < m-1 else n-1):
          tmp[m*m + ii*n + j] = - Pj[0][l,ii]
        
        if j == m-1:
          tmp[m*m + n*n - 1] = Amm * Pj[0][l,n-1]
      
        rsys.append(tmp)
    
    for l in range(n):
      for j in range(m if l < n-1 else m-1):
        tmp = [GFq(0)] * (m^2 + n^2)
        
        for ii in range(m):
          tmp[l*m + ii] = P0prime[1][ii,j]
        
        for ii in range(n if j < m-1 else n-1):
          tmp[m*m + ii*n + j] = - Pj[1][l,ii]
        
        if j == m-1:
          tmp[m*m + n*n - 1] = Amm * Pj[1][l,n-1]
      
        rsys.append(tmp)
    
    rsys = matrix(GFq, ncols=m^2 + n^2, entries=rsys)

    logging.debug(f"rsys:\n%s", rsys)

    rsys_rref = rsys.rref()

    if not all([rsys_rref[i][i] == 1 for i in range(rsys_rref.nrows())]):
      logging.debug("no sol")
      return None, None

    sol = rsys_rref.columns()[-1].list()
    
    logging.debug(f"sol:\n%s", sol)
    
    A = matrix(GFq, m, sol[:m*m])
    B_inv = matrix(GFq, n, sol[m*m:] + [Amm])

    return A, B_inv
 
  def crypto_sign_keypair(self):
    global write_rest

    q = self.params.q
    m = self.params.m
    n = self.params.n
    k = self.params.k
    s = self.params.s
    t = self.params.t
    w = self.params.w

    param = self.params

    GFq = self.GFq

    GF_BYTES = ceil(log(q, 2) / 8)
    GF_BITS = ceil(log(q, 2))

    root_seed = self.randombytes(param.sec_seed_bytes)

    A_inv = [None] * s
    B_inv = [None] * s
    G = [None] * s
  

    pub_seed0, sec_seed, root_seed = self.XOF(root_seed, [param.pub_seed_bytes, param.sec_seed_bytes, param.sec_seed_bytes])

    logging.debug(f"sec_seed:\n{[int(i) for i in sec_seed]}")
    logging.debug(f"pub_seed:\n{[int(i) for i in pub_seed0]}")

    if self.params.have_pk_opt:
      Gi_seed, pub_seed = self.XOF(pub_seed0, [param.pub_seed_bytes, param.pub_seed_bytes])
      G[0] = MEDSbase.rnd_sys_matrix(Gi_seed, GFq, k, m, n)
    else:
      G[0] = MEDSbase.rnd_sys_matrix(pub_seed0, GFq, k, m, n)
    
    logging.debug(f"G[0]:\n%s", G[0])
    
    for i in range(1, s):
      while True: # repeat until A[i] and B[i] are invertible and G[i] has systematic form
        if self.params.have_pk_opt:
          sec_seed_GFqs, sec_seed_T, sec_seed = self.XOF(sec_seed, [param.sec_seed_bytes, param.sec_seed_bytes, param.sec_seed_bytes])

          Tj = MEDSbase.rnd_inv_matrix(sec_seed_T, GFq, k)
          logging.debug(f"Tj:\n%s", Tj)
    
          sec_GFqs = MEDSbase.rnd_GFqs(sec_seed_GFqs, GFq, [1])

          Amm = sec_GFqs[0][0]
          logging.debug(f"Amm:\n{Amm}")

          A, B_inv[i] = self.solve(G[0], Tj, Amm)

          if not B_inv[i].is_invertible():
            logging.debug("no B")
            continue  # try agian for this index
    
          if not A.is_invertible():
            logging.debug("no A_inv")
            continue  # try agian for this index
          
          B = B_inv[i].inverse()
          A_inv[i] = A.inverse()
        else:
          logging.debug(f"sec_seed[s]:\n{[int(i) for i in sec_seed[i]]}")

          sec_seed_A, sec_seed_B, sec_seed = self.XOF(sec_seed, \
                  [param.sec_seed_bytes, param.sec_seed_bytes, param.sec_seed_bytes])

          logging.debug(f"sec_seed_A:\n{[int(i) for i in sec_seed_A]}")
          logging.debug(f"sec_seed_B:\n{[int(i) for i in sec_seed_B]}")

          A, A_inv[i] = MEDSbase.rnd_inv_matrix_pair(sec_seed_A, GFq, m)
          B, B_inv[i] = MEDSbase.rnd_inv_matrix_pair(sec_seed_B, GFq, n)
    
        logging.debug(f"A[{i}]:\n%s", A)
        logging.debug(f"A_inv[{i}]:\n%s", A_inv[i])
        logging.debug(f"B[{i}]:\n%s", B)
        logging.debug(f"B_inv[{i}]:\n%s", B_inv[i])
    
        G[i] = self.pi(A, B, G[0])
    
        G[i] = G[i].echelon_form()
    
        # check if we got systematic form
        if sum([G[i][j,j] for j in range(k)]) == k:
          logging.debug(f"G[{i}]:\n%s", G[i])
          break

        # if no systematic form loop to try again for this index
        logging.debug(f"redo G[{i}]")

    bs = bitstream.BitStream()

    for v in [j for A_inv_i in A_inv[1:] for row in A_inv_i for j in row]:
      bs.write(int(v), GF_BITS)

    bs.finalize()

    for v in [j for B_inv_i in B_inv[1:] for row in B_inv_i for j in row]:
      bs.write(int(v), GF_BITS)

    bs.finalize()

    for v in [j for row in G[0][:,k:] for j in row]:
      bs.write(int(v), GF_BITS)

    bs.finalize()

    comp = bs.data

    sk = root_seed + comp

    logging.debug(f"sk:\n0x{binascii.hexlify(sk).decode()}")


    G = [[j for v in Gi[:,k:].rows() for j in v] for Gi in G]

    if self.params.have_pk_opt:
      G = [Gi[(m*n - k + (m-1)*n-k):] for Gi in G]

    bs = bitstream.BitStream()

    for v in [j for Gi in G[1:] for j in Gi]:
      bs.write(int(v), GF_BITS)

    comp = bs.data

    pk = pub_seed0 + comp

    logging.debug(f"pub_seed (pk):\n{[int(i) for i  in pub_seed0]}")
    logging.debug(f"G[1:] (pk):\n{[int(j) for j in comp]}")

    logging.debug(f"pk:\n0x{binascii.hexlify(pk).decode()}")
  
    return sk, pk
  
  
  def crypto_sign(self, sk, msg):
    q = self.params.q
    m = self.params.m
    n = self.params.n
    k = self.params.k
    s = self.params.s
    t = self.params.t
    w = self.params.w

    param = self.params

    GFq = self.GFq
    GF_BYTES = ceil(log(q, 2) / 8)
    GF_BITS = ceil(log(q, 2))

    initial_seed = self.randombytes(param.sec_seed_bytes)

    # skip secret key seed
    sk = sk[param.sec_seed_bytes:]
  
    A_inv = [None]*s
    B_inv = [None]*s

    bs = bitstream.BitStream(sk)

    sk = []
    for i in range(floor(len(bs.data)*8/GF_BITS)):
      sk.append(GFq(bs.read(GF_BITS)))

    for i in range(1, s):
      A_inv[i] = matrix(GFq, m, m, sk[:m*m]); sk = sk[m*m:]
      logging.debug(f"A_inv[i]:\n%s", A_inv[i])
  
    for i in range(1, s):
      B_inv[i] = matrix(GFq, n, n, sk[:n*n]); sk = sk[n*n:]
      logging.debug(f"B_inv[i]:\n%s", B_inv[i])
  
    G_0 = matrix.identity(ring=GFq, n=k).augment(matrix(GFq, k, m*n-k, sk))
  
    logging.debug(f"G_0:\n%s", G_0)
    logging.debug(f"sig_seed:\n{[int(i) for i in initial_seed]}")
  
    if self.params.have_seed_tree:
      st_root, self.st_salt = self.XOF(initial_seed, [param.st_seed_bytes, param.st_salt_bytes])
  
      seeds = SeedTree(t, st_root, self.hash_pair)
    else:
      seeds = list(self.XOF(initial_seed, [param.st_seed_bytes] * t))
  
    A_tilde = [None]*t
    B_tilde = [None]*t
    G_tilde = [None]*t
  
    for i in range(t):
      seed = seeds[i]
  
      while True:
        pub_seed_A_tilde, pub_seed_B_tilde, seed = self.XOF(seed, [param.pub_seed_bytes, param.pub_seed_bytes, param.st_seed_bytes])

        A_tilde[i] = MEDSbase.rnd_inv_matrix(pub_seed_A_tilde, self.GFq, m)
        B_tilde[i] = MEDSbase.rnd_inv_matrix(pub_seed_B_tilde, self.GFq, n)
  
        logging.debug(f"A_tilde[{i}]:\n%s", A_tilde[i])
        logging.debug(f"B_tilde[{i}]:\n%s", B_tilde[i])
  
        G_tilde[i] = self.pi(A_tilde[i], B_tilde[i], G_0)
  
        logging.debug(f"G_tilde[{i}]:\n%s", G_tilde[i])
  
        G_tilde[i] = G_tilde[i].echelon_form()  # 'public' computation - no need for constant time
  
        # check if we got systematic form
        if sum([G_tilde[i][j,j] for j in range(k)]) == k:
          logging.debug(f"G_tilde[{i}]:\n%s", G_tilde[i])
  
          G_tilde[i] = G_tilde[i][:,k:]
          break

        if not self.params.have_seed_tree:
          seeds[i] = seed

    if type(msg) == str:
      msg = msg.encode('utf8')
  
    digest = self.hash(b"".join([MEDSbase.matrix_to_bytes(M) for M in G_tilde]) + msg)
  
    logging.debug(f"digest:\n{[int(i) for i in digest]}")

    h = self.parse_hash(digest) #, t, s, w)
  
    logging.debug(f"h:\n{h}")
  
    for i, h_i in enumerate(h):
      if h_i > 0:
        mu = A_tilde[i] * A_inv[h_i]
        nu = B_inv[h_i] * B_tilde[i]
  
        logging.debug(f"mu:\n%s", mu)
        logging.debug(f"nu:\n%s", nu)
  
        bs = BitStream()

        MEDSbase.matrix_to_bits(mu, bs)
        MEDSbase.matrix_to_bits(nu, bs)

        seeds[i] = bs.data

    if self.params.have_seed_tree:
      ret = b"".join(seeds.path())
      ret += bytes([0]) * (param.sig_size - param.digest_bytes - param.st_salt_bytes - len(ret))
      ret += digest
      ret += self.st_salt
    else:
      ret = b"".join(seeds) + digest

    ret += msg
  
    logging.debug(f"sm:\n0x{binascii.hexlify(ret).decode()}")
  
    return ret
  
  def crypto_sign_open(self, pk, sig):
    q = self.params.q
    m = self.params.m
    n = self.params.n
    k = self.params.k
    s = self.params.s
    t = self.params.t
    w = self.params.w

    param = self.params

    GFq = self.GFq
    GF_BYTES = ceil(log(q, 2) / 8)
    GF_BITS = ceil(log(q, 2))

    I_k = matrix.identity(ring=GFq, n=k)
  
    G = []
  
    pub_seed = pk[:param.pub_seed_bytes]; pk = pk[param.pub_seed_bytes:]

#    if self.params.have_pk_opt:
#      pub_G_seed = self.XOF(pub_seed, [param.pub_seed_bytes] * s)
#    else:
#      pub_G_seed = [pub_seed]

    bs = bitstream.BitStream(pk)

    pk = []
    for i in range(floor(len(bs.data)*8/GF_BITS)):
      pk.append(GFq(bs.read(GF_BITS)))

    if self.params.have_pk_opt:
      Gi_seed, pub_seed = self.XOF(pub_seed, [param.pub_seed_bytes, param.pub_seed_bytes])
      G.append(MEDSbase.rnd_sys_matrix(Gi_seed, GFq, k, m, n))
    else:
      G.append(MEDSbase.rnd_sys_matrix(pub_seed, GFq, k, m, n))
  
    #while len(pk) > 0:
    for Gi in range(1, s):
      if self.params.have_pk_opt:
        Graw = [GFq(0)]*(m*n - k + (m-1)*n-k) \
             + pk[:k*(m*n-k) - (m*n - k + (m-1)*n-k)]
        Gright = matrix(GFq, k, m*n-k, Graw)
        pk = pk[k*(m*n-k) - (m*n - k + (m-1)*n-k):]

        Gi_seed, pub_seed = self.XOF(pub_seed, [param.pub_seed_bytes, param.pub_seed_bytes])
        data = MEDSbase.rnd_GFqs(Gi_seed, GFq, [m*n-k, (m-1)*n-k])

        for i in range(1,m):
          for j in range(n):
            Gright[0,(i-1)*n+j] = GFq(1) if i==j else GFq(0)

#        for i in range(m*n - k):
#          Gright[0,i] = data[0][i]
            
        data = [0, 1]
        data[1] = [GFq(1) if i==j else GFq(0) for j in range(2,m) for i in range(n)] #+ data[1][-n:]

        for i in range((m-1)*n-n):
          #Gright[1,i] = GFq(1) if data[1][i] > (q>>1) else GFq(0)
          Gright[1,i] = data[1][i]

      else:
         Gright = matrix(GFq, k, m*n-k, pk[:k*(m*n-k)])
         pk = pk[k*(m*n-k):]

      G.append(I_k.augment(Gright))

 
    for i,g in enumerate(G):
      logging.debug(f"G[{i}]:\n%s", g)
  
    if self.params.have_seed_tree:
      digest = sig[param.sig_size-param.digest_bytes-param.st_salt_bytes:param.sig_size-param.st_salt_bytes]
      self.st_salt = sig[param.sig_size-param.st_salt_bytes:param.sig_size]
      path = sig[:self.params.sig_size-param.digest_bytes-param.st_salt_bytes]
    else:
      digest = sig[param.sig_size-param.digest_bytes:param.sig_size]
      path = sig[:self.params.sig_size-param.digest_bytes]
  
    logging.debug(f"digest:\n{[int(i) for i in digest]}")

    h = self.parse_hash(digest) #, t, s, w)
  
    seeds = self.seeds_from_path(h, path)
  
    G_hat = [None] * t
  
    for i, h_i in enumerate(h):
      logging.debug(f"seeds[{i}]:\n{[int(v) for v in seeds[i]]}")

      if h_i == 0:
  
        seed = seeds[i]
    
        while True:
          pub_seed_A_tilde, pub_seed_B_tilde, seed = self.XOF(seed, [param.pub_seed_bytes, param.pub_seed_bytes, param.st_seed_bytes])
  
          logging.debug(f"pub_seed_A_tilde:\n{[int(i) for i in pub_seed_A_tilde]}")
          A_tilde = MEDSbase.rnd_inv_matrix(pub_seed_A_tilde, self.GFq, m)

          logging.debug(f"pub_seed_B_tilde:\n{[int(i) for i in pub_seed_B_tilde]}")
          B_tilde = MEDSbase.rnd_inv_matrix(pub_seed_B_tilde, self.GFq, n)

          logging.debug(f"A_tilde[{i}]:\n%s", A_tilde)
          logging.debug(f"B_tilde[{i}]:\n%s", B_tilde)
    
          G_hat[i] = self.pi(A_tilde, B_tilde, G[0])
  
          logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])
  
          G_hat[i] = G_hat[i].echelon_form()
  
          logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])
  
          if not self.params.have_seed_tree:
            break
    
          # check if we got systematic form
          if sum([G_hat[i][j,j] for j in range(k)]) == k:
            break
    
      else:
        bs = BitStream(seeds[i])

        seeds[i] = []
        for _ in range(floor(len(bs.data)*8/GF_BITS)):
          seeds[i].append(GFq(bs.read(GF_BITS)))


        mu = matrix(GFq, m, m, seeds[i][:m*m])
        nu = matrix(GFq, n, n, seeds[i][m*m:])
  
        logging.debug(f"mu[{i}]:\n%s", mu)
        logging.debug(f"nu[{i}]:\n%s", nu)
   
        G_hat[i] = self.pi(mu, nu, G[h_i])
  
        logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])
  
        G_hat[i] = G_hat[i].echelon_form()
  
        logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])
  
  
    check = self.hash(b"".join([MEDSbase.matrix_to_bytes(M[:,k:]) for M in G_hat]) + \
                      sig[self.params.sig_size:])
  
    if not digest == check:
      raise BadSignatureError("Signature verification failed!")
  
    return sig[self.params.sig_size:]

class MEDS(MEDSbase):
  def __init__(self, param, rng = None):
    if type(param) == str:
      exec(f"MEDSbase.__init__(self, params.{param}, rng)")
    else:
      MEDSbase.__init__(self, param, rng)

  @property
  def sk(self):
    return self._sk

  @property
  def pk(self):
    return self._pk

  def crypto_sign_keypair(self):
    self._sk, self._pk = MEDSbase.crypto_sign_keypair(self)

  def crypto_sign(self, msg):
    return MEDSbase.crypto_sign(self, self._sk, msg)

  def crypto_sign_open(self, sig):
    return MEDSbase.crypto_sign_open(self, self._pk, sig)


for p in params.params:
  exec(f"""
class {p.name}(MEDS):
  def __init__(self):
    MEDS.__init__(self, params.{p.name})""")


if __name__ == "__main__":

  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument("-v", "--verbous", action='store_true', help = "verbous debugging output")
  parser.add_argument("-l", "--list", action='store_true', help = "list parameter sets")
  parser.add_argument('parset', nargs='?', default='toy', help = "selected parameter set")

  args = parser.parse_args()

  if args.list:
    print("Available parameter sets:\n")
    params.print_list()
    sys.exit()

  formatter = logging.Formatter("(%(funcName)s) %(message)s\n")

  # Log to pipe 3 if it exists, e.g., "KAT_gen.sage 3>&1".
  if path.exists("/proc/self/fd/3"):
    handler = logging.FileHandler(filename = "/proc/self/fd/3", mode='w')
    handler.setFormatter(formatter)
  
    logger = logging.getLogger()
    logger.setLevel("DEBUG")
    logger.addHandler(handler)

  if args.verbous:
    handler = logging.FileHandler(filename = "KAT/" + args.parset + ".log", mode='w')
    handler.setFormatter(formatter)
  
    logger = logging.getLogger()
    logger.setLevel("DEBUG")
    logger.addHandler(handler)

  import randombytes

  # Test and benchmark:
  try:
    randombytes.randombytes_init(b"deadbeeffeedbeeffeeddeadbeefbeefdeadbeeffeedbeef", None, 256)

    meds = MEDS(args.parset, randombytes)

    print(f"parameter set: {meds.params.name}\n")

    start_time = time.time()
    meds.crypto_sign_keypair()
    pk_time = time.time() - start_time
    
    print("pk:  {0} bytes".format(len(meds.pk)))

#assert(len(meds.pk) == meds.params.pk_size)
    
    msg = b"TestTestTestTest\0"
    
    start_time = time.time()
    sm = meds.crypto_sign(msg)
    sig_time = time.time() - start_time

    print("sig: {} bytes".format(len(sm)))

#assert(len(sm)-len(msg) == meds.params.sig_size)

    ## break signature for testing:
    # sig = bytearray(sig); sig[54] = 8
  
    start_time = time.time()
    tmp = meds.crypto_sign_open(sm)
    verify_time = time.time() - start_time

    print(tmp == msg)

    print()
    print("success")
    print()
 
    print("Time:")
    print("keygen: {0:1.4f}".format(pk_time))
    print("sign:   {0:1.4f}".format(sig_time))
    print("verify: {0:1.4f}".format(verify_time))

    print()

    print(f"FS security:   {meds.params.fiat_shamir:12.5f}")

  except BadSignatureError as e:
    print()
    print("!!! FAILURE !!!")
    print()
    print("Exception:", e)
    print()
 
