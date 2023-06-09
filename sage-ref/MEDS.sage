#!/usr/bin/env sage

from Crypto.Hash import SHAKE256
from Crypto.Hash import SHA256

import secrets, time, sys
import logging

import os.path as path

import binascii

import bitstream
from bitstream import BitStream


from util import *

from SeedTree import SeedTree

import params


class BadSignatureError(Exception):
  pass


class MEDSbase:
  def __init__(self, params, rng=secrets.token_bytes):
    self.params = params

    self.randombytes = rng

  def crypto_sign_keypair(self):
    q = self.params.q
    m = self.params.m
    n = self.params.n
    k = self.params.k
    s = self.params.s
    t = self.params.t
    w = self.params.w

    param = self.params

    GFq = GF(q)

    GF_BYTES = ceil(log(q, 2) / 8)
    GF_BITS = ceil(log(q, 2))

    root_seed = self.randombytes(param.sec_seed_bytes)

    A_inv = [None] * s
    B_inv = [None] * s
    G = [None] * s


    pub_seed, sec_seed, root_seed = XOF(root_seed, [param.pub_seed_bytes, param.sec_seed_bytes, param.sec_seed_bytes])

    logging.debug(f"sigma:\n%s", [int(i) for i in sec_seed])
    logging.debug(f"sigma_G0:\n%s", [int(i) for i in pub_seed])

    G[0] = rnd_sys_matrix(pub_seed, GFq, k, m, n)

    logging.debug(f"G[0]:\n%s", G[0])

    for i in range(1, s):
      while True: # repeat until A[i] and B[i] are invertible and G[i] has systematic form
        sec_seed_GFqs, sec_seed_T, sec_seed = XOF(sec_seed, [param.sec_seed_bytes, param.sec_seed_bytes, param.sec_seed_bytes])

        Ti = rnd_inv_matrix(sec_seed_T, GFq, k)
        logging.debug(f"Ti:\n%s", Ti)

        Amm = rnd_GFq(sec_seed_GFqs, GFq)
        logging.debug(f"Amm:\n{Amm}")

        G0prime = Ti * G[0]

        logging.debug(f"G0prime:\n%s", G0prime)


        A, B_inv[i] = solve_symb([matrix(GFq, m, n, G0prime.rows()[i]) for i in range(2)], Amm)

        if A == None and B_inv[i] == None:
          logging.debug("no sol")
          continue  # try agian for this index

        if not B_inv[i].is_invertible():
          logging.debug("no B")
          continue  # try agian for this index

        if not A.is_invertible():
          logging.debug("no A_inv")
          continue  # try agian for this index

        B = B_inv[i].inverse()
        A_inv[i] = A.inverse()

        logging.debug(f"A[{i}]:\n%s", A)
        logging.debug(f"A_inv[{i}]:\n%s", A_inv[i])
        logging.debug(f"B[{i}]:\n%s", B)
        logging.debug(f"B_inv[{i}]:\n%s", B_inv[i])

        G[i] = pi(A, B, G[0])

        G[i] = G[i].echelon_form()

        # check if we got systematic form
        if sum([G[i][j,j] for j in range(k)]) == k:
          logging.debug(f"G[{i}]:\n%s", G[i])
          break

        # if no systematic form loop to try again for this index
        logging.debug(f"redo G[{i}]")


    G = [[j for v in Gi[:,k:].rows() for j in v] for Gi in G]

    G = [Gi[(m*n - k + (m-1)*n-k):] for Gi in G]

    bs = bitstream.BitStream()

    for Gi in G[1:]:
      for v in [j  for j in Gi]:
        bs.write(int(v), GF_BITS)

      bs.finalize()

    comp = bs.data

    pk = pub_seed + comp

    logging.debug(f"sigma_G0 (pk):\n%s", [int(i) for i in pub_seed])
    logging.debug(f"G[1:] (pk):\n%s", [int(j) for j in comp])

    logging.debug(f"pk:\n0x%s", binascii.hexlify(pk).decode())


    bs = bitstream.BitStream()

    for A_inv_i in A_inv[1:]:
      for v in [j for row in A_inv_i for j in row]:
        bs.write(int(v), GF_BITS)

      bs.finalize()

    for B_inv_i in B_inv[1:]:
      for v in [j for row in B_inv_i for j in row]:
        bs.write(int(v), GF_BITS)

      bs.finalize()

    sk = root_seed + pub_seed + bs.data

    logging.debug(f"sk:\n0x%s", binascii.hexlify(sk).decode())

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

    GFq = GF(q)
    GF_BYTES = ceil(log(q, 2) / 8)
    GF_BITS = ceil(log(q, 2))

    initial_seed = self.randombytes(param.sec_seed_bytes)

    # skip secret key seed
    sk = sk[param.sec_seed_bytes:]

    pub_seed = sk[:param.pub_seed_bytes]
    sk = sk[param.pub_seed_bytes:]

    A_inv = [None]*s
    B_inv = [None]*s

    bs = bitstream.BitStream(sk)

    for i in range(1, s):
      data = [GFq(bs.read(GF_BITS)) for _ in range(m*m)]
      bs.finalize()

      A_inv[i] = matrix(GFq, m, m, data)
      logging.debug(f"A_inv[i]:\n%s", A_inv[i])

    for i in range(1, s):
      data = [GFq(bs.read(GF_BITS)) for _ in range(n*n)]
      bs.finalize()

      B_inv[i] = matrix(GFq, n, n, data)
      logging.debug(f"B_inv[i]:\n%s", B_inv[i])

    G_0 = rnd_sys_matrix(pub_seed, GFq, k, m, n)

    logging.debug(f"G_0:\n%s", G_0)
    logging.debug(f"delta:\n%s", [int(i) for i in initial_seed])

    st_root, self.st_salt = XOF(initial_seed, [param.st_seed_bytes, param.st_salt_bytes])
    seeds = SeedTree(t, st_root, G(self.params), self.st_salt)

    for i in range(t):
      logging.debug(f"sigma[{i}]:\n0x%s", binascii.hexlify(seeds[i]).decode())

    A_tilde = [None]*t
    B_tilde = [None]*t
    G_tilde = [None]*t

    for i in range(t):
      seed = seeds[i]

      while True:
        pub_seed_A_tilde, pub_seed_B_tilde, seed = XOF(seed, [param.st_seed_bytes, param.st_seed_bytes, param.st_seed_bytes])

        logging.debug(f"sigma_A_tilde[{i}]:\n0x%s", binascii.hexlify(pub_seed_A_tilde).decode())

        A_tilde[i] = rnd_inv_matrix(pub_seed_A_tilde, GFq, m)
        B_tilde[i] = rnd_inv_matrix(pub_seed_B_tilde, GFq, n)

        logging.debug(f"A_tilde[{i}]:\n%s", A_tilde[i])
        logging.debug(f"B_tilde[{i}]:\n%s", B_tilde[i])

        G_tilde[i] = pi(A_tilde[i], B_tilde[i], G_0)

        logging.debug(f"G_tilde[{i}]:\n%s", G_tilde[i])

        G_tilde[i] = G_tilde[i].echelon_form()  # 'public' computation - no need for constant time

        # check if we got systematic form
        if sum([G_tilde[i][j,j] for j in range(k)]) == k:
          logging.debug(f"G_tilde[{i}]:\n%s", G_tilde[i])

          break

    if type(msg) == str:
      msg = msg.encode('utf8')

    comp = bytes()

    for G_tilde_i in G_tilde:
      bs = BitStream()
      compress(G_tilde_i[:,k:], bs)

      comp += bs.data

    digest = H(self.params)(comp + msg)

    logging.debug(f"digest:\n%s", [int(i) for i in digest])

    h = parse_hash(digest, self.params)

    logging.debug(f"h:\n{h}")

    bs = BitStream()

    for i, h_i in enumerate(h):
      if h_i > 0:
        mu = A_tilde[i] * A_inv[h_i]
        nu = B_inv[h_i] * B_tilde[i]

        logging.debug(f"mu:\n%s", mu)
        logging.debug(f"nu:\n%s", nu)

        compress(mu, bs)
        bs.finalize()

        compress(nu, bs)
        bs.finalize()

        seeds[i] = None

    ret = bs.data

    ret += b"".join(seeds.path())
    ret += bytes([0]) * (param.sig_size - param.digest_bytes - param.st_salt_bytes - len(ret))
    ret += digest
    ret += self.st_salt

    ret += msg

    logging.debug(f"sm:\n0x%s", binascii.hexlify(ret).decode())

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

    GFq = GF(q)
    GF_BYTES = ceil(log(q, 2) / 8)
    GF_BITS = ceil(log(q, 2))

    logging.debug(f"pk:\n0x%s", binascii.hexlify(pk).decode())
    logging.debug(f"sm:\n0x%s", binascii.hexlify(sig).decode())

    I_k = matrix.identity(ring=GFq, n=k)

    G = []

    pub_seed = pk[:param.pub_seed_bytes]

    bs = bitstream.BitStream(pk[param.pub_seed_bytes:])

    G.append(rnd_sys_matrix(pub_seed, GFq, k, m, n))

    for Gi in range(1, s):
      data = [GFq(bs.read(GF_BITS)) for _ in range(n + (k-2)*(m*n-k))]
      bs.finalize()

      Graw = [GFq(0)]*(m*n - k + (m-1)*n-k) + data
      Gright = matrix(GFq, k, m*n-k, Graw)

      Gi_seed, pub_seed = XOF(pub_seed, [param.pub_seed_bytes, param.pub_seed_bytes])

      for i in range(1,m):
        for j in range(n):
          Gright[0,(i-1)*n+j] = GFq(1) if i==j else GFq(0)

      data = [GFq(1) if i==j else GFq(0) for j in range(2,m) for i in range(n)] #+ data[1][-n:]

      for i in range((m-1)*n-n):
        Gright[1,i] = data[i]

      G.append(I_k.augment(Gright))


    for i,g in enumerate(G):
      logging.debug(f"G[{i}]:\n%s", g)

    munu = sig[:(ceil(m*m * GF_BITS / 8) + ceil(n*n * GF_BITS / 8))*w]
    sig = sig[(ceil(m*m * GF_BITS / 8) + ceil(n*n * GF_BITS / 8))*w:]
    path = sig[:param.seed_tree_cost]; sig = sig[param.seed_tree_cost:]
    digest = sig[:param.digest_bytes]; sig = sig[param.digest_bytes:]
    self.st_salt = sig[:param.st_salt_bytes]; msg = sig[param.st_salt_bytes:]

    logging.debug(f"munu:\n0x%s", binascii.hexlify(munu).decode())
    logging.debug(f"path:\n0x%s", binascii.hexlify(path).decode())
    logging.debug(f"digest:\n0x%s", binascii.hexlify(digest).decode())
    logging.debug(f"alpha:\n0x%s", binascii.hexlify(self.st_salt).decode())

    h = parse_hash(digest, self.params)

    seeds = seeds_from_path(h, path, self.st_salt, self.params)

    G_hat = [None] * t

    bs = BitStream(munu)

    for i, h_i in enumerate(h):
      if h_i == 0:
        logging.debug(f"seeds[{i}]:\n%s", [int(v) for v in seeds[i]])

        seed = seeds[i]

        while True:
          pub_seed_A_tilde, pub_seed_B_tilde, seed = XOF(seed, [param.st_seed_bytes, param.st_seed_bytes, param.st_seed_bytes])

          logging.debug(f"sigma_A_tilde[{i}]:\n%s", [int(i) for i in pub_seed_A_tilde])
          A_tilde = rnd_inv_matrix(pub_seed_A_tilde, GFq, m)

          logging.debug(f"sigma_B_tilde[{i}]:\n%s", [int(i) for i in pub_seed_B_tilde])
          B_tilde = rnd_inv_matrix(pub_seed_B_tilde, GFq, n)

          logging.debug(f"A_tilde[{i}]:\n%s", A_tilde)
          logging.debug(f"B_tilde[{i}]:\n%s", B_tilde)

          G_hat[i] = pi(A_tilde, B_tilde, G[0])

          logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])

          G_hat[i] = G_hat[i].echelon_form()

          logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])

          # check if we got systematic form
          if sum([G_hat[i][j,j] for j in range(k)]) == k:
            break

      else:
        mu = matrix(GFq, m, m, [GFq(bs.read(GF_BITS)) for _ in range(m*m)])
        bs.finalize()

        nu = matrix(GFq, n, n, [GFq(bs.read(GF_BITS)) for _ in range(n*n)])
        bs.finalize()

        logging.debug(f"mu[{i}]:\n%s", mu)
        logging.debug(f"nu[{i}]:\n%s", nu)

        G_hat[i] = pi(mu, nu, G[h_i])

        logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])

        G_hat[i] = G_hat[i].echelon_form()

        logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])

    comp = bytes()

    for G_hat_i in G_hat:
      bs = BitStream()
      compress(G_hat_i[:,k:], bs)

      comp += bs.data

    check = H(self.params)(comp + msg)

    if not digest == check:
      raise BadSignatureError("Signature verification failed!")

    return msg

class MEDS(MEDSbase):
  def __init__(self, param, rng = secrets.token_bytes):
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
    randombytes.randombytes_init(bytes([0]*48), None, 256)

    meds = MEDS(args.parset, randombytes.randombytes)

    print(f"parameter set: {meds.params.name}\n")

    start_time = time.time()
    meds.crypto_sign_keypair()
    pk_time = time.time() - start_time

    print("pk:  {0} bytes".format(len(meds.pk)))

    assert(len(meds.pk) == meds.params.pk_size)

    msg = b"Test"

    start_time = time.time()
    sm = meds.crypto_sign(msg)
    sig_time = time.time() - start_time

    print("sig: {} bytes".format(len(sm)))

    assert(len(sm)-len(msg) == meds.params.sig_size)

    ## break signed message for testing:
    #sm = bytearray(sm); sm[54] = 8

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

