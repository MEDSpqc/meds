#!/usr/bin/env sage

import secrets, time, sys, logging, binascii

import os.path as path

from Crypto.Hash import SHAKE256

from seedtree import SeedTree

from util import *

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

    A_inv = [None] * s
    B_inv = [None] * s
    T_inv = [None] * s
    G = [None] * s

    delta = self.randombytes(param.sec_seed_bytes)

    sigma_G0, sigma = XOF(delta, [param.pub_seed_bytes, param.sec_seed_bytes])

    logging.debug(f"sigma:\n%s", [int(i) for i in sigma])
    logging.debug(f"sigma_G0:\n%s", [int(i) for i in sigma_G0])

    G[0] = ExpandSystMat(sigma_G0, GFq, k, m, n)

    logging.debug(f"G[0]:\n%s", G[0])

    for i in range(1, s):
      while True: # repeat until A[i] and B[i] are invertible and G[i] has systematic form
        sigma_T, sigma = XOF(sigma, [param.sec_seed_bytes, param.sec_seed_bytes])

        Ti = ExpandInvMat(sigma_T, GFq, k)
        logging.debug(f"Ti:\n%s", Ti)

        G0prime = Ti * G[0]

        logging.debug(f"G0prime:\n%s", G0prime)


        check_A_i, check_B_i = solve_symb([matrix(GFq, m, n, G0prime.rows()[j]) for j in range(2)])

        if check_A_i == None and check_B_i == None:
          logging.debug("no sol")
          continue  # try agian for this index

        if not check_B_i.is_invertible():
          logging.debug("no B")
          continue  # try agian for this index

        if not check_A_i.is_invertible():
          logging.debug("no A_inv")
          continue  # try agian for this index

        A_i, A_inv[i] = check_A_i, check_A_i.inverse()
        B_i, B_inv[i] = check_B_i.inverse(), check_B_i

        logging.debug(f"A[{i}]:\n%s", A_i)
        logging.debug(f"A_inv[{i}]:\n%s", A_inv[i])
        logging.debug(f"B[{i}]:\n%s", B_i)
        logging.debug(f"B_inv[{i}]:\n%s", B_inv[i])

        G[i] = pi(A_i, B_i, G[0])

        logging.debug(f"G[{i}]:\n%s", G[i])

        T_inv[i] = G[i].submatrix(0,0,k,k)

        if T_inv[i].is_invertible():
          G[i] = T_inv[i].inverse() * G[i]

          logging.debug(f"G[{i}]:\n%s", G[i])
          logging.debug(f"T_inv[{i}]:\n%s", T_inv[i])

          # G[i] is in systematic form; break while loop
          break

        logging.debug(f"redo G[{i}]")


    # pk
    pk = sigma_G0

    for Gi in G[1:]:
      pk += CompressG(Gi, k, m, n)

    logging.debug(f"sigma_G0 (pk):\n%s", [int(i) for i in pk[:param.pub_seed_bytes]])
    logging.debug(f"G[1:] (pk):\n%s", [int(j) for j in pk[param.pub_seed_bytes:]])

    logging.debug(f"pk:\n0x%s", binascii.hexlify(pk).decode())


    # sk
    sk = delta + sigma_G0

    for A_inv_i in A_inv[1:]:
      sk += Compress(A_inv_i)

    for B_inv_i in B_inv[1:]:
      sk += Compress(B_inv_i)

    for T_inv_i in T_inv[1:]:
      sk += Compress(T_inv_i)

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

    GF_BITS = ceil(log(q, 2))


    f_sk = param.sec_seed_bytes

    sigma_G0  = sk[f_sk : f_sk + param.pub_seed_bytes]
    f_sk += param.pub_seed_bytes

    G_0 = ExpandSystMat(sigma_G0, GFq, k, m, n)


    A_inv = [None]*s
    B_inv = [None]*s
    T_inv = [None]*s

    l_Fq_nn = ceil((GF_BITS * n * n) / 8)
    l_Fq_mm = ceil((GF_BITS * m * m) / 8)
    l_Fq_kk = ceil((GF_BITS * k * k) / 8)

    for i in range(1, s):
      A_inv[i] = Decompress(sk[f_sk : f_sk + l_Fq_mm], GFq, m, m)
      f_sk += l_Fq_mm

      logging.debug(f"A_inv[{i}]:\n%s", A_inv[i])

    for i in range(1, s):
      B_inv[i] = Decompress(sk[f_sk : f_sk + l_Fq_nn], GFq, n, n)
      f_sk += l_Fq_nn

      logging.debug(f"B_inv[{i}]:\n%s", B_inv[i])

    for i in range(1, s):
      T_inv[i] = Decompress(sk[f_sk : f_sk + l_Fq_kk], GFq, k, k)
      f_sk += l_Fq_kk

      logging.debug(f"T_inv[{i}]:\n%s", T_inv[i])


    logging.debug(f"G_0:\n%s", G_0)


    delta = self.randombytes(param.sec_seed_bytes)

    logging.debug(f"delta:\n%s", [int(i) for i in delta])

    rho, alpha = XOF(delta, [param.st_seed_bytes, param.st_salt_bytes])

    sigma = SeedTree(t, rho, G(self.params, alpha))

    for i in range(t):
      logging.debug(f"sigma[{i}]:\n0x%s", binascii.hexlify(sigma[i]).decode())

    A_tilde = [None]*t
    B_tilde = [None]*t
    G_tilde = [None]*t
    M_tilde = [None]*t

    for i in range(t):
      while True:

        seed_buf = alpha + sigma[i] + i.to_bytes(4, "little")

        logging.debug(f"seed_buf:\n0x%s", binascii.hexlify(seed_buf).decode())

        sigma_M_tilde, sigma[i] = XOF(seed_buf, [param.pub_seed_bytes, param.st_seed_bytes])


        logging.debug(f"sigma_M_tilde[{i}]:\n0x%s", binascii.hexlify(sigma_M_tilde).decode())

        M_tilde[i] = matrix(GFq, 2, k, ExpandFqs(sigma_M_tilde, 2*k, GFq))

        logging.debug(f"M_tilde[{i}]:\n%s", M_tilde[i])

        if (M_tilde[i].rank() != M_tilde[i].nrows()):
          continue


        G0_prime = M_tilde[i] * G_0

        logging.debug(f"G0_prime:\n%s", G0_prime)

        A_tilde[i], B_tilde_inv = solve_symb([matrix(GFq, m, n, G0_prime.rows()[j]) for j in range(2)])

        if A_tilde[i] == None and B_tilde_inv == None:
          logging.debug("no sol")
          continue  # try agian for this index

        if not B_tilde_inv.is_invertible():
          logging.debug("no B_tilde")
          continue  # try agian for this index

        if not A_tilde[i].is_invertible():
          logging.debug("no A_tilde_inv")
          continue  # try agian for this index


        B_tilde[i] = B_tilde_inv.inverse()


        logging.debug(f"A_tilde[{i}]:\n%s", A_tilde[i])
        logging.debug(f"B_tilde[{i}]:\n%s", B_tilde[i])

        G_tilde[i] = pi(A_tilde[i], B_tilde[i], G_0)

        logging.debug(f"G_tilde[{i}]:\n%s", G_tilde[i])

        G_tilde[i] = SF(G_tilde[i])

        logging.debug(f"G_tilde[{i}]:\n%s", G_tilde[i])

        # check if we got systematic form
        if G_tilde[i] == None:
          # if no systematic form loop to try again for this index
          logging.debug(f"redo G[{i}]")
          continue

        # G_tilde[i] is in systematic form; breal while loop
        break


    if type(msg) == str:
      msg = msg.encode('utf8')

    comp = bytes()

    for G_tilde_i in G_tilde:
      comp += Compress(G_tilde_i[:,k:])

    d = H(self.params)(comp + msg)

    logging.debug(f"digest:\n%s", [int(i) for i in d])

    h = PaseHash(d, self.params)

    logging.debug(f"h:\n{h}")

    ret = bytes()

    for i in range(t):
      if h[i] > 0:
        kappa = M_tilde[i] * T_inv[h[i]]

        logging.debug(f"kappa[{i}]:\n%s", kappa)

        ret += Compress(kappa)

    p = SeedTreeToPath(h, rho, alpha, self.params)

    ret += p
    ret += d
    ret += alpha
    ret += msg

    logging.debug(f"sm:\n0x%s", binascii.hexlify(ret).decode())

    return ret

  def crypto_sign_open(self, pk, ms):
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
    logging.debug(f"sm:\n0x%s", binascii.hexlify(ms).decode())

    I_k = matrix.identity(ring=GFq, n=k)

    sigma_G0 = pk[:param.pub_seed_bytes]

    G = [ExpandSystMat(sigma_G0, GFq, k, m, n)]

    logging.debug(f"G[0]:\n%s", G[0])


    l_Gi = ceil((k-2)*(m*n-k) * GF_BITS / 8)
    f_pk = param.pub_seed_bytes

    for i in range(1, s):
      G.append(DecompressG(pk[f_pk : f_pk + l_Gi], GFq, k, m, n))
      f_pk += l_Gi

      logging.debug(f"G[{i}]:\n%s", G[i])


    kappa = ms[:ceil(2*k * GF_BITS / 8)*w]
    ms = ms[ceil(2*k * GF_BITS / 8)*w:]
    path = ms[:param.seed_tree_cost]; ms = ms[param.seed_tree_cost:]
    d = ms[:param.digest_bytes]; ms = ms[param.digest_bytes:]
    alpha = ms[:param.st_salt_bytes]; msg = ms[param.st_salt_bytes:]

    logging.debug(f"kappa:\n0x%s", binascii.hexlify(kappa).decode())
    logging.debug(f"path:\n0x%s", binascii.hexlify(path).decode())
    logging.debug(f"digest:\n0x%s", binascii.hexlify(d).decode())
    logging.debug(f"alpha:\n0x%s", binascii.hexlify(alpha).decode())

    h = PaseHash(d, self.params)

    sigma = PathToSeedTree(h, path, alpha, self.params)

    G_hat = [None] * t

    l_Fq_kk = ceil((GF_BITS * 2 * k) / 8)

    f_ms = 0

    for i in range(t):
      if h[i] > 0:
        kappa_i = Decompress(kappa[f_ms : f_ms + l_Fq_kk], GFq, 2, k)
        f_ms += l_Fq_kk

        logging.debug(f"kappa[{i}]:\n%s", kappa_i);

        G0_prime = kappa_i * G[h[i]]

        logging.debug(f"G0_prime[{i}]:\n%s", G0_prime);

        A_hat_i, B_hat_i_inv = solve_symb([matrix(GFq, m, n, G0_prime.rows()[j]) for j in range(G0_prime.nrows())])

        if A_hat_i == None and B_hat_inv == None:
          logging.debug("no sol")
          continue  # try agian for this index

        if not B_hat_i_inv.is_invertible():
          logging.debug("no B_hat")
          continue  # try agian for this index

        if not A_hat_i.is_invertible():
          logging.debug("no A_hat_inv")
          continue  # try agian for this index


        B_hat_i = B_hat_i_inv.inverse()


        logging.debug(f"A_hat[{i}]:\n%s", A_hat_i)
        logging.debug(f"B_hat[{i}]:\n%s", B_hat_i)


        G_hat[i] = pi(A_hat_i, B_hat_i, G[h[i]])

        logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])

        G_hat[i] = SF(G_hat[i])

        if G_hat[i] == None:
          raise BadSignatureError("Signature verification failed!")

        logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])

      else:
        logging.debug(f"seeds[{i}]:\n%s", [int(v) for v in sigma[i]])

        while True:

          seed_buf = alpha + sigma[i] + i.to_bytes(4, "little")

          logging.debug(f"seed_buf:\n0x%s", binascii.hexlify(seed_buf).decode())

          sigma_M_hat_i, sigma[i] = XOF(seed_buf, [param.pub_seed_bytes, param.st_seed_bytes])

          logging.debug(f"sigma_M_hat[{i}]:\n0x%s", binascii.hexlify(sigma_M_hat_i).decode())

          M_hat_i = matrix(GFq, 2, k, ExpandFqs(sigma_M_hat_i, 2*k, GFq))

          logging.debug(f"M_hat[{i}]:\n%s", M_hat_i)

          if (M_hat_i.rank() != M_hat_i.nrows()):
            continue

          G0_prime = M_hat_i * G[0]

          logging.debug(f"G0_prime[{i}]:\n%s", G0_prime)

          A_hat_i, B_hat_i_inv = solve_symb([matrix(GFq, m, n, G0_prime.rows()[j]) for j in range(2)])

          if A_hat_i == None and B_hat_inv == None:
            logging.debug("no sol")
            continue  # try agian for this index

          if not B_hat_i_inv.is_invertible():
            logging.debug("no B_hat")
            continue  # try agian for this index

          if not A_hat_i.is_invertible():
            logging.debug("no A_hat_inv")
            continue  # try agian for this index

          B_hat_i = B_hat_i_inv.inverse()


          logging.debug(f"A_hat[{i}]:\n%s", A_hat_i)
          logging.debug(f"B_hat[{i}]:\n%s", B_hat_i)


          G_hat[i] = pi(A_hat_i, B_hat_i, G[0])

          logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])

          G_hat[i] = SF(G_hat[i])

          logging.debug(f"G_hat[{i}]:\n%s", G_hat[i])

          # check if we got systematic form
          if G_hat[i] == None:
            # if no systematic form loop to try again for this index
            logging.debug(f"redo G[{i}]")
            continue

          # G_hat[i] is in systematic form; breal while loop
          break


    comp = bytes()

    for G_hat_i in G_hat:
      comp += Compress(G_hat_i[:,k:])

    d_prime = H(self.params)(comp + msg)

    if not d == d_prime:
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

  from randombytes import *

  # Test and benchmark:
  try:
    entropy_input = [0]*48
    randombytes_init(bytes(entropy_input), None, 256)

    meds = MEDS(args.parset, randombytes)

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

