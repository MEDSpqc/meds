import sys
import argparse
import logging
import os.path as path
import binascii

import MEDS
import params

import randombytes

def remove_prefix(text, prefix):
  assert text.startswith(prefix), f"bad prefix - '{text[:16]}' does not start with {prefix}"
  return text[len(prefix):]


if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  
  parser.add_argument("-v", "--verbous", action='store_true', help = "verbous debugging output")
  parser.add_argument("-c", "--count", type=int, help = "number of KAT rounds", default=100)
  parser.add_argument('parset', nargs='?', help = "parameter set", default=params.params[-1].name)
  
  args = parser.parse_args()
  
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

  args.parset = params.par_list[args.parset]


  sys.stdout.write(f"PQCsignKAT_{args.parset.sk_size}\n\n")

 
  req_file = f"PQCsignKAT_{args.parset.sk_size}.req"
   
  with open(req_file, "w") as req:
    randombytes.randombytes_init(list(range(48)), None, 256)

    for count in range(args.count):
      req.write(f"count = {count}\n")

      seed = randombytes.randombytes(48)

      req.write("seed = " + "".join([f"{v:02X}" for v in seed]) + "\n")

      mlen = 33*(count+1)

      req.write(f"mlen = {mlen}\n")

      msg = randombytes.randombytes(mlen)

      req.write("msg = " + "".join([f"{v:02X}" for v in msg]) + "\n")

      req.write("pk =\n")
      req.write("sk =\n")
      req.write("smlen =\n")
      req.write("sm =\n\n")


  rsp_file = f"PQCsignKAT_{args.parset.sk_size}.rsp"



  with open(req_file, "r") as req, open(rsp_file, "w") as rsp:
    rsp.write(f"# {args.parset.name}\n\n")

    while (line := req.readline().strip()) != "":
      count = int(remove_prefix(line, "count = "))
      rsp.write(f"count = {count}\n")

      sys.stdout.write(f"count = {count}\n")

      line = req.readline().strip()
      seed = binascii.unhexlify(remove_prefix(line, "seed = "))
      rsp.write("seed = " + "".join([f"{v:02X}" for v in seed]) + "\n")

      randombytes.randombytes_init(seed, None, 256)

      line = req.readline().strip()
      mlen = int(remove_prefix(line, "mlen = "))
      rsp.write(f"mlen = {mlen}\n")

      line = req.readline().strip()
      msg = binascii.unhexlify(remove_prefix(line, "msg = "))
      rsp.write("msg = " + "".join([f"{v:02X}" for v in msg]) + "\n")

      meds = MEDS.MEDS(args.parset, randombytes)

      meds.crypto_sign_keypair()

      rsp.write("pk = " + "".join([f"{v:02X}" for v in meds.pk]) + "\n")
      rsp.write("sk = " + "".join([f"{v:02X}" for v in meds.sk]) + "\n")

      sm = meds.crypto_sign(msg)

      rsp.write(f"smlen = {len(sm)}\n")
      rsp.write("sm = " + "".join([f"{v:02X}" for v in sm]) + "\n\n")

      msg1 = meds.crypto_sign_open(sm)

      assert len(msg) == len(msg1), f"crypto_sign_open returned bad 'mlen': Got {len(msg1)}, expected {len(msg)}\n"
      assert msg == msg1, "crypto_sign_open returned bad 'm' value\n"

      remove_prefix(req.readline().strip(), "pk =")
      remove_prefix(req.readline().strip(), "sk =")
      remove_prefix(req.readline().strip(), "smlen =")
      remove_prefix(req.readline().strip(), "sm =")
      line = req.readline()

