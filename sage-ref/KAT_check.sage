#!/usr/bin/env sage

import sys, random, string, base64, secrets, binascii
import logging

import MEDS
import params

import argparse

import os.path as path

import randombytes

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  
  parser.add_argument("-v", "--verbous", action='store_true', help = "verbous debugging output")
  parser.add_argument('katfile', help = "KAT input file")
  parser.add_argument('parset', nargs='?', help = "parameter set", default=params.params[-1])
  
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
 

  with open(args.katfile, "r") as rsp:
    for line in rsp:
  
      if line.startswith("#"):
        continue
  
      if line.startswith("count "):
        print(line.strip())
  
      if line.startswith("seed "):
        seed = binascii.unhexlify(line[7:].strip())

        randombytes.randombytes_init(seed, None, 256)

        meds = MEDS.MEDS(args.parset, randombytes)

        meds.crypto_sign_keypair()

      if line.startswith("pk "):
        pk = binascii.unhexlify(line[5:].strip())

        assert pk == meds.pk, f"\n{binascii.hexlify(pk).decode()}\n{binascii.hexlify(meds.pk).decode()}"

      if line.startswith("sk "):
        sk = binascii.unhexlify(line[5:].strip())

        assert sk == meds.sk, f"\n{binascii.hexlify(sk).decode()}\n{binascii.hexlify(meds.sk).decode()}"

      if line.startswith("msg "):
        msg = binascii.unhexlify(line[6:].strip())

      if line.startswith("sm "):
        sm = binascii.unhexlify(line[5:].strip())

        tmp = meds.crypto_sign(msg)

        assert tmp == sm, f"\n{binascii.hexlify(tmp).decode()}\n{binascii.hexlify(sm).decode()}"

        tmp = meds.crypto_sign_open(sm)

        assert tmp == msg, f"\n{binascii.hexlify(tmp).decode()}\n{binascii.hexlify(msg).decode()}"

        print("ok\n")

