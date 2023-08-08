#!/usr/bin/env sage

import sys, random, string, base64, secrets, binascii
import logging

from multiprocessing import Pool

import meds
import params

import argparse

import os.path as path

from randombytes import *


class Task:
  def init(self):
    self.count = None
    self.seed = None
    self.pk = None
    self.sk = None
    self.msg = None
    self.sm = None


def check(t):
  ret = str(t.count)

  try:
    randombytes_init(t.seed, None, 256)

    MEDS = meds.MEDS(args.parset, randombytes)

    MEDS.crypto_sign_keypair()

    assert t.pk == MEDS.pk, f"\n{binascii.hexlify(pk).decode()}\n{binascii.hexlify(MEDS.pk).decode()}"

    assert t.sk == MEDS.sk, f"\n{binascii.hexlify(sk).decode()}\n{binascii.hexlify(MEDS.sk).decode()}"

    tmp = MEDS.crypto_sign(t.msg)

    assert tmp == t.sm, f"\n{binascii.hexlify(tmp).decode()}\n{binascii.hexlify(sm).decode()}"

    tmp = MEDS.crypto_sign_open(t.sm)

    assert tmp == t.msg, f"\n{binascii.hexlify(tmp).decode()}\n{binascii.hexlify(msg).decode()}"

    ret += " ok"

  except NameError:
    ret += " ERROR!"

  return ret

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

  t = None

  with Pool() as pool:
    procs = []

    with open(args.katfile, "r") as rsp:
      for line in rsp:

        if line.startswith("#"):
          continue

        if line.startswith("count "):
          t = Task()

          t.count = line.strip()

        if line.startswith("seed "):
          t.seed = binascii.unhexlify(line[7:].strip())

        if line.startswith("pk "):
          t.pk = binascii.unhexlify(line[5:].strip())


        if line.startswith("sk "):
          t.sk = binascii.unhexlify(line[5:].strip())


        if line.startswith("msg "):
          t.msg = binascii.unhexlify(line[6:].strip())

        if line.startswith("sm "):
          t.sm = binascii.unhexlify(line[5:].strip())

          procs.append(pool.apply_async(check, [t]))

    for proc in procs:
      value = proc.get()

      print(value)

    pool.close()

    pool.join()


