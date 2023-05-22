#!/usr/bin/python

import sys, argparse

import params

parser = argparse.ArgumentParser()

parser.add_argument('parset', nargs='?', help = "parameter set", default=params.params[-1])

args = parser.parse_args()
 
par_set = params.par_list[args.parset]

print(f"""#ifndef API_H
#define API_H

#define CRYPTO_SECRETKEYBYTES {par_set.sk_size}
#define CRYPTO_PUBLICKEYBYTES {par_set.pk_size}
#define CRYPTO_BYTES {par_set.sig_size}

#define CRYPTO_ALGNAME "{par_set.name}"

int crypto_sign_keypair(
    unsigned char *pk,
    unsigned char *sk
  );

int crypto_sign(
    unsigned char *sm, unsigned long long *smlen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *sk
  );

int crypto_sign_open(
    unsigned char *m, unsigned long long *mlen,
    const unsigned char *sm, unsigned long long smlen,
    const unsigned char *pk
  );

#endif
""")

