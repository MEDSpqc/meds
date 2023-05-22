#!/bin/bash

for p in `./params.py -l`; do
  v=ref

  mkdir -p MEDS/$v/$p

  cp *.c MEDS/$v/$p/
  cp *.h MEDS/$v/$p/

  rm MEDS/$v/$p/randombytes.*
  rm MEDS/$v/$p/KAT_test.c

  cp -r NIST MEDS/$v/$p/

  ln -s NIST/PQCgenKAT_sign.c MEDS/$v/$p/
  ln -s NIST/rng.c MEDS/$v/$p/randombytes.c
  ln -s NIST/rng.h MEDS/$v/$p/randombytes.h

  ./params.py -p $p > MEDS/$v/$p/params.h

  ./params.py -a $p > MEDS/$v/$p/api.h

  cp NIST.mk MEDS/$v/$p/Makefile
done

ln -s ref MEDS/opt

