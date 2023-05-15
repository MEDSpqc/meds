#ifndef MEDS_H
#define MEDS_H

#include "params.h"
#include "clog.h"


#define MEDS_t_mask ((1 << CLOG2(MEDS_t))-1)
#define MEDS_s_mask ((1 << CLOG2(MEDS_s))-1)

#ifdef HAVE_SEED_TREE
  #define MEDS_MAX_PATH_LEN ((1 << CLOG2(MEDS_w)) + MEDS_w * (CLOG2(MEDS_t) - CLOG2(MEDS_w) - 1))
#else
  #define MEDS_MAX_PATH_LEN (MEDS_t-MEDS_w)
#endif


//#define MEDS_SIG_BYTES ((MEDS_MAX_PATH_LEN*MEDS_sec_bytes + (MEDS_w*MEDS_m*MEDS_m + MEDS_w*MEDS_n*MEDS_n) * sizeof(GFq_t)) + MEDS_sec_bytes)

#define SK_BYTES (MEDS_sec_bytes + ((MEDS_s-1)*MEDS_m*MEDS_m + (MEDS_s-1)*MEDS_n*MEDS_n + MEDS_k*(MEDS_m*MEDS_n-MEDS_k)) * sizeof(GFq_t))

//#define PK_BYTES (MEDS_sec_bytes + ((MEDS_s-1)*MEDS_k*(MEDS_m*MEDS_n-MEDS_k) * sizeof(GFq_t)))


int keygen(uint8_t *seed, int seed_len, uint8_t *sk, int sk_len, uint8_t *pk, int pk_len);

int sign(uint8_t *seed, int seed_len, uint8_t *sk, int sk_len, const char *msg, int msg_len, uint8_t *sig, int sig_len);

char* verify(uint8_t *pk, int pk_len, char *msg, int msg_len, uint8_t *sig, int sig_len);

#endif

