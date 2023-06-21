#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"

#include "fips202.h"

#include "params.h"

#include "api.h"
#include "randombytes.h"

#include "meds.h"

#define solve solve_opt

#include "seed.h"
#include "util.h"
#include "bitstream.h"

#include "matrixmod.h"

#define CEILING(x,y) (((x) + (y) - 1) / (y))


int crypto_sign_keypair(
    unsigned char *pk,
    unsigned char *sk
  )
{
  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);

  pmod_mat_t G_data[MEDS_k * MEDS_m * MEDS_n * MEDS_s];
  pmod_mat_t *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
    G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];

  uint8_t sigma_G0[MEDS_pub_seed_bytes];
  uint8_t sigma[MEDS_sec_seed_bytes];

  XOF((uint8_t*[]){sigma_G0, sigma},
       (size_t[]){MEDS_pub_seed_bytes, MEDS_sec_seed_bytes},
       delta, MEDS_sec_seed_bytes,
       2);

  LOG_VEC(sigma, MEDS_sec_seed_bytes);
  LOG_VEC_FMT(sigma_G0, MEDS_pub_seed_bytes, "sigma_G0");


  rnd_sys_mat(G[0], MEDS_k, MEDS_m*MEDS_n, sigma_G0, MEDS_pub_seed_bytes);

  LOG_MAT(G[0], MEDS_k, MEDS_m*MEDS_n);

  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_n * MEDS_n];
  pmod_mat_t T_inv_data[MEDS_s * MEDS_k * MEDS_k];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];
  pmod_mat_t *T_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
  {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
    T_inv[i] = &T_inv_data[i * MEDS_k * MEDS_k];
  }

  for (int i = 1; i < MEDS_s; i++)
  {
    pmod_mat_t A[MEDS_m * MEDS_m] = {0};
    pmod_mat_t B[MEDS_n * MEDS_n] = {0};

    while (1 == 1) // redo generation for this index until success
    {
      uint8_t sigma_Ti[MEDS_sec_seed_bytes];

      XOF((uint8_t*[]){sigma_Ti, sigma},
          (size_t[]){MEDS_sec_seed_bytes, MEDS_sec_seed_bytes},
          sigma, MEDS_sec_seed_bytes,
          2);

      pmod_mat_t Ti[MEDS_k * MEDS_k];

      rnd_inv_matrix(Ti, MEDS_k, MEDS_k, sigma_Ti, MEDS_sec_seed_bytes);

      LOG_MAT(Ti, MEDS_k, MEDS_k);

      pmod_mat_t G0prime[MEDS_k * MEDS_m * MEDS_n];

      pmod_mat_mul(G0prime, MEDS_k, MEDS_m * MEDS_n,
          Ti, MEDS_k, MEDS_k,
          G[0], MEDS_k, MEDS_m * MEDS_n);

      LOG_MAT(G0prime, MEDS_k, MEDS_m * MEDS_n);


      if (solve(A, B_inv[i], G0prime, false) < 0)
      {
        LOG("no sol");
        continue;
      }

      if (pmod_mat_inv(B, B_inv[i], MEDS_n, MEDS_n) < 0)
      {
        LOG("no inv B");
        continue;
      }

      if (pmod_mat_inv(A_inv[i], A, MEDS_m, MEDS_m) < 0)
      {
        LOG("no inv A_inv");
        continue;
      }

      LOG_MAT_FMT(A, MEDS_m, MEDS_m, "A[%i]", i);
      LOG_MAT_FMT(A_inv[i], MEDS_m, MEDS_m, "A_inv[%i]", i);
      LOG_MAT_FMT(B, MEDS_n, MEDS_n, "B[%i]", i);
      LOG_MAT_FMT(B_inv[i], MEDS_n, MEDS_n, "B_inv[%i]", i);


      pi(G[i], A, B, G[0]);

      LOG_MAT_FMT(G[i], MEDS_k, MEDS_m*MEDS_n, "G[%i]", i);

      for (int r = 0; r < MEDS_k; r++)
        memcpy(&T_inv[i][r * MEDS_k], &G[i][r * MEDS_m * MEDS_n], MEDS_k * sizeof(GFq_t));


      if (SF(G[i], G[i]) != 0)
      {
        LOG("redo G[%i]", i);
        continue; // Not systematic; try again for index i.
      }

      LOG_MAT_FMT(G[i], MEDS_k, MEDS_m*MEDS_n, "G[%i]", i);
      LOG_MAT_FMT(T_inv[i], MEDS_k, MEDS_k, "T_inv[%i]", i);

      // successfull generated G[s]; break out of while loop
      break;
    }
  }


  // copy pk data
  {
    uint8_t *tmp_pk = pk;

    memcpy(tmp_pk, sigma_G0, MEDS_pub_seed_bytes);
    LOG_VEC(tmp_pk, MEDS_pub_seed_bytes, "sigma_G0 (pk)");
    tmp_pk += MEDS_pub_seed_bytes;

    bitstream_t bs;

    bs_init(&bs, tmp_pk, MEDS_PK_BYTES - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int r = 2; r < MEDS_k; r++)
        for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
          bs_write(&bs, G[si][r*MEDS_m*MEDS_n + j], GFq_bits);

      bs_finalize(&bs);
    }

    LOG_VEC(tmp_pk, MEDS_PK_BYTES - MEDS_pub_seed_bytes, "G[1:] (pk)");
    tmp_pk += MEDS_PK_BYTES - MEDS_pub_seed_bytes;

    LOG_HEX(pk, MEDS_PK_BYTES);

    if (MEDS_PK_BYTES != MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0))
    {
      fprintf(stderr, "ERROR: MEDS_PK_BYTES and actual pk size do not match! %i vs %i\n", MEDS_PK_BYTES, MEDS_pub_seed_bytes + bs.byte_pos+(bs.bit_pos > 0 ? 1 : 0));
      fprintf(stderr, "%i %i\n", MEDS_pub_seed_bytes + bs.byte_pos, MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0));
      return -1;
    }
  }

  // copy sk data
  {
    memcpy(sk, delta, MEDS_sec_seed_bytes);
    memcpy(sk + MEDS_sec_seed_bytes, sigma_G0, MEDS_pub_seed_bytes);

    bitstream_t bs;

    bs_init(&bs, sk + MEDS_sec_seed_bytes + MEDS_pub_seed_bytes, MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_m*MEDS_m; j++)
        bs_write(&bs, A_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_n*MEDS_n; j++)
        bs_write(&bs, B_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_k*MEDS_k; j++)
        bs_write(&bs, T_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    LOG_HEX(sk, MEDS_SK_BYTES);
  }

  return 0;
}

int crypto_sign(
    unsigned char *sm, unsigned long long *smlen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *sk
  )
{
  // skip secret seed
  sk += MEDS_sec_seed_bytes;

  pmod_mat_t G_0[MEDS_k * MEDS_m * MEDS_n];


  rnd_sys_mat(G_0, MEDS_k, MEDS_m*MEDS_n, sk, MEDS_pub_seed_bytes);

  sk += MEDS_pub_seed_bytes;


  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_n * MEDS_n];
  pmod_mat_t T_inv_data[MEDS_s * MEDS_k * MEDS_k];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];
  pmod_mat_t *T_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
  {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
    T_inv[i] = &T_inv_data[i * MEDS_k * MEDS_k];
  }
 
  // Load secret key matrices.
  {
    bitstream_t bs;

    bs_init(&bs, (uint8_t*)sk, MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_m*MEDS_m; j++)
        A_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_n*MEDS_n; j++)
        B_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_k*MEDS_k; j++)
        T_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    bs_finalize(&bs);
  }


  for (int i = 1; i < MEDS_s; i++)
    LOG_MAT_FMT(A_inv[i], MEDS_m, MEDS_m, "A_inv[%i]", i);

  for (int i = 1; i < MEDS_s; i++)
    LOG_MAT_FMT(B_inv[i], MEDS_n, MEDS_n, "B_inv[%i]", i);

  for (int i = 1; i < MEDS_s; i++)
    LOG_MAT_FMT(T_inv[i], MEDS_k, MEDS_k, "T_inv[%i]", i);

  LOG_MAT(G_0, MEDS_k, MEDS_m*MEDS_n);


  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);

  LOG_VEC(delta, MEDS_sec_seed_bytes);


  uint8_t stree[MEDS_st_seed_bytes * SEED_TREE_size] = {0};
  uint8_t alpha[MEDS_st_salt_bytes];

  uint8_t *rho = &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(0,0)];

  XOF((uint8_t*[]){rho, alpha},
      (size_t[]){MEDS_st_seed_bytes, MEDS_st_salt_bytes},
      delta, MEDS_sec_seed_bytes,
      2);

  t_hash(stree, alpha, 0, 0);

  uint8_t *sigma = &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(MEDS_seed_tree_height, 0)];

  for (int i = 0; i < MEDS_t; i++)
  {
     LOG_HEX_FMT((&sigma[i*MEDS_st_seed_bytes]), MEDS_st_seed_bytes, "sigma[%i]", i);
  }

  pmod_mat_t A_tilde_data[MEDS_t * MEDS_m * MEDS_m];
  pmod_mat_t B_tilde_data[MEDS_t * MEDS_n * MEDS_n];
  pmod_mat_t M_tilde_data[MEDS_t * 2 * MEDS_k];

  pmod_mat_t *A_tilde[MEDS_t];
  pmod_mat_t *B_tilde[MEDS_t];
  pmod_mat_t *M_tilde[MEDS_t];

  for (int i = 0; i < MEDS_t; i++)
  {
    A_tilde[i] = &A_tilde_data[i * MEDS_m * MEDS_m];
    B_tilde[i] = &B_tilde_data[i * MEDS_n * MEDS_n];
    M_tilde[i] = &M_tilde_data[i * 2 * MEDS_k];
  }

  keccak_state h_shake;
  shake256_init(&h_shake);


  uint8_t seed_buf[MEDS_st_salt_bytes + MEDS_st_seed_bytes + sizeof(uint32_t)] = {0};
  memcpy(seed_buf, alpha, MEDS_st_salt_bytes);

  uint8_t *addr_pos = seed_buf + MEDS_st_salt_bytes + MEDS_st_seed_bytes;


  for (int i = 0; i < MEDS_t; i++)
  {
    pmod_mat_t G_tilde_ti[MEDS_k * MEDS_m * MEDS_n];

    while (1 == 1)
    {
      uint8_t sigma_M_tilde_i[MEDS_pub_seed_bytes];


      for (int j = 0; j < 4; j++)
        addr_pos[j] = (i >> (j*8)) & 0xff;

      memcpy(seed_buf + MEDS_st_salt_bytes, &sigma[i*MEDS_st_seed_bytes], MEDS_st_seed_bytes);


      XOF((uint8_t*[]){sigma_M_tilde_i, &sigma[i*MEDS_st_seed_bytes]},
           (size_t[]){MEDS_pub_seed_bytes, MEDS_st_seed_bytes},
           seed_buf, MEDS_st_salt_bytes + MEDS_st_seed_bytes + sizeof(uint32_t),
           2);

      LOG_HEX_FMT(sigma_M_tilde_i, MEDS_pub_seed_bytes, "sigma_M_tilde[%i]", i);

      rnd_inv_matrix(M_tilde[i], 2, MEDS_k, sigma_M_tilde_i, MEDS_pub_seed_bytes);

      LOG_MAT_FMT(M_tilde[i], 2, MEDS_k, "M_tilde[%i]", i);


      pmod_mat_t G0_prime[2 * MEDS_m * MEDS_n];

      pmod_mat_mul(G0_prime, 2, MEDS_m * MEDS_n, M_tilde[i], 2, MEDS_k, G_0, MEDS_k, MEDS_m * MEDS_n);

      LOG_MAT(G0_prime, 2, MEDS_m * MEDS_n);


      pmod_mat_t A_tilde_inv[MEDS_m * MEDS_m];
      pmod_mat_t B_tilde_inv[MEDS_n * MEDS_n];

      if (solve(A_tilde[i], B_tilde_inv, G0_prime, false) < 0)
      {
        LOG("no sol");
        continue;
      }

      if (pmod_mat_inv(B_tilde[i], B_tilde_inv, MEDS_n, MEDS_n) < 0)
      {
        LOG("no B_tilde");
        continue;
      }

      if (pmod_mat_inv(A_tilde_inv, A_tilde[i], MEDS_m, MEDS_m) < 0)
      {
        LOG("no A_tilde_inv");
        continue;
      }

      LOG_MAT_FMT(A_tilde[i], MEDS_m, MEDS_m, "A_tilde[%i]", i);
      LOG_MAT_FMT(B_tilde[i], MEDS_n, MEDS_n, "B_tilde[%i]", i);


      pi(G_tilde_ti, A_tilde[i], B_tilde[i], G_0);


      LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, "G_tilde[%i]", i);

      if (SF(G_tilde_ti, G_tilde_ti) == 0)
        break;
    }

    LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, "G_tilde[%i]", i);

    bitstream_t bs;
    uint8_t bs_buf[CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8)];
    
    bs_init(&bs, bs_buf, CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8));

    for (int r = 0; r < MEDS_k; r++)
      for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
        bs_write(&bs, G_tilde_ti[r * MEDS_m*MEDS_n + j], GFq_bits);

    shake256_absorb(&h_shake, bs_buf, CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8));
  }

  shake256_absorb(&h_shake, (uint8_t*)m, mlen);

  shake256_finalize(&h_shake);

  uint8_t digest[MEDS_digest_bytes];

  shake256_squeeze(digest, MEDS_digest_bytes, &h_shake);

  LOG_VEC(digest, MEDS_digest_bytes);


  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);

  LOG_VEC(h, MEDS_t);


  bitstream_t bs;

  bs_init(&bs, sm, MEDS_w * CEILING(2*MEDS_k * GFq_bits, 8));

  uint8_t *path = sm + MEDS_w * CEILING(2*MEDS_k * GFq_bits, 8);

  t_hash(stree, alpha, 0, 0);

  stree_to_path(stree, h, path, alpha);

  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] > 0)
    {
      {
        pmod_mat_t kappa[2*MEDS_k];

        pmod_mat_mul(kappa, 2, MEDS_k, M_tilde[i], 2, MEDS_k, T_inv[h[i]], MEDS_k, MEDS_k);

        LOG_MAT_FMT(kappa, 2, MEDS_k, "kappa[%i]", i);

        for (int j = 0; j <2*MEDS_k; j++)
          bs_write(&bs, kappa[j], GFq_bits);

//        // Check if verifier has systematic system:
//        int64_t sum = 0;
//
//        // Compute bottom right value of A_tilde[i] * A_inv[h[i]].
//        for (int j = 0; j < MEDS_m; j++)
//          sum = (sum + A_tilde[i][(MEDS_m-1)*MEDS_m + j] * A_inv[h[i]][j * MEDS_m + MEDS_m - 1]) % MEDS_p;
//
//        if (sum == 0)
//        {
//          LOG("REDO A[-1,-1] == 0");
//          goto redo;
//        }
      }

      bs_finalize(&bs);
    }
  }

  memcpy(sm + MEDS_SIG_BYTES - MEDS_digest_bytes - MEDS_st_salt_bytes, digest, MEDS_digest_bytes);
  memcpy(sm + MEDS_SIG_BYTES - MEDS_st_salt_bytes, alpha, MEDS_st_salt_bytes);
  memcpy(sm + MEDS_SIG_BYTES, m, mlen);

  *smlen = MEDS_SIG_BYTES + mlen;

  LOG_HEX(sm, MEDS_SIG_BYTES + mlen);

  return 0;
}

int crypto_sign_open(
    unsigned char *m, unsigned long long *mlen,
    const unsigned char *sm, unsigned long long smlen,
    const unsigned char *pk
  )
{
  LOG_HEX(pk, MEDS_PK_BYTES);
  LOG_HEX(sm, smlen);

  pmod_mat_t G_data[MEDS_k*MEDS_m*MEDS_n * MEDS_s];
  pmod_mat_t *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
    G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];


  rnd_sys_mat(G[0], MEDS_k, MEDS_m*MEDS_n, pk, MEDS_pub_seed_bytes);

  {
    bitstream_t bs;

    bs_init(&bs, (uint8_t*)pk + MEDS_pub_seed_bytes, MEDS_PK_BYTES - MEDS_pub_seed_bytes);

    for (int i = 1; i < MEDS_s; i++)
    {
      for (int r = 0; r < MEDS_k; r++)
        for (int c = 0; c < MEDS_k; c++)
          if (r == c)
            pmod_mat_set_entry(G[i], MEDS_k, MEDS_m * MEDS_n, r, c, 1);
          else
            pmod_mat_set_entry(G[i], MEDS_k, MEDS_m * MEDS_n, r, c, 0);

      for (int r = 2; r < MEDS_k; r++)
        for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
          G[i][r*MEDS_m*MEDS_n + j] = bs_read(&bs, GFq_bits);

      for (int ii = 0; ii < MEDS_m; ii++)
        for (int j = 0; j < MEDS_n; j++)
          G[i][ii*MEDS_n + j] = ii == j ? 1 : 0;

      for (int ii = 0; ii < MEDS_m; ii++)
        for (int j = 0; j < MEDS_n; j++)
          G[i][MEDS_m*MEDS_n + ii*MEDS_n + j] = (ii+1) == j ? 1 : 0;

      bs_finalize(&bs);
    }
  }

  for (int i = 0; i < MEDS_s; i++)
    LOG_MAT_FMT(G[i], MEDS_k, MEDS_m*MEDS_n, "G[%i]", i);

  LOG_HEX_FMT(sm, MEDS_w * CEILING(2*MEDS_k * GFq_bits, 8), "kappa");
  LOG_HEX_FMT(sm + MEDS_w * CEILING(2*MEDS_k * GFq_bits, 8),
      MEDS_max_path_len * MEDS_st_seed_bytes, "path");

  uint8_t *digest = (uint8_t*)sm + (MEDS_SIG_BYTES - MEDS_digest_bytes - MEDS_st_salt_bytes);

  uint8_t *alpha = (uint8_t*)sm + (MEDS_SIG_BYTES - MEDS_st_salt_bytes);

  LOG_HEX(digest, MEDS_digest_bytes);
  LOG_HEX(alpha, MEDS_st_salt_bytes);

  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);


  bitstream_t bs;

  bs_init(&bs, (uint8_t*)sm, MEDS_w * CEILING(2*MEDS_k * GFq_bits, 8));

  uint8_t *path = (uint8_t*)sm + MEDS_w * CEILING(2*MEDS_k * GFq_bits, 8);

  uint8_t stree[MEDS_st_seed_bytes * SEED_TREE_size] = {0};

  path_to_stree(stree, h, path, alpha);

  uint8_t *sigma = &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(MEDS_seed_tree_height, 0)];

  pmod_mat_t G_hat_i[MEDS_k*MEDS_m*MEDS_n];

  pmod_mat_t kappa[2*MEDS_k];

  keccak_state shake;
  shake256_init(&shake);


  uint8_t seed_buf[MEDS_st_salt_bytes + MEDS_st_seed_bytes + sizeof(uint32_t)] = {0};
  memcpy(seed_buf, alpha, MEDS_st_salt_bytes);

  uint8_t *addr_pos = seed_buf + MEDS_st_salt_bytes + MEDS_st_seed_bytes;


  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] > 0)
    {
      for (int j = 0; j < 2*MEDS_k; j++)
        kappa[j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);

      LOG_MAT_FMT(kappa, 2, MEDS_k, "kappa[%i]", i);


      pmod_mat_t G0_prime[2 * MEDS_m * MEDS_n];

      pmod_mat_mul(G0_prime, 2, MEDS_m * MEDS_n, kappa, 2, MEDS_k, G[h[i]], MEDS_k, MEDS_m * MEDS_n);

      LOG_MAT_FMT(G0_prime, 2, MEDS_m * MEDS_n, "G0_prime[%i]", i);


      pmod_mat_t A_hat[MEDS_m * MEDS_m];
      pmod_mat_t B_hat[MEDS_n * MEDS_n];

      pmod_mat_t A_hat_inv[MEDS_m * MEDS_m];
      pmod_mat_t B_hat_inv[MEDS_n * MEDS_n];

      if (solve(A_hat, B_hat_inv, G0_prime, true) < 0)
      {
        LOG("crypto_sign_open - no sol");
        printf("no sol\n");
        return -1;
      }

      if (pmod_mat_inv(B_hat, B_hat_inv, MEDS_n, MEDS_n) < 0)
      {
        LOG("no B_hat");
        return -1;
      }

      if (pmod_mat_inv(A_hat_inv, A_hat, MEDS_m, MEDS_m) < 0)
      {
        LOG("no A_hat_inv");
        return -1;
      }

      LOG_MAT_FMT(A_hat, MEDS_m, MEDS_m, "A_hat[%i]", i);
      LOG_MAT_FMT(B_hat, MEDS_n, MEDS_n, "B_hat[%i]", i);


      pi(G_hat_i, A_hat, B_hat, G[h[i]]);


      LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);

      if (SF(G_hat_i, G_hat_i) < 0)
      {
        fprintf(stderr, "Signature verification failed!\n");

        return -1;
      }

      LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);
    }
    else
    {
      while (1 == 1)
      {
        LOG_VEC_FMT(&sigma[i*MEDS_st_seed_bytes], MEDS_st_seed_bytes, "seeds[%i]", i);

        uint8_t sigma_M_hat_i[MEDS_pub_seed_bytes];


        for (int j = 0; j < 4; j++)
          addr_pos[j] = (i >> (j*8)) & 0xff;

        memcpy(seed_buf + MEDS_st_salt_bytes, &sigma[i*MEDS_st_seed_bytes], MEDS_st_seed_bytes);


        XOF((uint8_t*[]){sigma_M_hat_i, &sigma[i*MEDS_st_seed_bytes]},
            (size_t[]){MEDS_pub_seed_bytes, MEDS_st_seed_bytes},
            seed_buf, MEDS_st_salt_bytes + MEDS_st_seed_bytes + sizeof(uint32_t),
            2);

        pmod_mat_t M_hat_i[2*MEDS_k];

        LOG_HEX_FMT(sigma_M_hat_i, MEDS_pub_seed_bytes, "sigma_M_hat[%i]", i);

        rnd_inv_matrix(M_hat_i, 2, MEDS_k, sigma_M_hat_i, MEDS_pub_seed_bytes);

        LOG_MAT_FMT(M_hat_i, 2, MEDS_k, "M_hat[%i]", i);


        pmod_mat_t G0_prime[2 * MEDS_m * MEDS_n];

        pmod_mat_mul(G0_prime, 2, MEDS_m * MEDS_n, M_hat_i, 2, MEDS_k, G[0], MEDS_k, MEDS_m * MEDS_n);

        LOG_MAT_FMT(G0_prime, 2, MEDS_m * MEDS_n, "G0_prime[%i]", i);


        pmod_mat_t A_hat_i[MEDS_m * MEDS_m];
        pmod_mat_t B_hat_i[MEDS_n * MEDS_n];

        pmod_mat_t A_hat_inv[MEDS_m * MEDS_m];
        pmod_mat_t B_hat_inv[MEDS_n * MEDS_n];

        if (solve(A_hat_i, B_hat_inv, G0_prime, false) < 0)
        {
          LOG("no sol");
          continue;
        }

        if (pmod_mat_inv(B_hat_i, B_hat_inv, MEDS_n, MEDS_n) < 0)
        {
          LOG("no B_hat");
          continue;
        }

        if (pmod_mat_inv(A_hat_inv, A_hat_i, MEDS_m, MEDS_m) < 0)
        {
          LOG("no A_hat_inv");
          continue;
        }

        LOG_MAT_FMT(A_hat_i, MEDS_m, MEDS_m, "A_hat[%i]", i);
        LOG_MAT_FMT(B_hat_i, MEDS_n, MEDS_n, "B_hat[%i]", i);


        pi(G_hat_i, A_hat_i, B_hat_i, G[0]);

        LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);


        if (SF(G_hat_i, G_hat_i) == 0)
        {
          LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);
          break;
        }

        LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);
      }
    }

    {
      bitstream_t bs;
      uint8_t bs_buf[CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8)];

      bs_init(&bs, bs_buf, CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8));

      for (int r = 0; r < MEDS_k; r++)
        for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
          bs_write(&bs, G_hat_i[r * MEDS_m*MEDS_n + j], GFq_bits);

      shake256_absorb(&shake, bs_buf, CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8));
    }
  }

  shake256_absorb(&shake, (uint8_t*)(sm + MEDS_SIG_BYTES), smlen - MEDS_SIG_BYTES);

  shake256_finalize(&shake);

  uint8_t check[MEDS_digest_bytes];

  shake256_squeeze(check, MEDS_digest_bytes, &shake);

  if (memcmp(digest, check, MEDS_digest_bytes) != 0)
  {
    fprintf(stderr, "Signature verification failed!\n");

    return -1;
  }

  memcpy(m, (uint8_t*)(sm + MEDS_SIG_BYTES), smlen - MEDS_SIG_BYTES);
  *mlen = smlen - MEDS_SIG_BYTES;

  return 0;
}

