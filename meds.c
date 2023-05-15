#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"

#include "fips202.h"

#include "params.h"

#include "meds.h"
#include "seed.h"
#include "util.h"
#include "bitstream.h"

#include "matrixmod.h"

#define CEILING(x,y) (((x) + (y) - 1) / (y))


int keygen(uint8_t *sec_seed, int sec_seed_len, uint8_t *sk, int sk_len, uint8_t *pk, int pk_len)
{
  if (sec_seed_len != MEDS_sec_bytes)
    return -1;

  if (sk_len < ((MEDS_s-1)*MEDS_m*MEDS_m + (MEDS_s-1)*MEDS_n*MEDS_n + MEDS_k*(MEDS_m*MEDS_n-MEDS_k) + MEDS_sec_bytes))
    return -1;

  if (pk_len < ((MEDS_s-1)*MEDS_k*(MEDS_m*MEDS_n-MEDS_k) + MEDS_sec_bytes))
    return -1;


  // for (int i = 0; i < MEDS_sec_bytes; i++)
  //   printf("%i, ", sec_seed[i]);
  // printf("\n");


  pmod_mat_t *G0sub[MEDS_k];

  pmod_mat_t G_data[MEDS_k * MEDS_m * MEDS_n * MEDS_s];
  pmod_mat_t *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
    G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];

  LOG_VEC(sec_seed, MEDS_sec_bytes);

  keccak_state sec_shake;

  shake256_init(&sec_shake);
  shake256_absorb(&sec_shake, sec_seed, sec_seed_len);
  shake256_finalize(&sec_shake);

redo:

  uint8_t pub_seed[MEDS_sec_bytes];

  shake256_squeeze(pub_seed, MEDS_sec_bytes, &sec_shake);

  LOG_VEC(pub_seed, MEDS_sec_bytes);

  G_mat_init(G[0], G0sub);

  keccak_state pub_shake;

  shake256_init(&pub_shake);
  shake256_absorb(&pub_shake, pub_seed, sizeof(pub_seed));
  shake256_finalize(&pub_shake);

  rnd_sys_mat(G[0], MEDS_k, MEDS_m*MEDS_n, &pub_shake);

  LOG_MAT(G[0], MEDS_k, MEDS_m*MEDS_n);

  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_m * MEDS_m];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
  {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
  }

  for (int s = 1; s < MEDS_s; s++)
  {
    pmod_mat_t *Gsub[MEDS_k];

    G_mat_init(G[s], Gsub);

    pmod_mat_t Pj0[MEDS_m * MEDS_n] = {1, 0};
    pmod_mat_t Pj1[MEDS_m * MEDS_n] = {0, 1, 0};
    pmod_mat_t *Pj[2] = {Pj0, Pj1};

    rnd_matrix(&Pj[0][MEDS_k], 1, MEDS_m * MEDS_n - MEDS_k, &pub_shake);
    LOG_MAT(Pj[0], MEDS_m, MEDS_n);

    rnd_matrix(&Pj[1][MEDS_k], 1, (MEDS_m-1) * MEDS_n - MEDS_k, &pub_shake);
    rnd_matrix(&Pj[1][(MEDS_m-1)*MEDS_n], 1, MEDS_n, &sec_shake);
    LOG_MAT(Pj[1], MEDS_m, MEDS_n);

    while (1 == 1) // redo generation for this index until success
    //for (int repeat = 0; repeat < 100; repeat++)
    {
      // if (repeat > 5)
      // {
      //   printf("redo\n");
      //
      //   goto redo;
      // }
      //
      //printf("s: %i\n", s);

      pmod_mat_t Tj[MEDS_k * MEDS_k];

      rnd_inv_matrix(Tj, MEDS_k, MEDS_k, NULL, &sec_shake);

      LOG_MAT(Tj, MEDS_k, MEDS_k);

      pmod_mat_t G0prime_data[MEDS_k * MEDS_m * MEDS_n];

      pmod_mat_t *G0prime = G0prime_data;

      pmod_mat_mul(G0prime, MEDS_k, MEDS_m * MEDS_n,
          Tj, MEDS_k, MEDS_k,
          G[0], MEDS_k, MEDS_m * MEDS_n);

      LOG_MAT(G0prime, MEDS_k, MEDS_m * MEDS_n);

      GFq_t A00 = rnd_GF(&sec_shake);

      LOG_VAL(A00);

      pmod_mat_t rsys[(MEDS_m * MEDS_m + MEDS_n * MEDS_n - 1) * (MEDS_m * MEDS_m + MEDS_n * MEDS_n)] = {0};

      pmod_mat_t *tmp = rsys;

      // set up lineat eq system
      for (int l = 0; l < MEDS_n; l++)
        for (int j = 0; j < MEDS_m; j++)
        {
          for (int ii = (l == 0) ? 1 : 0; ii < MEDS_m; ii++)
            tmp[l*MEDS_m + ii - 1] = G0prime[ii*MEDS_m + j];

          for (int ii = 0; ii < MEDS_n; ii++)
            tmp[MEDS_m*MEDS_m + ii*MEDS_n + j - 1] = (MEDS_p - Pj[0][l*MEDS_m + ii]) % MEDS_p;

          if (l == 0)
            tmp[MEDS_m*MEDS_m + MEDS_n*MEDS_n - 1] = ((MEDS_p - A00) * G0prime[j]) % MEDS_p;

          tmp += MEDS_m * MEDS_m + MEDS_n * MEDS_n;
        }

      for (int l = 0; l < MEDS_n; l++)
        for (int j = 0; j < ((l < MEDS_n-1) ? MEDS_m : MEDS_m-1); j++)
        {
          for (int ii = (l == 0) ? 1 : 0; ii < MEDS_m; ii++)
            tmp[l*MEDS_m + ii - 1] = G0prime[MEDS_m * MEDS_n + ii*MEDS_m + j];

          for (int ii = 0; ii < MEDS_n; ii++)
            tmp[MEDS_m*MEDS_m + ii*MEDS_n + j - 1] = (MEDS_p - Pj[1][l*MEDS_m + ii]) % MEDS_p;

          if (l == 0)
            tmp[MEDS_m*MEDS_m + MEDS_n*MEDS_n - 1] = ((MEDS_p - A00) * G0prime[MEDS_m * MEDS_n + j]) % MEDS_p;

          tmp += MEDS_m * MEDS_m + MEDS_n * MEDS_n;
        }

      LOG_MAT(rsys, (MEDS_m * MEDS_m + MEDS_n * MEDS_n - 1), (MEDS_m * MEDS_m + MEDS_n * MEDS_n));

      // solve system
      if (pmod_mat_syst_ct(rsys, (MEDS_m * MEDS_m + MEDS_n * MEDS_n - 1), (MEDS_m * MEDS_m + MEDS_n * MEDS_n)) < 0)
      {
        LOG("no sol");
        continue;
      }

      GFq_t sol[MEDS_m * MEDS_m + MEDS_n * MEDS_n - 1];

      for (int i = 1; i < MEDS_m * MEDS_m + MEDS_n * MEDS_n; i++)
        sol[i-1] = rsys[i * (MEDS_m * MEDS_m + MEDS_n * MEDS_n)-1];


      LOG_VEC(sol, MEDS_m * MEDS_m + MEDS_n * MEDS_n - 1);


      pmod_mat_t A[MEDS_m * MEDS_m] = {0};
      pmod_mat_t B[MEDS_n * MEDS_n] = {0};

      A[0] = A00;

      for (int i = 1; i < MEDS_m*MEDS_m; i++)
        A[i] = sol[i-1];

      for (int i = 0; i < MEDS_n*MEDS_n; i++)
        B_inv[s][i] = sol[i + MEDS_m*MEDS_m - 1];

      int nA = 0;
      int nB = 0;

      if (pmod_mat_inv(B, B_inv[s], MEDS_n, MEDS_n) < 0)
      {
        LOG("no inv B");

        nB = 1;
      }

      if (pmod_mat_inv(A_inv[s], A, MEDS_m, MEDS_m) < 0)
      {
        LOG("no inv A_inv");

        nA = 1;
      }

      // In some cases we need entirely new Pj[0] and Pj[1] - so redo all.
      // This is over-cautious; in some cases a new Tj and A00 suffices (-> continue).
      if ((nA==1) && (nB==1))
      {
        LOG("A and B_inv not invertible; redo key gen.");
        goto redo;
      }

      if ((nA==1) || (nB==1))
      {
        continue;
      }

      LOG_MAT_FMT(A, MEDS_m, MEDS_m, "A[%i]", s);
      LOG_MAT_FMT(A_inv[s], MEDS_m, MEDS_m, "A_inv[%i]", s);
      LOG_MAT_FMT(B, MEDS_n, MEDS_n, "B[%i]", s);
      LOG_MAT_FMT(B_inv[s], MEDS_n, MEDS_n, "B_inv[%i]", s);

      for (int i = 0; i < MEDS_k; i++)
      {
        pmod_mat_mul(Gsub[i], MEDS_m, MEDS_n, A, MEDS_m, MEDS_m, G0sub[i], MEDS_m, MEDS_n);
        pmod_mat_mul(Gsub[i], MEDS_m, MEDS_n, Gsub[i], MEDS_m, MEDS_n, B, MEDS_n, MEDS_n);
      }

#ifdef BASE_CHANGE
      pmod_mat_t M[MEDS_k * MEDS_k];

      rnd_matrix(M, MEDS_k, MEDS_k, &sec_shake);

      pmod_mat_mul(G[s], MEDS_k, MEDS_m*MEDS_n, M, MEDS_k, MEDS_k, G[s], MEDS_k, MEDS_m*MEDS_n);

      LOG_MAT_FMT(G[s], MEDS_k, MEDS_m*MEDS_n, "G[%i]", s);

      if (pmod_mat_syst_non_ct(G[s], MEDS_k, MEDS_m*MEDS_n) != 0)
      {
        LOG_MAT_FMT(G[s], MEDS_k, MEDS_m*MEDS_n, "G[%i]", s);
        continue; // Not systematic; try again for index s.
      }
#else  // #ifdef BASE_CHANGE
      if (pmod_mat_syst_ct(G[s], MEDS_k, MEDS_m*MEDS_n) != 0)
      {
        LOG("redo G[%i]", s);
        continue; // Not systematic; try again for index s.
      }
#endif

      LOG_MAT_FMT(G[s], MEDS_k, MEDS_m*MEDS_n, "G[%i]", s);

      // successfull generated G[s]; break out of while loop
      break;
    }
  }

  // copy sk data
  uint8_t *tmp_sk = sk;

  memcpy(tmp_sk, sec_seed, sec_seed_len);
  LOG_VEC(tmp_sk, sec_seed_len, "sec_seed (sk)");
  tmp_sk += sec_seed_len;


  memcpy(tmp_sk, A_inv[1], (MEDS_s-1) * MEDS_m * MEDS_m * sizeof(GFq_t));
  LOG_VEC(tmp_sk, (MEDS_s-1) * MEDS_m * MEDS_m * sizeof(GFq_t), "A_inv (sk)");
  tmp_sk += (MEDS_s-1) * MEDS_m * MEDS_m * sizeof(GFq_t);

  memcpy(tmp_sk, B_inv[1], (MEDS_s-1) * MEDS_n * MEDS_n * sizeof(GFq_t));
  LOG_VEC(tmp_sk, (MEDS_s-1) * MEDS_n * MEDS_n * sizeof(GFq_t), "B_inv (sk)");
  tmp_sk += (MEDS_s-1) * MEDS_n * MEDS_n * sizeof(GFq_t);

  // copy right part only
  for (int r = 0; r < MEDS_k; r++)
  {
    memcpy(tmp_sk + r * (MEDS_m * MEDS_n - MEDS_k) * sizeof(GFq_t),
       &G_data[r * MEDS_m * MEDS_n + MEDS_k], (MEDS_m * MEDS_n - MEDS_k) * sizeof(GFq_t));
  }
  LOG_VEC(tmp_sk, MEDS_k * (MEDS_m * MEDS_n - MEDS_k) * sizeof(GFq_t), "G[0] (sk)");
  tmp_sk += MEDS_k * (MEDS_m * MEDS_n - MEDS_k) * sizeof(GFq_t);


  LOG_HEX(sk, SK_BYTES);


  // copy pk data
  uint8_t *tmp_pk = pk;

  memcpy(tmp_pk, pub_seed, MEDS_sec_bytes);
  LOG_VEC(tmp_pk, MEDS_sec_bytes, "pub_seed (pk)");
  tmp_pk += MEDS_sec_bytes;

  bitstream_t bs;

  bs_init(&bs, tmp_pk, pk_len - MEDS_sec_bytes);

  for (int si = 1; si < MEDS_s; si++)
  {
    for (int j = (MEDS_m-1)*MEDS_n; j < MEDS_m*MEDS_n; j++)
      bs_write(&bs, G[si][MEDS_m*MEDS_n + j], GFq_bits);

    for (int r = 2; r < MEDS_k; r++)
      for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
        bs_write(&bs, G[si][r*MEDS_m*MEDS_n + j], GFq_bits);
  }

  LOG_VEC(tmp_pk, pk_len - MEDS_sec_bytes, "G[1:] (pk)");
  tmp_pk += pk_len - MEDS_sec_bytes;

  LOG_HEX(pk, MEDS_PK_BYTES);

  if (MEDS_PK_BYTES != MEDS_sec_bytes + bs.byte_pos+1)
  {
    fprintf(stderr, "ERROR: MEDS_PK_BYTES and actual pk size do not match!\n");
    return -1;
  }

  return 0;
}

int sign(uint8_t *seed, int seed_len, uint8_t *sk, int sk_len, const char *msg, int msg_len, uint8_t *sig, int sig_len)
{
  if (sk_len < ((MEDS_s-1)*MEDS_m*MEDS_m + (MEDS_s-1)*MEDS_n*MEDS_n + MEDS_k*(MEDS_m*MEDS_n-MEDS_k) + MEDS_sec_bytes))
    return -1;

  if (sig_len < MEDS_SIG_BYTES)
    return -1;

  sk += MEDS_sec_bytes;

  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_n * MEDS_n];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
  {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
  }

 
  pmod_mat_t G_0[MEDS_k * MEDS_m * MEDS_n];

  for (int i = 1; i < MEDS_s; i++)
  {
    // load A_inv from secret key
    memcpy(A_inv[i], sk, (MEDS_m * MEDS_m) * sizeof(GFq_t));
    sk += (MEDS_m * MEDS_m) * sizeof(GFq_t);

    LOG_MAT(A_inv[i], MEDS_m, MEDS_m);
  }

  for (int i = 1; i < MEDS_s; i++)
  {
    // load B_inv from secret key
    memcpy(B_inv[i], sk, (MEDS_n * MEDS_n) * sizeof(GFq_t));
    sk += (MEDS_n * MEDS_n) * sizeof(GFq_t);

    LOG_MAT(B_inv[i], MEDS_n, MEDS_n);
  }


  pmod_mat_t *G0sub[MEDS_k];

  G_mat_init(G_0, G0sub);

  // load G_0 from secret key
  for (int r = 0; r < MEDS_k; r++)
  {
    memcpy(&G_0[r * MEDS_m * MEDS_n + MEDS_k], sk, (MEDS_m * MEDS_n - MEDS_k) * sizeof(GFq_t));
    sk += (MEDS_m * MEDS_n - MEDS_k) * sizeof(GFq_t);

    // set front part of G_0 to identity
    memset(&G_0[r * MEDS_m * MEDS_n], 0, MEDS_k * sizeof(GFq_t));
    pmod_mat_set_entry(G_0, MEDS_k, MEDS_m * MEDS_n, r, r, 1);
  }

  LOG_MAT(G_0, MEDS_k, MEDS_m*MEDS_m);

  LOG_VEC(seed, MEDS_sec_bytes);

  keccak_state shake_s;

  shake256_init(&shake_s);
  shake256_absorb(&shake_s, seed, seed_len);
  shake256_finalize(&shake_s);

  uint8_t *seeds;

#ifdef HAVE_SEED_TREE
  uint8_t stree[MEDS_sec_bytes * ((1 << (CLOG2(MEDS_t-1)+2)) - 1)] = {0};

  shake256_squeeze(stree, MEDS_sec_bytes, &shake_s);

  seeds = &stree[MEDS_sec_bytes * TREE_ADR(CLOG2(MEDS_t-1), 0)];

  uint8_t salt[MEDS_salt_sec_bytes];

  shake256_squeeze(salt, MEDS_salt_sec_bytes, &shake_s);

  t_hash(stree, salt, 0, 0);
#else
  uint8_t seeds_data[MEDS_t*MEDS_sec_bytes];

  seeds = seeds_data;

  shake256_squeeze(seeds, sizeof(seeds_data), &shake_s);
#endif


  keccak_state h_shake;
  shake256_init(&h_shake);


  pmod_mat_t A_tilde_data[MEDS_t * MEDS_m * MEDS_m];
  pmod_mat_t B_tilde_data[MEDS_t * MEDS_m * MEDS_m];

  pmod_mat_t *A_tilde[MEDS_t];
  pmod_mat_t *B_tilde[MEDS_t];

  for (int i = 0; i < MEDS_t; i++)
  {
    A_tilde[i] = &A_tilde_data[i * MEDS_m * MEDS_m];
    B_tilde[i] = &B_tilde_data[i * MEDS_n * MEDS_n];
  }

  for (int i = 0; i < MEDS_t; i++)
  {
    pmod_mat_t G_tilde_ti[MEDS_k * MEDS_m * MEDS_m];

    while (1 == 1)
    {
      keccak_state shake;

      shake256_init(&shake);
      shake256_absorb(&shake, &seeds[i*MEDS_sec_bytes], MEDS_sec_bytes);
      shake256_finalize(&shake);

      rnd_inv_matrix(A_tilde[i], MEDS_m, MEDS_m, NULL, &shake);
      rnd_inv_matrix(B_tilde[i], MEDS_n, MEDS_n, NULL, &shake);

      LOG_MAT_FMT(A_tilde[i], MEDS_m, MEDS_m, "A_tilde[%i]", i);
      LOG_MAT_FMT(B_tilde[i], MEDS_n, MEDS_n, "B_tilde[%i]", i);


      pmod_mat_t *Gsub[MEDS_k];

      G_mat_init(G_tilde_ti, Gsub);

      for (int j = 0; j < MEDS_k; j++)
      {
        pmod_mat_mul(Gsub[j], MEDS_m, MEDS_n, A_tilde[i], MEDS_m, MEDS_m, G0sub[j], MEDS_m, MEDS_n);
        pmod_mat_mul(Gsub[j], MEDS_m, MEDS_n, Gsub[j], MEDS_m, MEDS_n, B_tilde[i], MEDS_n, MEDS_n);
      }

      LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, "G_tilde[%i]", i);

#ifdef BASE_CHANGE
      pmod_mat_t M_tilde[MEDS_k*MEDS_k];

#ifdef HAVE_SEED_TREE
      rnd_matrix(M_tilde, MEDS_k, MEDS_k, &shake);
#else
      rnd_matrix(M_tilde, MEDS_k, MEDS_k, &shake_s);
#endif // HAVE_SEED_TREE

      pmod_mat_mul(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, M_tilde, MEDS_k, MEDS_k, G_tilde_ti, MEDS_k, MEDS_m*MEDS_n);

      LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, "G_tilde[%i]", i);

      if (pmod_mat_syst_non_ct(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n) == 0)
        break;

#else // BASE_CHANGE
      if (pmod_mat_syst_ct(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n) == 0)
        break;
#endif // BASE_CHANGE

      shake256_squeeze(&seeds[i*MEDS_sec_bytes], MEDS_sec_bytes, &shake);
    }

    LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, "G_tilde[%i]", i);

    for (int r = 0; r < MEDS_k; r++)
      shake256_absorb(&h_shake, (uint8_t*)&G_tilde_ti[r * MEDS_m*MEDS_n + MEDS_k],
          (MEDS_m*MEDS_n - MEDS_k) * sizeof(GFq_t));
  }

  shake256_absorb(&h_shake, (uint8_t*)msg, msg_len);

  shake256_finalize(&h_shake);

  uint8_t digest[MEDS_sec_bytes];

  shake256_squeeze(digest, MEDS_sec_bytes, &h_shake);

  LOG_VEC(digest, MEDS_sec_bytes);


  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_sec_bytes, h, MEDS_t);

  LOG_VEC(h, MEDS_t);

#ifdef HAVE_SEED_TREE
  uint8_t *path = sig;

  unsigned int indices[MEDS_w] = {0};

  int idx = 0;

  for (int i = 0; i < MEDS_t; i++)
    if (h[i] > 0)
      indices[idx++] = i;

  unsigned int empty_slots[MEDS_w] = {0};

  stree_to_path(stree, indices, path, empty_slots, CEILING((MEDS_n*MEDS_n + MEDS_m*MEDS_m) * GFq_bits, 8), salt);

  int mid = 0;
#endif 

#ifdef HAVE_SEED_TREE
#else
  uint8_t *sig_h = sig;
#endif

  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] > 0)
    {
#ifdef HAVE_SEED_TREE
      // set sig pointer to next empty slot
      uint8_t *sig_h = &path[empty_slots[mid++]];
#endif
      bitstream_t bs;

      bs_init(&bs, sig_h, CEILING((MEDS_n*MEDS_n + MEDS_m*MEDS_m) * GFq_bits, 8)+10);

      {
        pmod_mat_t mu[MEDS_m*MEDS_m];

        pmod_mat_mul(mu, MEDS_m, MEDS_m, A_tilde[i], MEDS_m, MEDS_m, A_inv[h[i]], MEDS_m, MEDS_m);

        LOG_MAT(mu, MEDS_m, MEDS_m);

        for (int j = 0; j < MEDS_m*MEDS_m; j++)
          bs_write(&bs, mu[j], GFq_bits);
      }

      {
        pmod_mat_t nu[MEDS_n*MEDS_n];

        pmod_mat_mul(nu, MEDS_n, MEDS_n, B_inv[h[i]], MEDS_n, MEDS_n, B_tilde[i], MEDS_n, MEDS_n);

        LOG_MAT(nu, MEDS_n, MEDS_n);

        for (int j = 0; j < MEDS_n*MEDS_n; j++)
          bs_write(&bs, nu[j], GFq_bits);
      }

      sig_h += bs.byte_pos + (bs.bit_pos > 1 ? 1 : 0);
    }
    else
    {
#ifdef HAVE_SEED_TREE
#else
      memcpy(sig_h, &seeds[i*MEDS_sec_bytes], MEDS_sec_bytes);
      sig_h += MEDS_sec_bytes;
#endif
    }
  }

#ifdef HAVE_SEED_TREE
  memcpy(sig + MEDS_SIG_BYTES - MEDS_sec_bytes - MEDS_salt_sec_bytes, digest, MEDS_sec_bytes);
  memcpy(sig + MEDS_SIG_BYTES - MEDS_salt_sec_bytes, salt, MEDS_salt_sec_bytes);
#else
  memcpy(sig + MEDS_SIG_BYTES - MEDS_sec_bytes, digest, MEDS_sec_bytes);
#endif

  LOG_HEX(sig, MEDS_SIG_BYTES);

  return 0;
}

char* verify(uint8_t *pk, int pk_len, char *msg, int msg_len, uint8_t *sig, int sig_len)
{
  pmod_mat_t G_data[MEDS_k*MEDS_m*MEDS_n * MEDS_s];
  pmod_mat_t *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
    G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];

  pmod_mat_t *Gsub[MEDS_s][MEDS_k];

  G_mat_init(G[0], Gsub[0]);


  keccak_state pub_shake;

  shake256_init(&pub_shake);
  shake256_absorb(&pub_shake, pk, MEDS_sec_bytes);
  shake256_finalize(&pub_shake);

  bitstream_t bs;

  bs_init(&bs, pk + MEDS_sec_bytes, pk_len - MEDS_sec_bytes);

  rnd_sys_mat(G[0], MEDS_k, MEDS_m*MEDS_n, &pub_shake);

  for (int i = 1; i < MEDS_s; i++)
  {
    G_mat_init(G[i], Gsub[i]);

    for (int r = 0; r < MEDS_k; r++)
      for (int c = 0; c < MEDS_k; c++)
        if (r == c)
          pmod_mat_set_entry(G[i], MEDS_k, MEDS_m * MEDS_n, r, c, 1);
        else
          pmod_mat_set_entry(G[i], MEDS_k, MEDS_m * MEDS_n, r, c, 0);

    for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
      G[i][j] = rnd_GF(&pub_shake);

    for (int j = MEDS_k; j < (MEDS_m-1)*MEDS_n; j++)
      G[i][MEDS_m*MEDS_n + j] = rnd_GF(&pub_shake);


    for (int j = (MEDS_m-1)*MEDS_n; j < MEDS_m*MEDS_n; j++)
      G[i][MEDS_m*MEDS_n + j] = bs_read(&bs, GFq_bits);

    for (int r = 2; r < MEDS_k; r++)
      for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
        G[i][r*MEDS_m*MEDS_n + j] = bs_read(&bs, GFq_bits);
  }

  for (int i = 0; i < MEDS_s; i++)
    LOG_MAT_FMT(G[i], MEDS_k, MEDS_m*MEDS_n, "G[%i]", i);

#ifdef HAVE_SEED_TREE
  uint8_t *digest = (uint8_t*)sig + (MEDS_SIG_BYTES - MEDS_sec_bytes - MEDS_salt_sec_bytes);

  uint8_t *salt = (uint8_t*)sig + (MEDS_SIG_BYTES - MEDS_salt_sec_bytes);
#else
  uint8_t *digest = (uint8_t*)sig + (MEDS_SIG_BYTES - MEDS_sec_bytes);
#endif

  uint8_t *path = sig;

  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_sec_bytes, h, MEDS_t);

#ifdef HAVE_SEED_TREE
  unsigned int indices[MEDS_w] = {0};
  unsigned int empty_slots[MEDS_w] = {0};

  uint8_t stree[MEDS_sec_bytes * ((1 << (CLOG2(MEDS_t-1)+2)) - 1)] = {0};

  {
    int idx = 0;

    for (int i = 0; i < MEDS_t; i++)
      if (h[i] > 0)
        indices[idx++] = i;
  }

  path_to_stree(stree, indices, path, empty_slots, CEILING((MEDS_n*MEDS_n + MEDS_m*MEDS_m) * GFq_bits, 8), salt);

#else
  uint8_t *seeds = path;
#endif


  pmod_mat_t G_hat_data[MEDS_k*MEDS_m*MEDS_n * MEDS_t];
  pmod_mat_t *G_hat[MEDS_t];

  for (int i = 0; i < MEDS_t; i++)
    G_hat[i] = &G_hat_data[i * MEDS_k * MEDS_m * MEDS_n];


  pmod_mat_t mu[MEDS_m*MEDS_m];
  pmod_mat_t nu[MEDS_n*MEDS_n];

#ifdef HAVE_SEED_TREE
  int mid = 0;
#endif

  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] == 0)
    {
      keccak_state shake;

#ifdef HAVE_SEED_TREE
      uint8_t *seeds = &stree[MEDS_sec_bytes * (TREE_ADR(CLOG2(MEDS_t-1), 0) + i)];

      while (1 == 1)
#endif
      {
        LOG_VEC_FMT(seeds, MEDS_sec_bytes, "seeds[%i]", i);

        shake256_init(&shake);
        shake256_absorb(&shake, seeds, MEDS_sec_bytes);
        shake256_finalize(&shake);

        seeds += MEDS_sec_bytes;

        pmod_mat_t A_tilde[MEDS_m*MEDS_m];
        pmod_mat_t B_tilde[MEDS_n*MEDS_n];

        rnd_inv_matrix(A_tilde, MEDS_m, MEDS_m, NULL, &shake);
        rnd_inv_matrix(B_tilde, MEDS_n, MEDS_n, NULL, &shake);

        LOG_MAT_FMT(A_tilde, MEDS_m, MEDS_m, "A_tilde[%i]", i);
        LOG_MAT_FMT(B_tilde, MEDS_n, MEDS_n, "B_tilde[%i]", i);

        pmod_mat_t *G_hatsub[MEDS_k];

        G_mat_init(G_hat[i], G_hatsub);

        for (int j = 0; j < MEDS_k; j++)
        {
          pmod_mat_mul(G_hatsub[j], MEDS_m, MEDS_n, A_tilde, MEDS_m, MEDS_m, Gsub[0][j], MEDS_m, MEDS_n);
          pmod_mat_mul(G_hatsub[j], MEDS_m, MEDS_n, G_hatsub[j], MEDS_m, MEDS_n, B_tilde, MEDS_n, MEDS_n);
        }

        LOG_MAT_FMT(G_hat[i], MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);

#if defined(BASE_CHANGE) && defined(HAVE_SEED_TREE)
        pmod_mat_t M_tilde[MEDS_k*MEDS_k];

        rnd_matrix(M_tilde, MEDS_k, MEDS_k, &shake);

        pmod_mat_mul(G_hat[i], MEDS_k, MEDS_m*MEDS_n, M_tilde, MEDS_k, MEDS_k, G_hat[i], MEDS_k, MEDS_m*MEDS_n);

        LOG_MAT_FMT(G_hat[i], MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);
#endif // defined(BASE_CHANGE) && defined(HAVE_SEED_TREE)

        int syst_success = pmod_mat_syst_non_ct(G_hat[i], MEDS_k, MEDS_m*MEDS_n);

        LOG_MAT_FMT(G_hat[i], MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);

#ifdef HAVE_SEED_TREE
        if (syst_success == 0)
          break;

        // reset to current seed slot nd retry with new seed
        seeds -= MEDS_sec_bytes;

        shake256_squeeze(seeds, MEDS_sec_bytes, &shake);
#else
        (void)(syst_success);
#endif

        //printf("redo G_hat[%i]\n", i);
      }
    }
    else
    {
#ifdef HAVE_SEED_TREE
      uint8_t *seeds = &path[empty_slots[mid++]];
#endif

      LOG_VEC_FMT(seeds, CEILING((MEDS_n*MEDS_n + MEDS_m*MEDS_m) * GFq_bits, 8), "seeds[%i]", i);

      bitstream_t bs;

      bs_init(&bs, seeds, CEILING((MEDS_n*MEDS_n + MEDS_m*MEDS_m) * GFq_bits, 8));

#ifdef HAVE_SEED_TREE
#else
      seeds += CEILING((MEDS_n*MEDS_n + MEDS_m*MEDS_m) * GFq_bits, 8);
#endif

      for (int j = 0; j < MEDS_m*MEDS_m; j++)
        mu[j] = bs_read(&bs, GFq_bits);

      for (int j = 0; j < MEDS_n*MEDS_n; j++)
        nu[j] = bs_read(&bs, GFq_bits);


      LOG_MAT_FMT(mu, MEDS_m, MEDS_m, "mu[%i]", i);
      LOG_MAT_FMT(nu, MEDS_n, MEDS_n, "nu[%i]", i);

      pmod_mat_t *G_hatsub[MEDS_k];

      G_mat_init(G_hat[i], G_hatsub);

      for (int j = 0; j < MEDS_k; j++)
      {
        pmod_mat_mul(G_hatsub[j], MEDS_m, MEDS_n, mu, MEDS_m, MEDS_m, Gsub[h[i]][j], MEDS_m, MEDS_n);
        pmod_mat_mul(G_hatsub[j], MEDS_m, MEDS_n, G_hatsub[j], MEDS_m, MEDS_n, nu, MEDS_n, MEDS_n);
      }

      LOG_MAT_FMT(G_hat[i], MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);

      pmod_mat_syst_non_ct(G_hat[i], MEDS_k, MEDS_m*MEDS_n);

      LOG_MAT_FMT(G_hat[i], MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);
    }
  }

  keccak_state shake;
  shake256_init(&shake);

  for (int i = 0; i < MEDS_t; i++)
    for (int r = 0; r < MEDS_k; r++)
      shake256_absorb(&shake, (uint8_t*)&G_hat[i][r * MEDS_m*MEDS_n + MEDS_k],
          (MEDS_m*MEDS_n - MEDS_k) * sizeof(GFq_t));

  shake256_absorb(&shake, (uint8_t*)msg, msg_len);

  shake256_finalize(&shake);

  uint8_t check[MEDS_sec_bytes];

  shake256_squeeze(check, MEDS_sec_bytes, &shake);

  if (memcmp(digest, check, MEDS_sec_bytes) != 0)
  {
    fprintf(stderr, "Signature verification failed!\n");

    return NULL;
    //exit(-1);
  }

  return msg;
}

