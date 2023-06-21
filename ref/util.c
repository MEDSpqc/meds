#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fips202.h"

#include "log.h"

#include "util.h"
#include "meds.h"
#include "matrixmod.h"


void XOF(uint8_t **buf, size_t *length, const uint8_t *seed, size_t seed_len, int num)
{
  keccak_state shake;

  shake256_absorb_once(&shake, seed, seed_len);

  for (int i = 0; i < num; i++)
    shake256_squeeze(buf[i], length[i], &shake);
}

GFq_t rnd_GF(keccak_state *shake)
{
  GFq_t val = MEDS_p;

  while (val >= MEDS_p)
  {
    uint8_t data[sizeof(GFq_t)];
    
    shake256_squeeze(data, sizeof(GFq_t), shake);

    val = 0;
    for (int i = 0; i < sizeof(GFq_t); i++)
      val |= data[i] << (i*8);

    val = val & ((1 << GFq_bits) - 1);
  }

  return val;
}

void rnd_sys_mat(pmod_mat_t *M, int M_r, int M_c, const uint8_t *seed, size_t seed_len)
{
  keccak_state shake;
  shake256_absorb_once(&shake, seed, seed_len);

  for (int r = 0; r < M_r; r++)
    for (int c = M_r; c < M_c; c++)
      pmod_mat_set_entry(M, M_r, M_c, r, c, rnd_GF(&shake));

  for (int r = 0; r < M_r; r++)
    for (int c = 0; c < M_r; c++)
      if (r == c)
        pmod_mat_set_entry(M, M_r, M_c, r, c, 1);
      else
        pmod_mat_set_entry(M, M_r, M_c, r, c, 0);
}

void rnd_inv_matrix(pmod_mat_t *M, int M_r, int M_c, uint8_t *seed, size_t seed_len)
{
  keccak_state shake;
  shake256_absorb_once(&shake, seed, seed_len);

  while (0==0)
  {
redo:
    for (int r = 0; r < M_r; r++)
      for (int c = 0; c < M_c; c++)
        pmod_mat_set_entry(M, M_r, M_c, r, c, rnd_GF(&shake));

    pmod_mat_t tmp[M_r * M_c];

    memcpy(tmp, M, M_r * M_c * sizeof(GFq_t));

    {
      for (int r = 0; r < M_r; r++)
      {
        // swap
        for (int r2 = r+1; r2 < M_r; r2++)
        {
          uint64_t Mrr = pmod_mat_entry(tmp, M_r, M_c, r, r);

          for (int c = r; c < M_c; c++)
          {
            uint64_t val = pmod_mat_entry(tmp, M_r, M_c, r2, c);

            uint64_t Mrc = pmod_mat_entry(tmp, M_r, M_c, r, c);

            pmod_mat_set_entry(tmp, M_r, M_c, r, c, (Mrc + val * (Mrr == 0)) % MEDS_p);
          }
        }

        uint64_t val = pmod_mat_entry(tmp, M_r, M_c, r, r);

        if (val == 0)
          goto redo;

        val = GF_inv(val);

        // normalize
        for (int c = r; c < M_c; c++)
        {
          uint64_t tmp0 = ((uint64_t)pmod_mat_entry(tmp, M_r, M_c, r, c) * val) % MEDS_p;
          pmod_mat_set_entry(tmp, M_r, M_c, r, c, tmp0);
        }

        // eliminate
        for (int r2 = r+1; r2 < M_r; r2++)
        {
          uint64_t factor = pmod_mat_entry(tmp, M_r, M_c, r2, r);

          for (int c = r; c < M_c; c++)
          {
            uint64_t tmp0 = pmod_mat_entry(tmp, M_r, M_c, r, c);
            uint64_t tmp1 = pmod_mat_entry(tmp, M_r, M_c, r2, c);

            int64_t val = (tmp0 * factor) % MEDS_p;

            val = tmp1 - val;

            val += MEDS_p * (val < 0);

            pmod_mat_set_entry(tmp, M_r, M_c,  r2, c, val);
          }
        }
      }

      return;
    }
  }
}

int parse_hash(uint8_t *digest, int digest_len, uint8_t *h, int len_h)
{
  if (len_h < MEDS_t)
    return -1;

  //LOG_VEC(digest, MEDS_digest_bytes);
  //LOG_VAL(digest_len);

  keccak_state shake;

  shake256_init(&shake);
  shake256_absorb(&shake, digest, digest_len);
  shake256_finalize(&shake);

  for (int i = 0; i < MEDS_t; i++)
    h[i] = 0;

  int i = 0;
  while (i < MEDS_w)
  {
    uint64_t pos = 0;
    uint8_t buf[MEDS_t_bytes] = {0};

    shake256_squeeze(buf, MEDS_t_bytes, &shake);

    for (int j = 0; j < MEDS_t_bytes; j++)
      pos |= buf[j] << (j*8);

    pos = pos & MEDS_t_mask;

    if (pos >= MEDS_t)
      continue;

    if (h[pos] > 0)
      continue;

    //LOG("pos: %lu", pos);

    uint8_t val = 0;

    while ((val < 1) || (val > MEDS_s-1))
    {
      shake256_squeeze(&val, 1, &shake);
      val = val & MEDS_s_mask;
    }

    //LOG("p: %lu  v: %u", pos, val);

    h[pos] = val;

    i++;
  }

  return 0;
}

int solve_symb(pmod_mat_t *A, pmod_mat_t *B_inv, pmod_mat_t *G0prime, bool partial)
{
  _Static_assert (MEDS_n == MEDS_m+1, "solve_symb requires MEDS_n == MEDS_m+1");

  pmod_mat_t Pj0[MEDS_m * MEDS_n] = {0};
  pmod_mat_t Pj1[MEDS_m * MEDS_n] = {0};
  pmod_mat_t *Pj[2] = {Pj0, Pj1};

  for (int i = 0; i < MEDS_m; i++)
  {
    Pj0[i * (MEDS_m + 1)] = 1;
    Pj1[i * (MEDS_m + 1) + 1] = 1;
  }

  pmod_mat_t rsys[(2 * MEDS_m * MEDS_n) * (MEDS_m * MEDS_m + MEDS_n * MEDS_n)] = {0};

  pmod_mat_t *tmp = rsys;

  // set up lineat eq system
  for (int l = 0; l < MEDS_m; l++)
    for (int j = 0; j < MEDS_n; j++)
    {
      for (int ii = 0; ii < MEDS_n; ii++)
        tmp[ii*MEDS_n + j] = Pj[0][l*MEDS_m + ii];

      for (int ii = 0; ii < MEDS_m; ii++)
        tmp[MEDS_n*MEDS_n + l*MEDS_m + ii] = (MEDS_p - G0prime[ii*MEDS_n + j]) % MEDS_p;

      tmp += MEDS_m * MEDS_m + MEDS_n * MEDS_n;
    }

  for (int l = 0; l < MEDS_m; l++)
    for (int j = 0; j < MEDS_n; j++)
    {
      for (int ii = 0; ii < MEDS_n; ii++)
        tmp[ii*MEDS_n + j] = Pj[1][l*MEDS_m + ii];

      for (int ii = 0; ii < MEDS_m; ii++)
        tmp[MEDS_n*MEDS_n + l*MEDS_m + ii] = (MEDS_p - G0prime[MEDS_m * MEDS_n + ii*MEDS_n + j]) % MEDS_p;

      tmp += MEDS_m * MEDS_m + MEDS_n * MEDS_n;
    }

  //LOG_MAT(rsys, (2 * MEDS_m * MEDS_n), (MEDS_m * MEDS_m + MEDS_n * MEDS_n));

  int r = 0;

  // solve system
  if ((r = pmod_mat_syst_ct_partial(rsys, (2 * MEDS_m * MEDS_n), (MEDS_m * MEDS_m + MEDS_n * MEDS_n), partial)) != 0)
  {
    if (partial)
    {
      //LOG_MAT(rsys, (2 * MEDS_m * MEDS_n -1), (MEDS_m * MEDS_m + MEDS_n * MEDS_n));

      GFq_t sol[MEDS_m * MEDS_m + MEDS_n * MEDS_n];

      for (int i = 0; i < r; i++)
        sol[i] = rsys[i * (MEDS_m * MEDS_m + MEDS_n * MEDS_n) + r];

      sol[r] = MEDS_p - 1;

      for (int i = r+1; i < MEDS_m * MEDS_m + MEDS_n * MEDS_n; i++)
        sol[i] = 0;

      LOG_VEC(sol, MEDS_m * MEDS_m + MEDS_n * MEDS_n);


      for (int i = 0; i < MEDS_m*MEDS_m; i++)
        A[i] = sol[i + MEDS_n*MEDS_n];

      for (int i = 0; i < MEDS_n*MEDS_n; i++)
        B_inv[i] = sol[i];

      LOG_MAT(A, MEDS_m, MEDS_m);
      LOG_MAT(B_inv, MEDS_n, MEDS_n);

      return 0;
    }
    else
    {
      LOG("no sol");
      return -1;
    }
  }

  //LOG_MAT(rsys, (2 * MEDS_m * MEDS_n), (MEDS_m * MEDS_m + MEDS_n * MEDS_n));

  GFq_t sol[MEDS_m * MEDS_m + MEDS_n * MEDS_n];

  for (int i = 0; i < 2 * MEDS_m * MEDS_n; i++)
    sol[i] = rsys[i * (MEDS_m * MEDS_m + MEDS_n * MEDS_n) + MEDS_m * MEDS_m + MEDS_n * MEDS_n - 1];

  sol[MEDS_m * MEDS_m + MEDS_n * MEDS_n - 1] = MEDS_p - 1;


  LOG_VEC(sol, MEDS_m * MEDS_m + MEDS_n * MEDS_n);


  for (int i = 0; i < MEDS_m*MEDS_m; i++)
    A[i] = sol[i + MEDS_n*MEDS_n];

  for (int i = 0; i < MEDS_n*MEDS_n; i++)
    B_inv[i] = sol[i];

  LOG_MAT(A, MEDS_m, MEDS_m);
  LOG_MAT(B_inv, MEDS_n, MEDS_n);

  return 0;
}

int solve_opt(pmod_mat_t *A, pmod_mat_t *B_inv, pmod_mat_t *G0prime, bool partial)
{
  _Static_assert (MEDS_n == MEDS_m+1, "solve_opt requires MEDS_n == MEDS_m+1");

  // set up core sub-system
  pmod_mat_t N[MEDS_n * (2 * MEDS_m)] = {0};

  for (int i = 0; i < MEDS_m; i++)
    for (int j = 0; j < MEDS_n; j++)
    {
      pmod_mat_set_entry(N, MEDS_n, 2 * MEDS_m, j, i,
          (MEDS_p - pmod_mat_entry((G0prime + (MEDS_m*MEDS_n)), MEDS_m, MEDS_n, i, j)) % MEDS_p);

      pmod_mat_set_entry(N, MEDS_n, 2 * MEDS_m, j, MEDS_m+i,
          pmod_mat_entry(G0prime, MEDS_m, MEDS_n, i, j));
    }

  LOG_MAT(N, MEDS_n, 2 * MEDS_m);


  // Systemize core sub-system.
  // We can tolarate partial solving on verifier side as long as we can pivot all except the last row.
  int piv;
  if ((piv = pmod_mat_syst_ct_partial(N, MEDS_n, 2 * MEDS_m, partial)) != 0)
  {
    if (partial)
    {
      if (piv != MEDS_n-1)
      {
        return -1;  // This should not happen for a valid signaure.
      }
    }
    else
      return -1;
  }

  LOG_MAT(N, MEDS_n, 2 * MEDS_m);


  // Extract 2nd sub-system.
  pmod_mat_t N1[(MEDS_m-1) * MEDS_m] = {0};

  for (int i = 0; i < MEDS_m; i++)
    {
      pmod_mat_set_entry(N1, MEDS_n, MEDS_m, 0, i,
          pmod_mat_entry(N, MEDS_n, 2 * MEDS_m, MEDS_n-1, MEDS_m + i));
    }

  LOG_MAT(N1, MEDS_m-1, 2 * MEDS_m);


  // Remove front diagonale of core sub-system.
  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_m; j++)
        pmod_mat_set_entry(N, MEDS_n-1, MEDS_m, i, j,
            pmod_mat_entry(N, MEDS_n-1, 2 * MEDS_m, i, MEDS_m + j));

  // Consider N MEDS_n-1 x MEDS_m, now...
  LOG_MAT(N, MEDS_n-1, MEDS_m);


  // Fully reduce 2nd sub-system.
  for (int row = 1; row < MEDS_m-1; row++)
    for (int i = 0; i < MEDS_m; i++)
      for (int j = 0; j < MEDS_m; j++)
      {
        uint64_t tmp0 = pmod_mat_entry(N1, MEDS_m-1, MEDS_m, row-1, i);
        uint64_t tmp1 = pmod_mat_entry(N, MEDS_n-1, MEDS_m, i, j);

        uint64_t prod = (tmp0 * tmp1) % MEDS_p;

        uint64_t tmp2 = pmod_mat_entry(N1, MEDS_m-1, MEDS_m, row, j);

        int64_t diff = (MEDS_p + tmp2 - prod) % MEDS_p;

        pmod_mat_set_entry(N1, MEDS_m-1, MEDS_m, row, j, diff);
      }


  LOG_MAT(N1, MEDS_m-1, MEDS_m);


  GFq_t sol[MEDS_m * MEDS_m + MEDS_n * MEDS_n] = {0};


  int N1_r;

  // Sytemize 2nd sub-system.
  if ((N1_r = pmod_mat_syst_ct_partial(N1, MEDS_m-1, MEDS_m, partial)) != 0)
  {
    // The verifier can continune if the 2nd sub-system only partially systemizes.
    if (partial)
    {
      LOG_MAT(N1, MEDS_m-1, MEDS_m);

      // Fill in solutions from the 2nd sub-system.
      for (int i = 0; i < N1_r; i++)
        sol[2*MEDS_m*MEDS_n - (MEDS_m-1) + i] = pmod_mat_entry(N1, MEDS_m-1, MEDS_m, i, N1_r);

      // Fix solution at pivot location to -1.
      sol[2*MEDS_m*MEDS_n - (MEDS_m-1) + N1_r] = MEDS_p - 1;

      // Remaining spots resolve to zero.
      for (int i = N1_r+1; i < MEDS_m-1; i++)
        sol[2*MEDS_m*MEDS_n - (MEDS_m-1) + i] = 0;
    }
    else
      return -1;
  }
  else
  {
    LOG_MAT(N1, MEDS_m-1, MEDS_m);

    // Fill in solutions from the 2nd sub-system.
    for (int i = 0; i < MEDS_m-1; i++)
      sol[2*MEDS_m*MEDS_n - (MEDS_m-1) + i] = pmod_mat_entry(N1, MEDS_m-1, MEDS_m, i, MEDS_m-1);

    // Fix to -1.
    sol[MEDS_m*MEDS_m + MEDS_n*MEDS_n - 1] = MEDS_p - 1;

    N1_r = MEDS_m-1;
  }

  LOG_MAT(N, MEDS_n-1, MEDS_m);

  // Fill in solutions from core sub-system.
  for (int i = 0; i < MEDS_m; i++)
    sol[2*MEDS_m*MEDS_n - MEDS_m - (MEDS_m-1) + i] = pmod_mat_entry(N, MEDS_n-1, MEDS_m, i, N1_r);

  LOG_VEC(sol, MEDS_m*MEDS_m + MEDS_n*MEDS_n);


  // Back-substitue second last block.
  for (int c = N1_r - 1; c >= 0; c--)
    for (int r = 0; r < MEDS_m; r++)
    {
      uint64_t tmp1 = pmod_mat_entry(N, MEDS_n-1, MEDS_m, r, c);
      uint64_t tmp2 = sol[2*MEDS_m*MEDS_n - (MEDS_m-1) + c];

      uint64_t prod = (tmp1 * tmp2) % MEDS_p;

      uint64_t val = sol[2*MEDS_m*MEDS_n - MEDS_m - (MEDS_m-1) + r];

      val = ((MEDS_p + val) - prod) % MEDS_p;

      sol[2*MEDS_m*MEDS_n - MEDS_m - (MEDS_m-1) + r] = val;
    }

  LOG_VEC(sol, MEDS_m*MEDS_m + MEDS_n*MEDS_n);


  // Prepare substitution-source of first block in bottom part.
  pmod_mat_t P01nt[MEDS_n * MEDS_m] = {0};

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_m; j++)
      pmod_mat_set_entry(P01nt, MEDS_n, MEDS_m, i, j,
          (MEDS_p - pmod_mat_entry((G0prime + (MEDS_m*MEDS_n)), MEDS_m, MEDS_n, j, i)) % MEDS_p);

  LOG_MAT(P01nt, MEDS_n, MEDS_m);

  // Transfer last column to solution vector.
  for (int i = 0; i < MEDS_n; i++)
    sol[MEDS_m*MEDS_n + i] = pmod_mat_entry(P01nt, MEDS_n, MEDS_m, i, N1_r);

  LOG_VEC(sol, MEDS_m*MEDS_m + MEDS_n*MEDS_n);

  // Perfrom back-substution for first block in bottom part.
  for (int c = N1_r - 1; c >= 0; c--)
    for (int r = 0; r < MEDS_n; r++)
    {
      uint64_t tmp1 = pmod_mat_entry(P01nt, MEDS_n, MEDS_m, r, c);
      uint64_t tmp2 = pmod_mat_entry(N1, MEDS_m-1, MEDS_m, c, N1_r);

      uint64_t prod = (tmp1 * tmp2) % MEDS_p;

      uint64_t val = sol[MEDS_m*MEDS_n + r];

      val = ((MEDS_p + val) - prod) % MEDS_p;

      sol[MEDS_m*MEDS_n + r] = val;
    }

  LOG_VEC(sol, MEDS_m*MEDS_m + MEDS_n*MEDS_n);


  // Prepare substitution-source for top part.
  pmod_mat_t P00nt[MEDS_n * MEDS_m] = {0};

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_m; j++)
      pmod_mat_set_entry(P00nt, MEDS_n, MEDS_m, i, j,
          (MEDS_p - pmod_mat_entry(G0prime, MEDS_m, MEDS_n, j, i)) % MEDS_p);

  LOG_MAT(P00nt, MEDS_n, MEDS_m);


  // For the last block in the top part, transfer last column to solution vector.
  for (int i = 0; i < MEDS_n; i++)
    sol[(MEDS_m-1)*MEDS_n + i] = pmod_mat_entry(P00nt, MEDS_n, MEDS_m, i, N1_r);

  LOG_VEC(sol, MEDS_m*MEDS_m + MEDS_n*MEDS_n);

  // Perfrom back-substution for last block in top part.
  for (int c = N1_r - 1; c >= 0; c--)
    for (int r = 0; r < MEDS_n; r++)
    {
      uint64_t tmp1 = pmod_mat_entry(P00nt, MEDS_n, MEDS_m, r, c);
      uint64_t tmp2 = pmod_mat_entry(N1, MEDS_m-1, MEDS_m, c, N1_r);

      uint64_t prod = (tmp1 * tmp2) % MEDS_p;

      uint64_t val = sol[(MEDS_m-1)*MEDS_n + r];

      val = ((MEDS_p + val) - prod) % MEDS_p;

      sol[(MEDS_m-1)*MEDS_n + r] = val;
    }

  LOG_VEC(sol, MEDS_m*MEDS_m + MEDS_n*MEDS_n);


  // Perfrom back-substution for all remaining blocks in bottom part.
  for (int b = MEDS_m-3; b >= 0; b--)
    for (int c = MEDS_m - 1; c >= 0; c--)
      for (int r = 0; r < MEDS_m; r++)
      {
        uint64_t tmp1 = pmod_mat_entry(N, MEDS_n-1, MEDS_m, r, c);
        uint64_t tmp2 = sol[(MEDS_m+1)*MEDS_n + b*MEDS_m + MEDS_m + c];


        uint64_t prod = (tmp1 * tmp2) % MEDS_p;

        uint64_t val = sol[(MEDS_m+1)*MEDS_n + b*MEDS_m + r];

        val = ((MEDS_p + val) - prod) % MEDS_p;

        sol[(MEDS_m+1)*MEDS_n + b*MEDS_m + r] = val;
      }

  LOG_VEC(sol, MEDS_m*MEDS_m + MEDS_n*MEDS_n);


  // Perfrom back-substution for all remaining blocks in top part.
  for (int b = MEDS_m-2; b >= 0; b--)
    for (int c = MEDS_m - 1; c >= 0; c--)
      for (int r = 0; r < MEDS_n; r++)
      {
        uint64_t tmp1 = pmod_mat_entry(P00nt, MEDS_n, MEDS_m, r, c);
        uint64_t tmp2 = sol[2*MEDS_m*MEDS_n - (MEDS_m - 1) - (MEDS_m-1-b)*MEDS_m + c];


        uint64_t prod = (tmp1 * tmp2) % MEDS_p;

        uint64_t val = sol[b*MEDS_n + r];

        val = ((MEDS_p + val) - prod) % MEDS_p;

        sol[b*MEDS_n + r] = val;
      }

  LOG_VEC(sol, MEDS_m*MEDS_m + MEDS_n*MEDS_n);


  for (int i = 0; i < MEDS_m*MEDS_m; i++)
    A[i] = sol[i + MEDS_n*MEDS_n];

  for (int i = 0; i < MEDS_n*MEDS_n; i++)
    B_inv[i] = sol[i];

  LOG_MAT(A, MEDS_m, MEDS_m);
  LOG_MAT(B_inv, MEDS_n, MEDS_n);

  return 0;
}

/*int solve_opt(pmod_mat_t *A, pmod_mat_t *B_inv, pmod_mat_t *G0prime, GFq_t Amm)
{
  _Static_assert (MEDS_n == MEDS_m, "solve_opt requires MEDS_n == MEDS_m");

  pmod_mat_t P0prime0[MEDS_m*MEDS_n];
  pmod_mat_t P0prime1[MEDS_m*MEDS_n];

  for (int i = 0; i < MEDS_m*MEDS_n; i++)
  {
    P0prime0[i] = G0prime[i];
    P0prime1[i] = G0prime[i + MEDS_m * MEDS_n];
  }

  pmod_mat_t N[MEDS_n * MEDS_m];

  for (int i = 0; i < MEDS_m; i++)
    for (int j = 0; j < MEDS_n; j++)
      N[j*MEDS_m + i] = (MEDS_p - P0prime0[i*MEDS_n + j]) % MEDS_p;

  //LOG_MAT(N, MEDS_n, MEDS_m);


  pmod_mat_t M[MEDS_n*(MEDS_m + MEDS_m + 2)] = {0};

  for (int i = 0; i < MEDS_m; i++)
    for (int j = 0; j < MEDS_n; j++)
      M[j*(MEDS_m + MEDS_m + 2) + i] = (MEDS_p - P0prime1[i*MEDS_n + j]) % MEDS_p;

  for (int i = 0; i < MEDS_m; i++)
    for (int j = 0; j < MEDS_n; j++)
      M[j*(MEDS_m + MEDS_m + 2) + i + MEDS_n] = P0prime0[i*MEDS_n + j];

  for (int j = 0; j < MEDS_n; j++)
    M[j*(MEDS_m + MEDS_m + 2) + MEDS_m + MEDS_n] = ((uint32_t)P0prime0[(MEDS_m-1)*MEDS_n + j] * (MEDS_p - (uint32_t)Amm)) % MEDS_p;

  for (int j = 0; j < MEDS_n; j++)
    M[j*(MEDS_m + MEDS_m + 2) + MEDS_m + MEDS_n + 1] = ((uint32_t)P0prime1[(MEDS_m-1)*MEDS_n + j] * (uint32_t)Amm) % MEDS_p;


  //LOG_MAT(M, MEDS_n, MEDS_m + MEDS_m + 2);

  if (pmod_mat_syst_ct(M, MEDS_n-1, MEDS_m + MEDS_m + 2) < 0)
    return -1;

  //LOG_MAT_FMT(M, MEDS_n, MEDS_m + MEDS_m + 2, "M part");

  // eliminate last row
  for (int r = 0; r < MEDS_n-1; r++)
  {
    uint64_t factor = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n-1, r);

    // ignore last column
    for (int c = MEDS_n-1; c < MEDS_m + MEDS_m + 1; c++)
    {
      uint64_t tmp0 = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n-1, c);
      uint64_t tmp1 = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, r, c);

      int64_t val = (tmp1 * factor) % MEDS_p;

      val = tmp0 - val;

      val += MEDS_p * (val < 0);

      pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2,  MEDS_n-1, c, val);
    }

    pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n-1, r, 0);
  }

  // normalize last row
  {
    uint64_t val = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n-1, MEDS_n-1);

    if (val == 0)
      return -1;

    val = GF_inv(val);

    // ignore last column
    for (int c = MEDS_n; c < MEDS_m + MEDS_m + 1; c++)
    {
      uint64_t tmp = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n-1, c);

      tmp = (tmp * val) % MEDS_p;

      pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n-1, c, tmp);
    }
  }

  pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n-1, MEDS_n-1, 1);

  M[MEDS_n*(MEDS_m + MEDS_m + 2)-1] = 0;

  //LOG_MAT_FMT(M, MEDS_n, MEDS_m + MEDS_m + 2, "M red");

  // back substitute
  for (int r = 0; r < MEDS_n-1; r++)
  {
    uint64_t factor = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, r, MEDS_n-1);

    // ignore last column
    for (int c = MEDS_n; c < MEDS_m + MEDS_m + 1; c++)
    {
        uint64_t tmp0 = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, MEDS_n-1, c);
        uint64_t tmp1 = pmod_mat_entry(M, MEDS_n, MEDS_m + MEDS_m + 2, r, c);

        int64_t val = (tmp0 * factor) % MEDS_p;

        val = tmp1 - val;

        val += MEDS_p * (val < 0);

        pmod_mat_set_entry(M, MEDS_n, MEDS_m + MEDS_m + 2,  r, c, val);
    }

    pmod_mat_set_entry(M, M_r, MEDS_m + MEDS_m + 2, r, MEDS_n-1, 0);
  }


  //LOG_MAT_FMT(M, MEDS_n, MEDS_m + MEDS_m + 2, "M done");


  GFq_t sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m] = {0};

  sol[MEDS_n*MEDS_n + MEDS_m * MEDS_m - 1] = Amm;

  for (int i = 0; i < MEDS_n - 1; i++)
    sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - MEDS_n + i] = M[(i+1) * (MEDS_m + MEDS_m + 2) - 1];

  for (int i = 0; i < MEDS_n; i++)
    sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - 2*MEDS_n + i] = M[(i+1) * (MEDS_m + MEDS_m + 2) - 2];

  for (int i = 0; i < MEDS_n; i++)
    sol[MEDS_n*MEDS_n - MEDS_n + i] = ((uint32_t)P0prime0[(MEDS_m-1)*MEDS_n + i] * (uint32_t)Amm) % MEDS_p;

  //LOG_VEC_FMT(sol, MEDS_n*MEDS_n + MEDS_m*MEDS_m, "initial sol");


  // incomplete blocks:

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n-1; j++)
    {
      uint32_t tmp = (uint32_t)sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - 2*MEDS_n + i] + MEDS_p - 
        ((uint32_t)M[i * (MEDS_m + MEDS_m + 2) + MEDS_n + MEDS_n-2 - j] * (uint32_t)sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - 2 - j]) % MEDS_p;
      sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - 2*MEDS_n + i] = tmp % MEDS_p;
    }

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n-1; j++)
    {
      uint32_t tmp = (uint32_t)sol[MEDS_n*MEDS_n - MEDS_n + i] + MEDS_p - 
        ((uint32_t)N[i * (MEDS_n) + MEDS_m-2 - j] * (uint32_t)sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - 2 - j]) % MEDS_p;
      sol[MEDS_n*MEDS_n - MEDS_n + i] = tmp % MEDS_p;
    }

  //LOG_VEC_FMT(sol, MEDS_n*MEDS_n + MEDS_m*MEDS_m, "incomplete blocks");


  // complete blocks:

  for (int block = 3; block <= MEDS_n; block++)
    for (int i = 0; i < MEDS_n; i++)
      for (int j = 0; j < MEDS_n; j++)
      {
        uint32_t tmp = sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - block*MEDS_n + i] + MEDS_p - 
          ((uint32_t)M[i * (MEDS_m + MEDS_m + 2) + MEDS_n + MEDS_n-1 - j] * (uint32_t)sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - 1 - (block-2)*MEDS_n - j]) % MEDS_p;
        sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - block*MEDS_n + i] = tmp % MEDS_p;
      }

  for (int block = 2; block <= MEDS_n; block++)
    for (int i = 0; i < MEDS_n; i++)
      for (int j = 0; j < MEDS_n; j++)
      {
        uint32_t tmp = sol[MEDS_n*MEDS_n - block*MEDS_n + i] + MEDS_p - 
          ((uint32_t)N[i * (MEDS_n) + MEDS_m-1 - j] * (uint32_t)sol[MEDS_n*MEDS_n + MEDS_m*MEDS_m - 1 - (block-1)*MEDS_n - j]) % MEDS_p;
        sol[MEDS_n*MEDS_n - block*MEDS_n + i] = tmp % MEDS_p;
      }

  //LOG_VEC_FMT(sol, MEDS_n*MEDS_n + MEDS_m*MEDS_m, "complete blocks");


  for (int i = 0; i < MEDS_m * MEDS_m; i++)
    A[i] = sol[i + MEDS_n * MEDS_n];

  for (int i = 0; i < MEDS_n * MEDS_n; i++)
    B_inv[i] = sol[i];

  //LOG_MAT(A, MEDS_m, MEDS_m);
  //LOG_MAT(B_inv, MEDS_n, MEDS_n);

  return 0;
}
*/

void pi(pmod_mat_t *Gout, pmod_mat_t *A, pmod_mat_t *B, pmod_mat_t *G)
{
  for (int i = 0; i < MEDS_k; i++)
  {
    pmod_mat_mul(&Gout[i*MEDS_m*MEDS_n], MEDS_m, MEDS_n, A, MEDS_m, MEDS_m, &G[i*MEDS_m*MEDS_n], MEDS_m, MEDS_n);
    pmod_mat_mul(&Gout[i*MEDS_m*MEDS_n], MEDS_m, MEDS_n, &Gout[i*MEDS_m*MEDS_n], MEDS_m, MEDS_n, B, MEDS_n, MEDS_n);
  }
}

int SF(pmod_mat_t *Gprime, pmod_mat_t *G)
{
  pmod_mat_t M[MEDS_k * MEDS_k];

  for (int r = 0; r < MEDS_k; r++)
    memcpy(&M[r * MEDS_k], &G[r * MEDS_m * MEDS_n], MEDS_k * sizeof(GFq_t));

  if (pmod_mat_inv(M, M, MEDS_k, MEDS_k) == 0)
  {
     pmod_mat_mul(Gprime, MEDS_k, MEDS_m * MEDS_n, M, MEDS_k, MEDS_k, G, MEDS_k, MEDS_m * MEDS_n);

     return 0;
  }

  return -1;
}

