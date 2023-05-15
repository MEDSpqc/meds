#include <stdint.h>
#include <string.h>

#include "fips202.h"

#include "util.h"
#include "meds.h"
#include "matrixmod.h"

//#define DEBUG

void G_mat_init(pmod_mat_t *G, pmod_mat_t *Gsub[MEDS_k])
{
  for (int i = 0; i < MEDS_k; i++)
    Gsub[i] = G + i*MEDS_m*MEDS_n;
}

GFq_t rnd_GF(keccak_state *shake)
{
  GFq_t val = MEDS_p;

  while (val >= MEDS_p)
  {
    shake256_squeeze((void*)&val, sizeof(GFq_t), shake);

    val = val & ((1 << GFq_bits) - 1);

    // uint8_t data[sizeof(GFq_t)];
    //
    //shake256_squeeze(data, sizeof(GFq_t), shake);

    //for (int i = 0; i < sizeof(GFq_t); i++)
    //  val = (val << 8) | data[i];
  }

  return val;
}

GFq_t rnd_GF_min(keccak_state *shake, int min)
{
  GFq_t val = MEDS_p;

  while ((val < min) || (val >= MEDS_p))
  {
    shake256_squeeze((void*)&val, sizeof(GFq_t), shake);

    val = val & ((1 << GFq_bits) - 1);

    // uint8_t data[sizeof(GFq_t)];
    //
    //shake256_squeeze(data, sizeof(GFq_t), shake);

    //for (int i = 0; i < sizeof(GFq_t); i++)
    //  val = (val << 8) | data[i];
  }

  return val;
}

void sys_mat_from_bytes(pmod_mat_t *M, int M_r, int M_c, uint8_t *data, int data_len)
{
  for (int r = 0; r < M_r; r++)
    for (int c = M_r; c < M_c; c++)
    {
      pmod_mat_set_entry(M, M_r, M_c, r, c, *((GFq_t*)data));

      data += sizeof(GFq_t);
    }

  for (int r = 0; r < M_r; r++)
    for (int c = 0; c < M_r; c++)
      if (r == c)
        pmod_mat_set_entry(M, M_r, M_c, r, c, 1);
      else
        pmod_mat_set_entry(M, M_r, M_c, r, c, 0);
}

void rnd_sys_mat(pmod_mat_t *M, int M_r, int M_c, keccak_state *shake)
{
  for (int r = 0; r < M_r; r++)
    for (int c = M_r; c < M_c; c++)
      pmod_mat_set_entry(M, M_r, M_c, r, c, rnd_GF(shake));

  for (int r = 0; r < M_r; r++)
    for (int c = 0; c < M_r; c++)
      if (r == c)
        pmod_mat_set_entry(M, M_r, M_c, r, c, 1);
      else
        pmod_mat_set_entry(M, M_r, M_c, r, c, 0);
}

void rnd_matrix(pmod_mat_t *M, int M_r, int M_c, keccak_state *shake)
{
  for (int r = 0; r < M_r; r++)
    for (int c = 0; c < M_c; c++)
      pmod_mat_set_entry(M, M_r, M_c, r, c, rnd_GF(shake));
}

void rnd_invtl_matrix(pmod_mat_t *M, int M_r, int M_c, keccak_state *shake)
{
  pmod_mat_t L[M_r*M_c];

  for (int i = 0; i < M_r; i++)
  {
    for (int j = 0; j < i; j++)
      pmod_mat_set_entry(L, M_r, M_c, j, i, 0);

    pmod_mat_set_entry(L, M_r, M_c, i, i, 1);

    for (int j = i+1; j < M_c; j++)
      pmod_mat_set_entry(L, M_r, M_c, j, i, rnd_GF(shake));
  }

  pmod_mat_t U[M_r*M_c];

  for (int i = 0; i < M_r; i++)
  {
    for (int j = 0; j < i; j++)
      pmod_mat_set_entry(U, M_r, M_c, i, j, 0);

    // get only non-zero values for the diagonal
    pmod_mat_set_entry(U, M_r, M_c, i, i, rnd_GF_min(shake, 1));

    for (int j = i+1; j < M_c; j++)
      pmod_mat_set_entry(U, M_r, M_c, i, j, rnd_GF(shake));
  }

  for (int c = 0; c < M_c; c++)
    for (int r = 0; r < M_r; r++)
    {
      int val = 0;

      //for (int i = 0; i <= c; i++)
      for (int i = 0; i < M_r; i++)
        val += pmod_mat_entry(L, M_r, M_c, r, i) * pmod_mat_entry(U, M_r, M_c, i, c);

      M[r*M_c + c] = val % MEDS_p;
    }

}

void rnd_inv_matrix(pmod_mat_t *M, int M_r, int M_c, pmod_mat_t *M_inv, keccak_state *shake)
{
  if (M_inv != NULL)
    while (0==0)
    {
      rnd_matrix(M, M_r, M_c, shake);

      if (pmod_mat_inv(M_inv, M, M_r, M_c) == 0)
        return;
    }

  pmod_mat_t L[M_r*M_c];

  for (int i = 0; i < M_r; i++)
  {
    for (int j = 0; j < i; j++)
      pmod_mat_set_entry(L, M_r, M_c, j, i, 0);

    pmod_mat_set_entry(L, M_r, M_c, i, i, 1);

    for (int j = i+1; j < M_c; j++)
      pmod_mat_set_entry(L, M_r, M_c, j, i, rnd_GF(shake));
  }

#ifdef DEBUG
  fprintf(stderr, "(%s) L:\n", __func__);
  pmod_mat_fprint(stderr, L, M_r, M_c);
  fprintf(stderr, "\n");
#endif //DEBUG

  pmod_mat_t U[M_r*M_c];

  for (int i = 0; i < M_r; i++)
  {
    for (int j = 0; j < i; j++)
      pmod_mat_set_entry(U, M_r, M_c, i, j, 0);

    // get only non-zero values for the diagonal
    pmod_mat_set_entry(U, M_r, M_c, i, i, rnd_GF_min(shake, 1));

    for (int j = i+1; j < M_c; j++)
      pmod_mat_set_entry(U, M_r, M_c, i, j, rnd_GF(shake));
  }

#ifdef DEBUG
  fprintf(stderr, "(%s) U:\n", __func__);
  pmod_mat_fprint(stderr, U, M_r, M_c);
  fprintf(stderr, "\n");
#endif //DEBUG

  //pmod_mat_mul(M, M_r, M_c, L, M_r, M_c, U, M_r, M_c);

  for (int c = 0; c < M_c; c++)
    for (int r = 0; r < M_r; r++)
    {
      uint64_t val = 0;

      //for (int i = 0; i <= c; i++)
      for (int i = 0; i < M_r; i++)
        val = (val + 
            (uint64_t)pmod_mat_entry(L, M_r, M_c, r, i) * 
            (uint64_t)pmod_mat_entry(U, M_r, M_c, i, c)
            ) % MEDS_p;

      M[r*M_c + c] = val; // % MEDS_p;
    }

  if (M_inv == NULL)
    return;
  {
    //pmod_mat_inv(L, L, M_r, M_c);
    {
      pmod_mat_t M[M_r * M_c*2];

      for (int r = 0; r < M_r; r++)
      {
        memcpy(&M[r * M_c*2], &L[r * M_c], M_c * sizeof(GFq_t));

        for (int c = 0; c < M_c; c++)
          pmod_mat_set_entry(M, M_r, M_c*2, r, M_c + c, r==c ? 1 : 0);
      }

      for (int r = M_r-2; r >= 0; r--)
      {
        // eliminate
        for (int r2 = r+1; r2 < M_r; r2++)
        {
          uint64_t factor = pmod_mat_entry(M, M_r, M_c*2, r2, r);

          for (int c = 0; c < M_c*2; c++)
          {
            uint64_t tmp0 = pmod_mat_entry(M, M_r, M_c*2, r, c);
            uint64_t tmp1 = pmod_mat_entry(M, M_r, M_c*2, r2, c);

            int64_t val = (tmp0 * factor) % MEDS_p;

            val = tmp1 - val;

            val += MEDS_p * (val < 0);

            pmod_mat_set_entry(M, M_r, M_c*2, r2, c, val);
          }
        }
      }

      for (int r = 0; r < M_r; r++)
        memcpy(&L[r * M_c], &M[r * M_c*2 + M_c], M_c * sizeof(GFq_t));
    }

    //pmod_mat_inv(U, U, M_r, M_c);
    {
      pmod_mat_t M[M_r * M_c*2];

      for (int r = 0; r < M_r; r++)
      {
        memcpy(&M[r * M_c*2], &U[r * M_c], M_c * sizeof(GFq_t));

        for (int c = 0; c < M_c; c++)
          pmod_mat_set_entry(M, M_r, M_c*2, r, M_c + c, r==c ? 1 : 0);
      }

      for (int r = 0; r < M_r; r++)
      {
        uint64_t val = pmod_mat_entry(M, M_r, M_c*2, r, r);

        val = GF_inv(val);

        // normalize
        for (int c = r; c < M_c*2; c++)
        {
          uint64_t tmp = ((uint64_t)pmod_mat_entry(M, M_r, M_c*2, r, c) * val) % MEDS_p;
          pmod_mat_set_entry(M, M_r, M_c*2, r, c, tmp);
        }
      }

      // back substitution
      for (int r = M_r - 1; r >= 0; r--)
        for (int r2 = 0; r2 < r; r2++)
        {
          uint64_t factor = pmod_mat_entry(M, M_r, M_c*2, r2, r);

          uint64_t tmp0 = pmod_mat_entry(M, M_r, M_c*2, r, r);
          uint64_t tmp1 = pmod_mat_entry(M, M_r, M_c*2, r2, r);

          int64_t val = (tmp0 * factor) % MEDS_p;

          val = tmp1 - val;

          val += MEDS_p * (val < 0);

          pmod_mat_set_entry(M, M_r, M_c*2,  r2, r, val);

          for (int c = M_r; c < M_c*2; c++)
          {
            uint64_t tmp0 = pmod_mat_entry(M, M_r, M_c*2, r, c);
            uint64_t tmp1 = pmod_mat_entry(M, M_r, M_c*2, r2, c);

            int64_t val = (tmp0 * factor) % MEDS_p;

            val = tmp1 - val;

            val += MEDS_p * (val < 0);

            pmod_mat_set_entry(M, M_r, M_c*2,  r2, c, val);
          }
        }

      //pmod_mat_syst_ct(M, M_r, M_c*2);

      for (int r = 0; r < M_r; r++)
        memcpy(&U[r * M_c], &M[r * M_c*2 + M_c], M_c * sizeof(GFq_t));
    }

    pmod_mat_mul(M_inv, M_r, M_c, U, M_r, M_c, L, M_r, M_c);
  }
}

int parse_hash(uint8_t *digest, int digest_len, uint8_t *h, int len_h)
{
  if (len_h < MEDS_t)
    return -1;

#ifdef DEBUG
  fprintf(stderr, "(%s) digest: [", __func__);
  for (int i = 0; i < MEDS_sec_bytes; i++)
  {
    fprintf(stderr, "%i", digest[i]);
    if (i < MEDS_sec_bytes-1)
      fprintf(stderr, ", ");
  }
  fprintf(stderr, "]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "(%s) digest len: %i\n", __func__, digest_len);
  fprintf(stderr, "\n");
#endif

  keccak_state shake;

  shake256_init(&shake);
  shake256_absorb(&shake, digest, digest_len);
  shake256_finalize(&shake);

  uint64_t pos[MEDS_w];
  uint64_t val[MEDS_w];

  int len_pos = 0;

  while (len_pos < MEDS_w)
  {
    uint64_t tmp_pos;

    shake256_squeeze((uint8_t*)&tmp_pos, CLOG8(MEDS_t-1), &shake);

    tmp_pos = tmp_pos & MEDS_t_mask;

    if (tmp_pos >= MEDS_t)
      continue;

    int j;

    for (j = 0; j < len_pos; j++)
      if (pos[j] == tmp_pos)
        break;

    if (j < len_pos)
      continue;

#ifdef DEBUG
    fprintf(stderr, "(%s) pos: %lu\n", __func__, tmp_pos);
    fprintf(stderr, "\n");
#endif

    uint8_t tmp_val = 0;

    while ((tmp_val < 1) || (tmp_val > MEDS_s-1))
    {
      shake256_squeeze(&tmp_val, 1, &shake);
      tmp_val = tmp_val & MEDS_s_mask;
//#ifdef DEBUG
//      fprintf(stderr, "1 %i %i\n", MEDS_s_mask, tmp_val);
//#endif
    }

    pos[len_pos] = tmp_pos;
    val[len_pos] = tmp_val;

    len_pos += 1;
  }

  for (int i = 0; i < MEDS_t; i++)
    h[i] = 0;

  for (int i = 0; i < MEDS_w; i++)
  {
#ifdef DEBUG
    fprintf(stderr, "(%s) p: %lu  v: %lu\n", __func__, pos[i], val[i]);
    fprintf(stderr, "\n");
#endif
    h[pos[i]] = val[i];
  }

  return 0;
}

