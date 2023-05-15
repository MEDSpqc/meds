#ifndef UTILS_H
#define UTILS_H

#include "params.h"
#include "matrixmod.h"

void G_mat_init(pmod_mat_t *G, pmod_mat_t *Gsub[MEDS_k]);

GFq_t rnd_GF(keccak_state *shake);

GFq_t rnd_GF_min(keccak_state *shake, int min);

void sys_mat_from_bytes(pmod_mat_t *M, int M_r, int M_c, uint8_t *data, int data_len);

void rnd_sys_mat(pmod_mat_t *M, int M_r, int M_c, keccak_state *shake);

void rnd_matrix(pmod_mat_t *M, int M_r, int M_c, keccak_state *shake);

void rnd_inv_matrix(pmod_mat_t *M, int M_r, int M_c, pmod_mat_t *M_inv, keccak_state *shake);

int parse_hash(uint8_t *digest, int digest_len, uint8_t *h, int len_h);

#endif

