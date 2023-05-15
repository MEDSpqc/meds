#include <stdio.h>
#include <string.h>

#include "fips202.h"

#ifndef MEDS_t
#include "params.h"
#endif

#include "clog.h"
#include "seed.h"

void print_tree(uint8_t *stree)
{
  int h = 0;
  int i = 0;

  int start = i;
  int end = i+1;

  for (; h < CLOG2(MEDS_t-1)+1; h++)
  {
    for (int i = 0; i < (1 << (CLOG2(MEDS_t-1) - h))-1; i++)
      printf("  ");

    for (int i = start; i < end; i++)
    {
      // print incomplete tree for non-power-of-two number of leafs
      if ((i << (CLOG2((MEDS_t - 1)) - h)) >= MEDS_t)
        break;

      printf("%02x", stree[MEDS_sec_bytes * TREE_ADR(h, i)]);

      for (int j = 0; j < ((1 << (CLOG2(MEDS_t-1) - h)) - 1)*2 + 1; j++)
        printf("  ");
    }

    printf("\n");

    start = start<<1;
    end = end<<1;
  }
}

void t_hash(uint8_t *stree, uint8_t *salt, int h, int i)
{
  keccak_state shake;

  int start = i;
  int end = i+1;

  uint8_t buf[MEDS_salt_sec_bytes + MEDS_sec_bytes + sizeof(uint32_t)];

  memcpy(buf, salt, MEDS_salt_sec_bytes);

  uint32_t *pos; pos = (uint32_t*)(buf + MEDS_salt_sec_bytes + MEDS_sec_bytes);

  for (h = h+1; h < CLOG2(MEDS_t-1)+1; h++)
  {
    start = start<<1;
    end = end<<1;

    if ((start << (CLOG2((MEDS_t - 1)) - h)) >= MEDS_t)
      break;

    for (int i = start; i < end; i+=2)
    {
      if ((i << (CLOG2((MEDS_t - 1)) - h)) >= MEDS_t)
        break;

      *pos = TREE_ADR(h-1, i>>1);

      memcpy(buf + MEDS_salt_sec_bytes,
          &stree[MEDS_sec_bytes * (*pos)],
          MEDS_sec_bytes);

      shake256_init(&shake);
      shake256_absorb(&shake, buf, MEDS_salt_sec_bytes + MEDS_sec_bytes + sizeof(uint32_t));
      shake256_finalize(&shake);

      int len = 2*MEDS_sec_bytes;

      if (((i+1) << (CLOG2((MEDS_t - 1)) - h)) >= MEDS_t)
        len = MEDS_sec_bytes;

      shake256_squeeze(&stree[MEDS_sec_bytes * TREE_ADR(h, i)], len, &shake);
    }
  }
}

void stree_delete_path(uint8_t *stree, int index)
{
  for (int h = 0; h < CLOG2(MEDS_t-1)+1; h++)
    memset(&stree[MEDS_sec_bytes * TREE_ADR(h, index >> (CLOG2((MEDS_t - 1)) - h))], 0, MEDS_sec_bytes);
}

void stree_to_path_to_stree(uint8_t *stree, unsigned int *indices, uint8_t *path, uint32_t *empty_slots, int empty_slot_size, uint8_t *salt, int mode)
{
  int h = 0;
  int i = 0;

  int offset = 0;

  unsigned int id = 0;

  while (1 == 1)
  {
    int index_leaf = 0;

    while ((indices[id] >> (CLOG2((MEDS_t - 1)) - h)) == i)
    {
      h += 1;
      i <<= 1;

      if (h > CLOG2((MEDS_t-1)))
      {
        h -= 1;
        i >>= 1;

        *(empty_slots++) = offset;

        if (id+1 < MEDS_w)
          id++;

        //printf("%i %i %02x -- hit end\n", h, i, stree[MEDS_sec_bytes * TREE_ADR(h, i)]);

        offset += empty_slot_size;

        index_leaf = 1;

        break;
      }
    }

    if (index_leaf == 0)
    {
      //printf("%i %i %02x\n", h, i, stree[MEDS_sec_bytes * TREE_ADR(h, i)]);

      if (mode == STREE_TO_PATH)
        for (int j = 0; j < MEDS_sec_bytes; j++)
          path[offset++] = stree[MEDS_sec_bytes * TREE_ADR(h, i) + j];
      else
      {
        for (int j = 0; j < MEDS_sec_bytes; j++)
          stree[MEDS_sec_bytes * TREE_ADR(h, i) + j] = path[offset++];

        //printf("%i %i %02x\n", h, i, stree[MEDS_sec_bytes * TREE_ADR(h, i)]);

        t_hash(stree, salt, h, i);
      }
    }

    // backtrack
    while ((i & 1) == 1)
    {
      h -= 1;
      i >>= 1;
    }

    // next leaf
    i+=1;

    if ((i << (CLOG2((MEDS_t - 1)) - h)) >= MEDS_t)
      return;
  }
}

