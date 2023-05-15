#ifndef SEED_H
#define SEED_H

//#ifndef CLOG2
//
//#define CLOG22(n) ((n&2)?1:0)
//#define CLOG24(n) ((n&(0xC))?(2+CLOG22(n>>2)):(CLOG22(n)))
//#define CLOG28(n) ((n&0xF0)?(4+CLOG24(n>>4)):(CLOG24(n)))
//#define CLOG216(n) ((n&0xFF00)?(8+CLOG28(n>>8)):(CLOG28(n)))
//#define CLOG232(n) ((n&0xFFFF0000)?(16+CLOG216(n>>16)):(CLOG216(n)))
//#define CLOG2(n) ((n)==0?0:CLOG232((n))+1)
//
//#endif


#define TREE_ADR(h, i) ((1 << (h+1)) - 2 + (i))

void print_tree(uint8_t *stree);

void t_hash(uint8_t *stree, uint8_t *salt, int h, int i);

void stree_delete_path(uint8_t *stree, int index);

#define STREE_TO_PATH 0
#define PATH_TO_STREE 1

void stree_to_path_to_stree(uint8_t *stree, unsigned int *indices, uint8_t *path, uint32_t *empty_slots, int empty_slot_size, uint8_t *salt, int mode);

#define stree_to_path(stree, indices, path, empty_slots, empty_slot_size, salt) stree_to_path_to_stree(stree, indices, path, empty_slots, empty_slot_size, salt, STREE_TO_PATH)

#define path_to_stree(stree, indices, path, empty_slots, empty_slot_size, salt) stree_to_path_to_stree(stree, indices, path, empty_slots, empty_slot_size, salt, PATH_TO_STREE)

#endif

