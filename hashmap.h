/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef HASHMAP_H
#define HASHMAP_H

#include "flint.h"

#define HASHMAP_START_SIZE 1024
#define HASHMAP_START_MASK (HASHMAP_START_SIZE - 1)

/******************************************************************************

   Hashmap types with one word key

******************************************************************************/

typedef struct hashmap1_elem_s
{
   ulong key;
   void * value;
   int in_use;
} hashmap1_elem_s;

typedef struct hashmap1_s
{
   slong alloc;
   slong num_used;
   ulong mask;
   hashmap1_elem_s * data;
} hashmap1_s;

typedef hashmap1_s hashmap1_t[1];

/******************************************************************************

   Hash functions

******************************************************************************/

/* from lookup3.c, by Bob Jenkins, May 2006, Public Domain. */

#define hash_rot(x, k) (((x) << (k)) | ((x) >> (32 - (k))))

#define hash_mix(a,b,c) \
{ \
  c ^= b; c -= hash_rot(b, 14); \
  a ^= c; a -= hash_rot(c, 11); \
  b ^= a; b -= hash_rot(a, 25); \
  c ^= b; c -= hash_rot(b, 16); \
  a ^= c; a -= hash_rot(c,  4); \
  b ^= a; b -= hash_rot(a, 14); \
  c ^= b; c -= hash_rot(b, 24); \
}

/* End of Public Domain code. */

#if FLINT64

static __inline__
ulong hash_word(ulong val)
{
   int * ptr = (int * ) &val;
   int a = ptr[0], b = ptr[1], c = 0;

   hash_mix(a, b, c);
   
   ptr[0] = b;
   ptr[1] = c;

   return val;
}

#else

static __inline__
ulong hash_word(ulong a)
{
   int b = 0, c = 0;

   hash_mix(a, b, c);
   
   return c;
}

#endif

/******************************************************************************

   Hashmap functions with one word key

******************************************************************************/

FLINT_DLL void hashmap1_init(hashmap1_t h);

FLINT_DLL void hashmap1_init2(hashmap1_t h, slong size);

FLINT_DLL void hashmap1_clear(hashmap1_t h);

static __inline__
ulong hashmap1_hash_key(ulong key, hashmap1_t h)
{
   return hash_word(key) & h->mask;
}

FLINT_DLL slong hashmap1_hash(ulong key, hashmap1_t h);

FLINT_DLL void hashmap1_rehash(hashmap1_t h);

FLINT_DLL void hashmap1_insert(ulong key, void * value, hashmap1_t h);

FLINT_DLL int hashmap1_find(void ** ptr, ulong key, hashmap1_t h);

#endif
