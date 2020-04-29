/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef HASHMAP_H
#define HASHMAP_H

#include "flint.h"

#define MIN_HASHMAP_BITS 5

typedef struct
{
   ulong size;
   ulong mask;
   ulong num;
   ulong *table;
   slong *keys;
   void **vals;
} hashmap_struct;

typedef hashmap_struct hashmap_t[1];

/* Initialize hash table to accomodate size elements */
void hashmap_init(hashmap_t h, slong size);

/* Clear hash table */
void hashmap_clear(hashmap_t h);

/* Get value associated with given key */
void * hashmap_get(hashmap_t h, slong key);

/* Assign value to a given key */
void hashmap_put(hashmap_t h, slong key, void *val);

/* Remove value with a given key */
void hashmap_rem(hashmap_t h, slong key);

#endif
