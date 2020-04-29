/*
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef HEAP_H
#define HEAP_H

#include "flint.h"

#define MIN_HEAP_CAP 32

typedef struct 
{
    slong *idx; /* Heap of indices */
    slong *pos; /* Inverse of idx (map from index to heap position) */
    slong *val; /* Map from index to values */
    slong num; /* Number of elements in heap */
    slong max; /* Largest index in heap */
    slong cap; /* Capacity of heap */
} heap_struct;

typedef heap_struct heap_t[1];

/* Initialize heap with initial capacity cap */
void heap_init(heap_t h, slong cap);

/* Free heap contents */
void heap_clear(heap_t h);

/* Push new score onto heap, returning assigned index */
/* Note: if want to re-insert a current element with */
/* new score, use heap_adjust instead */
slong heap_push(heap_t h, slong val);

/* Pop element with minimal score from heap */
slong heap_pop(heap_t h, slong *val);

/* Modify score of existing element in heap */
slong heap_adjust(heap_t h, slong idx, slong val);

#endif