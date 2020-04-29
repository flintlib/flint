/*
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "flint.h"
#include "heap.h"

static void _heap_up(heap_t h, slong pos)
{
    const slong idx = h->idx[pos];
    slong nidx, npos;
    for (; pos > 0; pos = npos)
    {
        npos = (pos-1)/2;
        nidx = h->idx[npos];
        if (h->val[idx] >= h->val[nidx]) break;

        h->idx[pos] = nidx;
        h->pos[nidx] = pos;
    }
    h->idx[pos] = idx;
    h->pos[idx] = pos;
}

void _heap_down(heap_t h, slong pos)
{
    const slong idx = h->idx[pos];
    slong nidx, npos;
    for (; (npos = 2*pos+1) < h->num; pos = npos)
    {
        if (npos+1 < h->num && h->val[h->idx[npos]] > h->val[h->idx[npos+1]]) ++npos;
        nidx = h->idx[npos];
        if (h->val[idx] <= h->val[nidx]) break;

        h->idx[pos] = nidx;
        h->pos[nidx] = pos;
    }
    h->idx[pos] = idx;
    h->pos[idx] = pos;
}

void heap_init(heap_t h, slong cap)
{
    h->cap = FLINT_MAX(cap, MIN_HEAP_CAP);
    h->idx = flint_malloc(h->cap*sizeof(h->idx));
    h->pos = flint_malloc(h->cap*sizeof(h->pos));
    h->val = flint_malloc(h->cap*sizeof(h->val));
    h->num = h->max = 0;
    h->cap = cap;
}

void heap_clear(heap_t h)
{
    flint_free(h->idx);
    flint_free(h->pos);
    flint_free(h->val);
    memset(h, 0, sizeof(*h));
}

slong heap_push(heap_t h, slong val)
{
    if(h->max == h->cap)
    {
        h->cap *= 2;
        h->idx = flint_realloc(h->idx, h->cap*sizeof(h->idx));
        h->pos = flint_realloc(h->pos, h->cap*sizeof(h->pos));
        h->val = flint_realloc(h->val, h->cap*sizeof(h->val));
    }
    h->idx[h->max] = h->max;
    h->val[h->max] = val;
    ++h->num;
    _heap_up(h, h->max);
    return h->max++;
}

slong heap_pop(heap_t h, slong *val)
{
    slong idx = h->idx[0];
    if(val) *val = h->val[idx];
    h->pos[idx] = -1;
    h->idx[0] = h->idx[--h->num];
    _heap_down(h, 0);
    return idx;
}

slong heap_adjust(heap_t h, slong idx, slong val)
{
    slong oval = h->val[idx], pos = h->pos[idx];
    h->val[idx] = val;
    if (pos != -1)
    {
        if (oval < val) _heap_down(h, pos);
        if (oval > val) _heap_up(h, pos);
    }
    return oval;
}
