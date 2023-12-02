/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>

#include "flint.h"
#include "ulong_extras.h"

typedef struct apow {
    ulong k;
    ulong ak;
} apow_t;

typedef struct {
    ulong n;
    double ninv;
    ulong m;
    ulong am;
    apow_t * table;
} bsgs_struct;

typedef bsgs_struct bsgs_t[1];

static int
apow_cmp(const apow_t * x, const apow_t * y)
{
    return (x->ak < y->ak) ? -1 : (x->ak > y->ak);
}

/* set size of table m=sqrt(nk) to compute k logs in a group of size n */
void
bsgs_table_init(bsgs_t t, ulong a, ulong n, ulong m)
{
    ulong k, ak;
    t->table = (apow_t *)flint_malloc(m * sizeof(apow_t));

    t->n = n;
    t->ninv = n_precompute_inverse(n);
    t->m = m;

    for (k = 0, ak = 1; k < m; k++)
    {
        t->table[k].k = k;
        t->table[k].ak = ak;
        ak = n_mulmod_precomp(ak, a, n, t->ninv);
    }

    t->am = n_invmod(ak, n);
    qsort(t->table, m, sizeof(apow_t), (int(*)(const void*,const void*))apow_cmp);
}

void
bsgs_table_clear(bsgs_t t)
{
    flint_free(t->table);
}

ulong
n_discrete_log_bsgs_table(const bsgs_t t, ulong b)
{
    ulong i;
    apow_t c, * x;

    c.k = 0;
    c.ak = b;
    for (i = 0; i < t->m; i++)
    {
        x = bsearch(&c, t->table, t->m, sizeof(apow_t),
            (int(*)(const void*,const void*))apow_cmp);
        if (x != NULL)
            return i * t->m + x->k;
        c.ak = n_mulmod_precomp(c.ak, t->am, t->n, t->ninv);
    }

    flint_throw(FLINT_ERROR, "Exception (n_discrete_log_bsgs).  discrete log not found.\n");
}

ulong
n_discrete_log_bsgs(ulong b, ulong a, ulong n)
{
    ulong m;
    bsgs_t table;

    m = ceil(sqrt((double) n));
    bsgs_table_init(table, a, n, m);
    m = n_discrete_log_bsgs_table(table, b);
    bsgs_table_clear(table);

    return m;
}
