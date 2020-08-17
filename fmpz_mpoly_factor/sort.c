/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz_mpoly_factor.h"

typedef struct {
    slong idx;
    fmpz exp;
    const fmpz_mpoly_struct * polys;
    const fmpz_mpoly_ctx_struct * ctx;
} sort_struct;

static int _sort(const void * a_, const void * b_)
{
    int cmp;
    const sort_struct * a = (const sort_struct *) a_;
    const sort_struct * b = (const sort_struct *) b_;
    const fmpz_mpoly_struct * apoly = a->polys + a->idx;
    const fmpz_mpoly_struct * bpoly = b->polys + b->idx;

    cmp = fmpz_cmp(&a->exp, &b->exp);
    if (cmp != 0)
        return cmp;

    return fmpz_mpoly_cmp(apoly, bpoly, a->ctx);
}

void fmpz_mpoly_factor_sort(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    sort_struct * data;
    fmpz_mpoly_struct * fc;

    if (f->num < 1)
        return;

    data = (sort_struct *) flint_malloc(f->num*sizeof(sort_struct));
    for (i = 0; i < f->num; i++)
    {
        data[i].idx = i;
        data[i].exp = f->exp[i];
        data[i].polys = f->poly;
        data[i].ctx = ctx;
    }

    qsort(data, f->num, sizeof(sort_struct), _sort);

    /* we will not permute in place */
    fc = (fmpz_mpoly_struct *) flint_malloc(f->num*sizeof(fmpz_mpoly_struct));
    memcpy(fc, f->poly, f->num*sizeof(fmpz_mpoly_struct));

    for (i = 0; i < f->num; i++)
    {
        f->exp[i] = data[i].exp;
        f->poly[i] = fc[data[i].idx];
    }
    
    flint_free(fc);
    flint_free(data);
}
