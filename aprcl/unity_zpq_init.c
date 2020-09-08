/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zpq_init(unity_zpq f, ulong q, ulong p, const fmpz_t n)
{
    slong i;

    f->p = p;
    f->q = q;

    fmpz_mod_ctx_init(f->ctx, n);
    f->polys = (fmpz_mod_poly_t *) flint_malloc(p * sizeof(fmpz_mod_poly_t));

    for (i = 0; i < p; i++)
    {
        fmpz_mod_poly_init(f->polys[i], f->ctx);
    }
}

