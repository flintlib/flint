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
unity_zpq_scalar_mul_ui(unity_zpq f, const unity_zpq g, ulong s)
{
    slong i;

    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));

    for (i = 0; i < f->p; i++)
        fmpz_mod_poly_scalar_mul_ui(f->polys[i], g->polys[i], s, f->ctx);
}

