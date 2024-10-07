/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "aprcl.h"

void
unity_zpq_clear(unity_zpq f)
{
    ulong i;

    for (i = 0; i < f->p; i++)
    {
        fmpz_mod_poly_clear(f->polys[i], f->ctx);
    }

    f->p = 0;
    f->q = 0;

    fmpz_mod_ctx_clear(f->ctx);
    flint_free(f->polys);
}
