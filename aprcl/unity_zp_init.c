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
unity_zp_init(unity_zp f, ulong p, ulong exp, const fmpz_t n)
{
    f->p = p;
    f->exp = exp;
    fmpz_mod_ctx_init(f->ctx, n);
    fmpz_mod_poly_init(f->poly, f->ctx);
}

void
unity_zp_clear(unity_zp f)
{
    fmpz_mod_poly_clear(f->poly, f->ctx);
    fmpz_mod_ctx_clear(f->ctx);
}

