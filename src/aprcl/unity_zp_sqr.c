/*
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_poly.h"
#include "aprcl.h"

void
unity_zp_sqr(unity_zp f, const unity_zp g)
{
    slong i;

    if (g->poly->length == 0)
    {
        fmpz_mod_poly_zero(f->poly, f->ctx);
        return;
    }

    {
        UNITY_ZP_MUL_BEGIN(t, g->poly->length, g->poly->length)
        _unity_zp_mul_reduce(f, g, g, t);
        UNITY_ZP_MUL_END(t)
    }
}

void
unity_zp_sqr_inplace(unity_zp f, const unity_zp g, fmpz_t * t)
{
    if (2 * g->poly->length - 1 <= SQUARING_SPACE)
    {
        if (!_unity_zp_sqr_special(f, g, (fmpz *) t))
            _unity_zp_mul_reduce(f, g, g, (fmpz *) t);
    }
    else
        unity_zp_sqr(f, g);
}
