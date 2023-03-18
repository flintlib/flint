/*
    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "gr_poly.h"

mp_limb_t
_nmod_poly_resultant_hgcd(mp_srcptr poly1, slong len1,
                               mp_srcptr poly2, slong len2, nmod_t mod)
{
    gr_ctx_t ctx;
    slong cutoff = NMOD_BITS(mod) <= 8 ? NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;
    mp_limb_t res;

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_resultant_hgcd(&res, poly1, len1, poly2, len2, NMOD_POLY_HGCD_CUTOFF, cutoff, ctx));

    return res;
}

mp_limb_t
nmod_poly_resultant_hgcd(const nmod_poly_t f, const nmod_poly_t g)
{
    const slong len1 = f->length;
    const slong len2 = g->length;
    mp_limb_t r;

    if (len1 == 0 || len2 == 0)
    {
        r = 0;
    }
    else
    {
        if (len1 >= len2)
        {
            r = _nmod_poly_resultant_hgcd(f->coeffs, len1, g->coeffs, len2, f->mod);
        }
        else
        {
            r = _nmod_poly_resultant_hgcd(g->coeffs, len2, f->coeffs, len1, f->mod);
            if (((len1 | len2) & WORD(1)) == WORD(0))
                r = nmod_neg(r, f->mod);
        }
    }

    return r;
}

