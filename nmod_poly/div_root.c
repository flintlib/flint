/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

mp_limb_t
_nmod_poly_div_root(mp_ptr Q, mp_srcptr A, slong len, mp_limb_t c, nmod_t mod)
{
    mp_limb_t r, t;
    slong i;

    if (len < 2)
        return UWORD(0);

    t = A[len - 2];
    r = Q[len - 2] = A[len - 1];

    for (i = len - 2; i > 0; i--)
    {
        r = nmod_add(nmod_mul(r, c, mod), t, mod);
        t = A[i-1];
        Q[i-1] = r;
    }

    r = nmod_add(nmod_mul(r, c, mod), t, mod);
    return r;
}

mp_limb_t
nmod_poly_div_root(nmod_poly_t Q, 
                 const nmod_poly_t A, mp_limb_t c)
{
    mp_limb_t rem;

    slong len = A->length;

    if (len == 0)
    {
        nmod_poly_zero(Q);
        return UWORD(0);
    }

    if (len == 1)
    {
        rem = A->coeffs[0];
        nmod_poly_zero(Q);
        return rem;
    }

    if (c == 0)
    {
        rem = A->coeffs[0];
        nmod_poly_shift_right(Q, A, 1);
        return rem;
    }

    nmod_poly_fit_length(Q, len - 1);
    rem = _nmod_poly_div_root(Q->coeffs, A->coeffs, len, c, Q->mod);
    Q->length = len - 1;
    return rem;
}
