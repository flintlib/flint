/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"

void
_nmod_poly_taylor_shift(mp_ptr poly, mp_limb_t c, slong len, nmod_t mod)
{
    if (len < 100 || len > mod.n)
        _nmod_poly_taylor_shift_horner(poly, c, len, mod);
    else if ((c == 1 || c == mod.n - 1) && len < 1000)
        _nmod_poly_taylor_shift_horner(poly, c, len, mod);
    else
        _nmod_poly_taylor_shift_convolution(poly, c, len, mod);
}

void
nmod_poly_taylor_shift(nmod_poly_t g, const nmod_poly_t f, mp_limb_t c)
{
    if (f != g)
        nmod_poly_set(g, f);

    _nmod_poly_taylor_shift(g->coeffs, c, g->length, g->mod);
}
