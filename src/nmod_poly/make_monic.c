/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_make_monic(nn_ptr res,
                            nn_srcptr poly, slong len, nmod_t mod)
{
    ulong inv;

    inv = n_invmod(poly[len - 1], mod.n);
    _nmod_vec_scalar_mul_nmod(res, poly, len, inv, mod);
}

void nmod_poly_make_monic(nmod_poly_t res, const nmod_poly_t poly)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_DIVZERO, "Exception (nmod_poly_make_monic). Division by zero.\n");
    }

    nmod_poly_fit_length(res, poly->length);
    _nmod_poly_make_monic(res->coeffs, poly->coeffs, poly->length, poly->mod);
    res->length = poly->length;
}
