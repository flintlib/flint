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
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"


void
_nmod_poly_interpolate_nmod_vec(mp_ptr poly,
                            mp_srcptr xs, mp_srcptr ys, slong n, nmod_t mod)
{
    if (n < 6)
        _nmod_poly_interpolate_nmod_vec_newton(poly, xs, ys, n, mod);
    else if (n < 16)
        _nmod_poly_interpolate_nmod_vec_barycentric(poly, xs, ys, n, mod);
    else
        _nmod_poly_interpolate_nmod_vec_fast(poly, xs, ys, n, mod);
}

void
nmod_poly_interpolate_nmod_vec(nmod_poly_t poly,
                                    mp_srcptr xs, mp_srcptr ys, slong n)
{
    if (n == 0)
    {
        nmod_poly_zero(poly);
    }
    else
    {
        nmod_poly_fit_length(poly, n);
        poly->length = n;
        _nmod_poly_interpolate_nmod_vec(poly->coeffs,
            xs, ys, n, poly->mod);
        _nmod_poly_normalise(poly);
    }
}
