/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"

void
_nmod_poly_evaluate_nmod_vec(ulong_ptr ys, ulong_srcptr coeffs, slong len,
    ulong_srcptr xs, slong n, nmod_t mod)
{
    if (len < 32)
        _nmod_poly_evaluate_nmod_vec_iter(ys, coeffs, len, xs, n, mod);
    else
        _nmod_poly_evaluate_nmod_vec_fast(ys, coeffs, len, xs, n, mod);
}

void
nmod_poly_evaluate_nmod_vec(ulong_ptr ys,
        const nmod_poly_t poly, ulong_srcptr xs, slong n)
{
    _nmod_poly_evaluate_nmod_vec(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}
