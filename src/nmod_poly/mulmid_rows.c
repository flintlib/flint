/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"
#include "impl.h"

/*
    Coefficients nlo to nhi - 1 of the product, one scalar product of poly1
    per coefficient of poly2 (see _nmod_poly_want_rows in impl.h). Assumes
    len1, len2 >= 1 and 0 <= nlo < nhi <= len1 + len2 - 1, does not allow
    aliasing.
*/
void
_nmod_poly_mulmid_rows(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    slong j, lo, hi;

    /* the first coefficient of poly2 whose row reaches coefficient nlo */
    j = FLINT_MAX(0, nlo - len1 + 1);
    hi = FLINT_MIN(j + len1, nhi);
    _nmod_vec_scalar_mul_nmod(res, poly1 + nlo - j, hi - nlo, poly2[j], mod);
    _nmod_vec_zero(res + hi - nlo, nhi - hi);

    for (j++; j < FLINT_MIN(len2, nhi); j++)
    {
        lo = FLINT_MAX(j, nlo);
        hi = FLINT_MIN(j + len1, nhi);
        _nmod_vec_scalar_addmul_nmod(res + lo - nlo, poly1 + lo - j, hi - lo,
                                     poly2[j], mod);
    }
}

/*
    The same when the sums of len2 products fit in one word
    (_nmod_poly_one_word): the rows accumulate without reduction, with
    _nmod_vec_nored_scalar_addmul_halflimb, and the outputs are reduced once
    at the end.
*/
void
_nmod_poly_mulmid_rows_one_word(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    slong j, lo, hi;

    _nmod_vec_zero(res, nhi - nlo);

    for (j = FLINT_MAX(0, nlo - len1 + 1); j < FLINT_MIN(len2, nhi); j++)
    {
        lo = FLINT_MAX(j, nlo);
        hi = FLINT_MIN(j + len1, nhi);
        _nmod_vec_nored_scalar_addmul_halflimb(res + lo - nlo, poly1 + lo - j,
                                               hi - lo, poly2[j]);
    }

    _nmod_vec_reduce(res, res, nhi - nlo, mod);
}
