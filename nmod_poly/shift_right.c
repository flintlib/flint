/*
    Copyright (C) 2010 William Hart

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

void _nmod_poly_shift_right(mp_ptr res, mp_srcptr poly, slong len, slong k)
{
    flint_mpn_copyi(res, poly + k, len);
}

void nmod_poly_shift_right(nmod_poly_t res, const nmod_poly_t poly, slong k)
{
    if (k >= poly->length) /* shift all coeffs out */
        res->length = 0;
    else
    {
        const slong len = poly->length - k;
        nmod_poly_fit_length(res, len);

        _nmod_poly_shift_right(res->coeffs, poly->coeffs, len, k);

        res->length = len;
    }
}
