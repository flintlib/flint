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

void _nmod_poly_shift_left(mp_ptr res, mp_srcptr poly, slong len, slong k)
{
    flint_mpn_copyd(res + k, poly, len);
    flint_mpn_zero(res, k);
}

void nmod_poly_shift_left(nmod_poly_t res, const nmod_poly_t poly, slong k)
{
    nmod_poly_fit_length(res, poly->length + k);
   
    _nmod_poly_shift_left(res->coeffs, poly->coeffs, poly->length, k);
   
    res->length = poly->length + k;
}

