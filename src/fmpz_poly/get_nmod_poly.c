/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
fmpz_poly_get_nmod_poly(nmod_poly_t res, const fmpz_poly_t poly)
{
    slong len = poly->length;

    if (len == 0)
    {
        nmod_poly_zero(res);
    }
    else
    {
        nmod_poly_fit_length(res, len);
        _fmpz_vec_get_nmod_vec(res->coeffs, poly->coeffs, len, res->mod);
        _nmod_poly_set_length(res, len);
        _nmod_poly_normalise(res);
    }
}
