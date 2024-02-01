/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"

void
nmod_poly_set_trunc(nmod_poly_t res, const nmod_poly_t poly, slong len)
{
    if (poly == res)
    {
        if (res->length > len)
        {
            res->length = len;
            _nmod_poly_normalise(res);
        }
    }
    else
    {
        slong rlen;

        rlen = FLINT_MIN(len, poly->length);
        while (rlen > 0 && poly->coeffs[rlen - 1] == 0)
            rlen--;

        nmod_poly_fit_length(res, rlen);
        _nmod_vec_set(res->coeffs, poly->coeffs, rlen);
        _nmod_poly_set_length(res, rlen);
    }
}
