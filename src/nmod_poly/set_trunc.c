/*
    Copyright (C) 2014 Fredrik Johansson

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
nmod_poly_set_trunc(nmod_poly_t res, const nmod_poly_t poly, slong n)
{
    if (poly == res)
    {
        nmod_poly_truncate(res, n);
    }
    else
    {
        slong rlen;

        rlen = FLINT_MIN(n, poly->length);
        while (rlen > 0 && poly->coeffs[rlen - 1] == 0)
            rlen--;

        nmod_poly_fit_length(res, rlen);
        _nmod_vec_set(res->coeffs, poly->coeffs, rlen);
        _nmod_poly_set_length(res, rlen);
    }
}

