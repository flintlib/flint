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
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_set_trunc(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                             slong n, const fmpz_mod_ctx_t ctx)
{
    if (poly == res)
    {
        fmpz_mod_poly_truncate(res, n, ctx);
    }
    else
    {
        slong rlen;

        rlen = FLINT_MIN(n, poly->length);
        while (rlen > 0 && fmpz_is_zero(poly->coeffs + rlen - 1))
            rlen--;

        fmpz_mod_poly_fit_length(res, rlen, ctx);
        _fmpz_vec_set(res->coeffs, poly->coeffs, rlen);
        _fmpz_mod_poly_set_length(res, rlen);
    }
}

