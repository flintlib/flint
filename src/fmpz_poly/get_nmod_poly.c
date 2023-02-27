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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

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
        slong i;
        nmod_poly_fit_length(res, len);
        for (i = 0; i < len; i++)
            res->coeffs[i] = fmpz_fdiv_ui(poly->coeffs + i, res->mod.n);
        _nmod_poly_set_length(res, len);
        _nmod_poly_normalise(res);
    }
}
