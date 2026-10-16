/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_mulmid(fmpz *res, const fmpz *poly1, slong len1,
                                      const fmpz *poly2, slong len2,
                                      slong nlo, slong nhi, const fmpz_mod_ctx_t ctx)
{
    _fmpz_poly_mulmid(res, poly1, len1, poly2, len2, nlo, nhi);
    _fmpz_mod_vec_set_fmpz_vec(res, res, nhi - nlo, ctx);
}

void fmpz_mod_poly_mulmid(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                const fmpz_mod_poly_t poly2, slong nlo, slong nhi, const fmpz_mod_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len;

    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi >= 0);

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;

    if ((res == poly1) || (res == poly2))
    {
        fmpz *t = _fmpz_vec_init(len);
        _fmpz_mod_poly_mulmid(t, poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi, ctx);
        _fmpz_vec_clear(res->coeffs, res->alloc);
        res->coeffs = t;
        res->alloc  = len;
        res->length = len;
        _fmpz_mod_poly_normalise(res);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, len, ctx);
        _fmpz_mod_poly_mulmid(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, nlo, nhi, ctx);
        _fmpz_mod_poly_set_length(res, len);
        _fmpz_mod_poly_normalise(res);
    }
}
