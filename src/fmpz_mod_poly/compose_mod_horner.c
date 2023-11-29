/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_compose_mod_horner(fmpz * res, const fmpz * f, slong lenf, const fmpz * g,
                                              const fmpz * h, slong lenh, const fmpz_mod_ctx_t ctx)
{
    slong i, len;
    fmpz * t;

    if (lenh == 1)
        return;

    if (lenf == 1)
    {
        fmpz_set(res, f);
        return;
    }

    if (lenh == 2)
    {
        _fmpz_mod_poly_evaluate_fmpz(res, f, lenf, g, ctx);
        return;
    }

    len = lenh - 1;
    i = lenf - 1;
    t = _fmpz_vec_init(2 * lenh - 3);

    _fmpz_mod_poly_scalar_mul_fmpz(res, g, len, f + i, ctx);
    i--;
    if (i >= 0)
        fmpz_mod_add(res, res, f + i, ctx);

    while (i > 0)
    {
        i--;
        _fmpz_mod_poly_mulmod(t, res, len, g, len, h, lenh, ctx);
        _fmpz_mod_poly_add(res, t, len, f + i, 1, ctx);
    }

    _fmpz_vec_clear(t, 2 * lenh - 3);
}

void fmpz_mod_poly_compose_mod_horner(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                         const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx)
{
    fmpz_t inv3;
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;
    slong vec_len = FLINT_MAX(len3 - 1, len2);

    fmpz * ptr2;

    if (len3 == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_compose_mod_horner). Division by zero \n");
    }

    if (len1 == 0 || len3 == 1)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    if (len1 == 1)
    {
        fmpz_mod_poly_set(res, poly1, ctx);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        fmpz_mod_poly_t tmp;
        fmpz_mod_poly_init(tmp, ctx);
        fmpz_mod_poly_compose_mod_horner(tmp, poly1, poly2, poly3, ctx);
        fmpz_mod_poly_swap(tmp, res, ctx);
        fmpz_mod_poly_clear(tmp, ctx);
        return;
    }

    ptr2 = _fmpz_vec_init(vec_len);

    if (len2 <= len3 - 1)
    {
        _fmpz_vec_set(ptr2, poly2->coeffs, len2);
        _fmpz_vec_zero(ptr2 + len2, vec_len - len2);
    }
    else
    {
        fmpz_init(inv3);
        fmpz_invmod(inv3, poly3->coeffs + len, fmpz_mod_ctx_modulus(ctx));
        _fmpz_mod_poly_rem(ptr2, poly2->coeffs, len2,
                         poly3->coeffs, len3, inv3, ctx);
        fmpz_clear(inv3);
    }

    fmpz_mod_poly_fit_length(res, len, ctx);
    _fmpz_mod_poly_compose_mod_horner(res->coeffs, poly1->coeffs, len1,
                         ptr2, poly3->coeffs, len3, ctx);
    _fmpz_mod_poly_set_length(res, len);
    _fmpz_mod_poly_normalise(res);

    _fmpz_vec_clear(ptr2, vec_len);
}
