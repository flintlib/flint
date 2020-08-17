/*
    Copyright (C) 2011, 2016 Fredrik Johansson

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
#include "fmpq_poly.h"

void
_fmpq_poly_sin_cos_series_basecase_can(fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden, slong Alen, slong n, int can);

void
_fmpq_poly_cos_series_basecase(fmpz * g, fmpz_t gden, 
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    fmpz * tmp;

    if (hlen == 1 || n == 1)
    {
        fmpz_one(g);
        _fmpz_vec_zero(g + 1, n - 1);
        fmpz_one(gden);
        return;
    }

    tmp = _fmpz_vec_init(n + 1);
    _fmpq_poly_sin_cos_series_basecase_can(tmp, tmp + 1, g, gden, h, hden, hlen, n, 2);
    _fmpz_vec_clear(tmp, n + 1);
}

void
_fmpq_poly_cos_series_tangent(fmpz * g, fmpz_t gden, 
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    fmpz * t;
    fmpz * u;
    fmpz_t tden;
    fmpz_t uden;

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);

    /* cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2) */
    fmpz_mul_ui(uden, hden, 2);
    _fmpq_poly_tan_series(t, tden, h, uden, hlen, n);
    _fmpq_poly_mullow(u, uden, t, tden, n, t, tden, n, n);

    _fmpq_poly_canonicalise(u, uden, n);
    _fmpz_vec_neg(t, u, n);
    fmpz_set(tden, uden);
    fmpz_set(t, tden);
    fmpz_set(u, uden);

    _fmpq_poly_div_series(g, gden, t, tden, n, u, uden, n, n);
    _fmpq_poly_canonicalise(g, gden, n);

    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    fmpz_clear(tden);
    fmpz_clear(uden);
}

void
_fmpq_poly_cos_series(fmpz * g, fmpz_t gden, 
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    if (hlen < 20 || n < 20)
        _fmpq_poly_cos_series_basecase(g, gden, h, hden, hlen, n);
    else
        _fmpq_poly_cos_series_tangent(g, gden, h, hden, hlen, n);
}

void fmpq_poly_cos_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (poly->length == 0 || n == 1)
    {
        fmpq_poly_one(res);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        flint_printf("Exception (fmpq_poly_cos_series). Constant term != 0.\n");
        flint_abort();
    }

    fmpq_poly_fit_length(res, n);
    _fmpq_poly_cos_series(res->coeffs, res->den,
        poly->coeffs, poly->den, poly->length, n);
    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
