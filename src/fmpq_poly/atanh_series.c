/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void
_fmpq_poly_atanh_series(fmpz * g, fmpz_t gden,
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    fmpz * t;
    fmpz * u;
    fmpz_t tden;
    fmpz_t uden;
    slong h2len;

    hlen = FLINT_MIN(hlen, n);
    h2len = FLINT_MIN(2 * hlen - 1, n);

    if (hlen == 1)
    {
        _fmpz_vec_zero(g, n);
        fmpz_one(gden);
        return;
    }

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);

    /* atanh(h(x)) = integral(h'(x)/(1-h(x)^2)) */
    _fmpq_poly_mullow(u, uden, h, hden, hlen, h, hden, hlen, h2len);
    _fmpq_poly_canonicalise(u, uden, h2len);
    _fmpz_vec_neg(u, u, h2len);
    fmpz_set(u, uden);  /* u += 1 */
    _fmpq_poly_derivative(t, tden, h, hden, hlen);
    _fmpq_poly_div_series(g, gden, t, tden, hlen - 1, u, uden, h2len, n);
    _fmpq_poly_canonicalise(g, gden, n - 1);
    _fmpq_poly_integral(g, gden, g, gden, n);

    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    fmpz_clear(tden);
    fmpz_clear(uden);
}

void fmpq_poly_atanh_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length && !fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_atanh_series). Constant term != 0.\n");
    }

    if (poly->length == 0 || n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_atanh_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_atanh_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
