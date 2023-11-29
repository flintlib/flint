/*
    Copyright (C) 2010 Sebastian Pancratz
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
_fmpq_poly_tan_series(fmpz * g, fmpz_t gden,
                    const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    slong m;
    fmpz * t, * u, * v;
    fmpz_t tden, uden, vden;

    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        fmpz_zero(g);
        fmpz_one(gden);
        _fmpz_vec_zero(g + hlen, n - hlen);
        return;
    }

    if (n <= 3)
    {
        fmpz_zero(g);
        if (n >= 2)
            fmpz_set(g + 1, h + 1);
        if (hlen == 3)
            fmpz_set(g + 2, h + 2);
        else if (n == 3)
            fmpz_zero(g + 2);
        fmpz_set(gden, hden);
        _fmpq_poly_canonicalise(g, gden, n);
        return;
    }

    m = (n + 1) / 2;

    _fmpq_poly_tan_series(g, gden, h, hden, hlen, m);
    _fmpz_vec_zero(g + m, n - m);

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    v = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);
    fmpz_init(vden);

    _fmpq_poly_mul(u, uden, g, gden, m, g, gden, m);
    fmpz_set(u, uden); /* u += 1 */
    if (2*m - 1 < n)
        fmpz_zero(u + n - 1);

    _fmpq_poly_atan_series(t, tden, g, gden, n, n);
    _fmpq_poly_sub(t, tden, t, tden, n, h, hden, hlen);
    _fmpq_poly_mullow(v + m, vden, u, uden, n, t + m, tden, n - m, n - m);
    _fmpq_poly_sub(g, gden, g, gden, m, v, vden, n);
    _fmpq_poly_canonicalise(g, gden, n);

    fmpz_clear(tden);
    fmpz_clear(uden);
    fmpz_clear(vden);
    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    _fmpz_vec_clear(v, n);
}

void fmpq_poly_tan_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length && !fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_tan_series). Constant term != 0.\n");
    }

    if (poly->length == 0 || n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_tan_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_tan_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
