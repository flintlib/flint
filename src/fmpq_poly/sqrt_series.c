/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"


/*
    Karp-Markstein: with g = h^(-1/2) + O(x^m) and v = h g mod x^m,
    sqrt(h) = v + g (h - v^2) / 2 + O(x^(2m)).
*/
void
_fmpq_poly_sqrt_series(fmpz * rpoly, fmpz_t rden,
                      const fmpz * poly, const fmpz_t den, slong len, slong n)
{
    fmpz * g, * t, * u;
    fmpz_t gden, vden, uden, d;
    slong m, L, tlen, r1, r2, rlen;

    len = FLINT_MIN(len, n);

    if (n <= 2 || len == 1)
    {
        g = _fmpz_vec_init(n);
        fmpz_init(gden);
        _fmpq_poly_invsqrt_series(g, gden, poly, den, len, n);
        _fmpq_poly_mullow(rpoly, rden, g, gden, n, poly, den, len, n);
        _fmpq_poly_canonicalise(rpoly, rden, n);
        _fmpz_vec_clear(g, n);
        fmpz_clear(gden);
        return;
    }

    m = (n + 1) / 2;
    L = n - m;

    g = _fmpz_vec_init(m + 2 * n);
    t = g + m;
    u = t + n;

    fmpz_init(gden);
    fmpz_init(vden);
    fmpz_init(uden);
    fmpz_init(d);

    _fmpq_poly_invsqrt_series(g, gden, poly, den, len, m);

    /* v = h g mod x^m, stored in rpoly */
    _fmpz_poly_mullow(rpoly, g, m, poly, FLINT_MIN(len, m), m);
    fmpz_mul(vden, gden, den);
    _fmpq_poly_canonicalise(rpoly, vden, m);

    /* u / uden = (h - v^2)[m, n) */
    tlen = FLINT_MIN(2 * m - 1, n);
    r1 = FLINT_MAX(0, FLINT_MIN(len - m, L));
    r2 = FLINT_MAX(0, tlen - m);
    rlen = FLINT_MAX(r1, r2);

    /* common denominator of h and v^2 */
    fmpz_mul(uden, vden, vden);
    fmpz_lcm(uden, uden, den);

    _fmpz_vec_zero(u, rlen);
    if (r1 > 0)
    {
        fmpz_divexact(d, uden, den);
        _fmpz_vec_scalar_mul_fmpz(u, poly + m, r1, d);
    }
    if (r2 > 0)
    {
        _fmpz_poly_mulmid(t, rpoly, m, rpoly, m, m, tlen);
        fmpz_mul(d, vden, vden);
        fmpz_divexact(d, uden, d);
        if (!fmpz_is_one(d))
            _fmpz_vec_scalar_mul_fmpz(t, t, r2, d);
        _fmpz_vec_sub(u, u, t, r2);
    }

    /* high part = g u / (2 gden uden) */
    _fmpz_poly_mullow(t, g, L, u, rlen, L);
    fmpz_mul(d, gden, uden);
    fmpz_mul_2exp(d, d, 1);
    _fmpq_poly_canonicalise(t, d, L);

    /* rpoly / rden = v + x^m t / d, canonical */
    {
        fmpz_t e, f;
        fmpz_init(e);
        fmpz_init(f);
        fmpz_gcd(e, vden, d);
        fmpz_divexact(f, d, e);
        if (!fmpz_is_one(f))
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, m, f);
        fmpz_divexact(e, vden, e);
        _fmpz_vec_scalar_mul_fmpz(rpoly + m, t, L, e);
        fmpz_mul(rden, vden, f);
        fmpz_clear(e);
        fmpz_clear(f);
    }

    fmpz_clear(gden);
    fmpz_clear(vden);
    fmpz_clear(uden);
    fmpz_clear(d);
    _fmpz_vec_clear(g, m + 2 * n);
}

void fmpq_poly_sqrt_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length < 1 || !fmpz_equal(poly->coeffs, poly->den))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_sqrt_series). Constant term != 1.\n");
    }

    if (n < 1)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_sqrt_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_sqrt_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    fmpq_poly_canonicalise(res);
}
