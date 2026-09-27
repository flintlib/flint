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

/* Set (rpoly, rden) to (rpoly, rden) + x^m (w, wden), where w has
   length L and (w, wden) is overwritten. The result is canonical if both
   inputs are. */
static void
_fmpq_poly_append_shifted(fmpz * rpoly, fmpz_t rden, slong m,
                          fmpz * w, fmpz_t wden, slong L)
{
    fmpz_t d, t;

    fmpz_init(d);
    fmpz_init(t);

    fmpz_gcd(d, rden, wden);

    /* scale low part by wden / d, high part by rden / d */
    fmpz_divexact(t, wden, d);
    if (!fmpz_is_one(t))
        _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, m, t);
    fmpz_divexact(d, rden, d);
    if (!fmpz_is_one(d))
        _fmpz_vec_scalar_mul_fmpz(rpoly + m, w, L, d);
    else
        _fmpz_vec_set(rpoly + m, w, L);
    fmpz_mul(rden, rden, t);

    fmpz_clear(d);
    fmpz_clear(t);
}

/*
    Third-order Newton iteration: given g = h^(-1/2) + O(x^m), let
    e = 1 - h g^2 = O(x^m). Then

        h^(-1/2) = g (1 + e/2 + 3e^2/8) + O(x^(3m)).
*/
void
_fmpq_poly_invsqrt_series(fmpz * rpoly, fmpz_t rden,
                      const fmpz * poly, const fmpz_t den, slong len, slong n)
{
    slong a[FLINT_BITS];
    slong i, m, k, L, tlen, hnlen, plen;
    fmpz * t, * u, * v;
    fmpz_t uden, wden;

    fmpz_one(rpoly);
    fmpz_one(rden);

    if (n == 1)
        return;

    len = FLINT_MIN(len, n);

    if (len == 1)
    {
        _fmpz_vec_zero(rpoly + 1, n - 1);
        return;
    }

    a[i = 0] = n;
    while (n > 1)
        a[++i] = (n = (n + 2) / 3);

    t = _fmpz_vec_init(3 * a[0]);
    u = t + a[0];
    v = u + a[0];

    fmpz_init(uden);
    fmpz_init(wden);

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];
        k = n - 2 * m;
        L = n - m;

        /* u / uden = (h g^2)[m, n) = -e / x^m */
        tlen = FLINT_MIN(2 * m - 1, n);
        hnlen = FLINT_MIN(len, n);
        plen = FLINT_MIN(n, tlen + hnlen - 1);
        _fmpz_poly_sqrlow(t, rpoly, m, tlen);
        _fmpz_poly_mulmid(u, t, tlen, poly, hnlen, m, plen);
        _fmpz_vec_zero(u + plen - m, n - plen);
        fmpz_mul(uden, rden, rden);
        fmpz_mul(uden, uden, den);

        if (k > 0)
        {
            /* e/2 + 3e^2/8 = -(4 uden u - 3 x^m u^2) / (8 uden^2) */
            _fmpz_poly_sqrlow(v, u, k, k);
            _fmpz_vec_scalar_mul_ui(v, v, k, 3);
            _fmpz_vec_scalar_mul_fmpz(u, u, L, uden);
            _fmpz_vec_scalar_mul_2exp(u, u, L, 2);
            _fmpz_vec_sub(u + m, u + m, v, k);
            fmpz_mul(wden, uden, uden);
            fmpz_mul_2exp(wden, wden, 3);
        }
        else
        {
            /* e/2 = -u / (2 uden) */
            fmpz_mul_2exp(wden, uden, 1);
        }

        _fmpz_poly_mullow(t, rpoly, FLINT_MIN(m, L), u, L, L);
        _fmpz_vec_neg(t, t, L);
        fmpz_mul(wden, wden, rden);
        _fmpq_poly_canonicalise(t, wden, L);
        _fmpq_poly_append_shifted(rpoly, rden, m, t, wden, L);
    }

    fmpz_clear(uden);
    fmpz_clear(wden);
    _fmpz_vec_clear(t, 3 * a[0]);
}

void fmpq_poly_invsqrt_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length < 1 || !fmpz_equal(poly->coeffs, poly->den))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_invsqrt_series). Constant term != 1.\n");
    }

    if (n < 1)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_invsqrt_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_invsqrt_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    fmpq_poly_canonicalise(res); /* XXX: necessary? */
}
