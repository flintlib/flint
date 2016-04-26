/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011, 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

static __inline__ ulong
fmpz_gcd_ui(const fmpz_t x, ulong c)
{
    return n_gcd(c, fmpz_fdiv_ui(x, c));
}

void
_fmpq_poly_exp_series_basecase(fmpz * B, fmpz_t Bden,
    const fmpz * A, const fmpz_t Aden, slong Alen, slong n)
{
    fmpz_t t, u;
    slong j, k;

    fmpz_init(t);
    fmpz_init(u);

    fmpz_fac_ui(t, n - 1);
    fmpz_pow_ui(u, Aden, n - 1);

    fmpz_mul(Bden, t, u);
    fmpz_set(B, Bden);

    for (k = 1; k < n; k++)
    {
        fmpz_zero(t);

        for (j = 1; j < FLINT_MIN(Alen, k + 1); j++)
        {
            fmpz_mul_ui(u, A + j, j);
            fmpz_addmul(t, u, B + k - j);
        }

        fmpz_mul_ui(u, Aden, k);
        fmpz_divexact(B + k, t, u);
    }

    _fmpq_poly_canonicalise(B, Bden, n);

    fmpz_clear(t);
    fmpz_clear(u);
}

void
_fmpq_poly_exp_series_newton(fmpz * g, fmpz_t gden,
                    const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    slong m;
    fmpz * t, * u;
    fmpz_t tden, uden;

    hlen = FLINT_MIN(hlen, n);

    if (hlen < 10)
    {
        _fmpq_poly_exp_series_basecase(g, gden, h, hden, hlen, n);
        return;
    }

    m = (n + 1) / 2;
    _fmpq_poly_exp_series(g, gden, h, hden, hlen, m);
    _fmpz_vec_zero(g + m, n - m);

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);

    _fmpq_poly_log_series(t, tden, g, gden, m, n);
    _fmpq_poly_sub(t, tden, t, tden, n, h, hden, hlen);
    /* TODO: half of product is redundant! */
    _fmpq_poly_mullow(u, uden, t, tden, n, g, gden, m, n);
    _fmpq_poly_sub(g, gden, g, gden, m, u, uden, n);
    _fmpq_poly_canonicalise(g, gden, n);

    fmpz_clear(tden);
    fmpz_clear(uden);
    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
}

void
_fmpq_poly_exp_series(fmpz * B, fmpz_t Bden,
    const fmpz * A, const fmpz_t Aden, slong Alen, slong n)
{
    Alen = FLINT_MIN(Alen, n);

    if (Alen == 1)
    {
        fmpz_one(B);
        fmpz_one(Bden);
        _fmpz_vec_zero(B + 1, n - 1);
        return;
    }

    /* A is a monomial (p/q) * x^d */
    if (_fmpz_vec_is_zero(A + 1, Alen - 2))
    {
        fmpz * R;
        ulong v;
        slong i, d, m;

        d = Alen - 1;       /* Degree of input monomial. */
        m = (n - 1) / d;    /* m*d is highest degree in output. */
        R = _fmpz_vec_init(m + 1);

        /* A[d]/Aden could be non-canonical due to truncation */
        fmpz_gcd(R, A + d, Aden);
        fmpz_divexact(B + d, A + d, R);
        fmpz_divexact(R, Aden, R);  /* store q in R[0] */

        fmpz_set(R + 1, R);
        fmpz_set(Bden, R);

        for (i = 2; i <= m; i++)
        {
            /* Computing (p/q)^i / i! from the previous term, we only
               need to remove the gcd between the numerator and i. */
            fmpz_mul(B + i * d, B + (i - 1) * d, B + d);
            fmpz_mul(Bden, Bden, R);

            v = fmpz_gcd_ui(B + i * d, i);
            fmpz_divexact_ui(B + i * d, B + i * d, v);
            fmpz_mul_ui(Bden, Bden, i / v);
            fmpz_mul_ui(R + i, R, i / v);
        }

        /* Put all terms on the same denominator as the last term. */
        for (i = m - 1; i > 0; i--)
        {
            fmpz_mul(B + i * d, B + i * d, R + m);
            fmpz_mul(R + m, R + m, R + i);
        }

        /* Constant term = 1. */
        fmpz_set(B, Bden);

        for (i = 0; d != 1 && i < n; i++)
            if (i % d != 0)
                fmpz_zero(B + i);

        _fmpz_vec_clear(R, m + 1);
        return;
    }

    if (Alen < 15)
    {
        _fmpq_poly_exp_series_basecase(B, Bden, A, Aden, Alen, n);
    }
    else
    {
        _fmpq_poly_exp_series_newton(B, Bden, A, Aden, Alen, n);
    }
}

void fmpq_poly_exp_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
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
        flint_printf("Exception (fmpq_poly_exp_series). Constant term != 0.\n");
        flint_abort();
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_exp_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_exp_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}

