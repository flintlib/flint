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

static void
_set_vec(fmpz * rnum, fmpz_t den,
                const fmpz * xnum, const fmpz * xden, slong len)
{
    slong j;
    fmpz_t t;
    fmpz_init(t);
    fmpz_one(den);

    for (j = 0; j < len; j++)
        fmpz_lcm(den, den, xden + j);

    for (j = 0; j < len; j++)
    {
        fmpz_divexact(t, den, xden + j);
        fmpz_mul(rnum + j, xnum + j, t);
    }

    fmpz_clear(t);
}

void
_fmpq_poly_revert_series_lagrange(fmpz * Qinv, fmpz_t den,
                        const fmpz * Q, const fmpz_t Qden, slong Qlen, slong n)
{
    slong i;
    fmpz *R, *S, *T, *dens, *tmp;
    fmpz_t Rden, Sden, Tden;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen <= 2)
    {
        fmpz_zero(Qinv);

        if (Qlen == 2)
        {
            fmpz_set(Qinv + 1, Qden);
            fmpz_set(den, Q + 1);
            _fmpq_poly_canonicalise(Qinv, den, 2);
        }

        _fmpz_vec_zero(Qinv + 2, n - 2);
    }
    else
    {
        dens = _fmpz_vec_init(n);
        R = _fmpz_vec_init(n - 1);
        S = _fmpz_vec_init(n - 1);
        T = _fmpz_vec_init(n - 1);
        fmpz_init(Rden);
        fmpz_init(Sden);
        fmpz_init(Tden);

        fmpz_zero(Qinv);
        fmpz_one(dens);
        fmpz_set(Qinv + 1, Qden);
        fmpz_set(dens + 1, Q + 1);

        _fmpq_poly_inv_series(R, Rden, Q + 1, Qden, Qlen - 1, n - 1);
        _fmpq_poly_canonicalise(R, Rden, n - 1);

        _fmpz_vec_set(S, R, n - 1);
        fmpz_set(Sden, Rden);

        for (i = 2; i < n; i++)
        {
            _fmpq_poly_mullow(T, Tden, S, Sden, n - 1, R, Rden, n - 1, n - 1);
            _fmpq_poly_canonicalise(T, Tden, n - 1);
            fmpz_set(Qinv + i, T + i - 1);
            fmpz_mul_ui(dens + i, Tden, i);
            tmp = S; S = T; T = tmp;
            fmpz_swap(Sden, Tden);
        }

        _set_vec(Qinv, den, Qinv, dens, n);
        _fmpq_poly_canonicalise(Qinv, den, n);

        _fmpz_vec_clear(R, n - 1);
        _fmpz_vec_clear(S, n - 1);
        _fmpz_vec_clear(T, n - 1);
        _fmpz_vec_clear(dens, n);
        fmpz_clear(Rden);
        fmpz_clear(Sden);
        fmpz_clear(Tden);
    }
}

void
fmpq_poly_revert_series_lagrange(fmpq_poly_t res,
            const fmpq_poly_t poly, slong n)
{
    if (poly->length < 2 || !fmpz_is_zero(poly->coeffs)
                         || fmpz_is_zero(poly->coeffs + 1))
    {
        flint_throw(FLINT_ERROR, "(fmpq_poly_revert_series_lagrange): "
                "Input must have zero constant term and nonzero coefficient of x^1.\n");
    }

    if (n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_revert_series_lagrange(res->coeffs,
                res->den, poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_revert_series_lagrange(t->coeffs,
                t->den, poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
