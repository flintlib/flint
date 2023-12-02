/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpq.h"
#include "fmpq_poly.h"

void
_fmpq_poly_power_sums(fmpz * res, fmpz_t rden, const fmpz * poly, slong len,
                      slong n)
{
    if (fmpz_is_one(poly + len - 1))
    {
        /* monic polynomial */
        _fmpz_poly_power_sums_naive(res, poly, len, n);
        fmpz_one(rden);
    }
    else if (len == 2)
    {
        /* degree 1 */
        fmpz_t a;
        slong i;

        fmpz_init(a);
        fmpz_set(a, poly + 1);
        fmpz_set(rden, poly);

        if (fmpz_sgn(a) < 0)
            fmpz_neg(a, a);
        else
            fmpz_neg(rden, rden);

        fmpz_one(res);

        for (i = 1; i < n; i++)
            fmpz_mul(res + i, res + i - 1, rden);
        fmpz_one(rden);

        for (i = n - 2; i >= 0; i--)
        {
            fmpz_mul(rden, rden, a);
            fmpz_mul(res + i, res + i, rden);
        }
        fmpz_set(res, rden);
        fmpz_clear(a);
    }
    else
    {
        /* general case */
        slong i, k;

        fmpz_one(rden);

        for (k = 1; k < FLINT_MIN(n, len); k++)
        {
            fmpz_mul_ui(res + k, poly + len - 1 - k, k);
            fmpz_mul(res + k, res + k, rden);

            for (i = 1; i < k - 1; i++)
                fmpz_mul(res + i, res + i, poly + len - 1);
            for (i = 1; i < k; i++)
                fmpz_addmul(res + k, poly + len - 1 - k + i, res + i);
            fmpz_neg(res + k, res + k);
            fmpz_mul(rden, rden, poly + len - 1);
        }

        for (k = len; k < n; k++)
        {
            fmpz_zero(res + k);
            for (i = k - len + 1; i < k - 1; i++)
                fmpz_mul(res + i, res + i, poly + len - 1);
            for (i = k - len + 1; i < k; i++)
                fmpz_addmul(res + k, poly + len - 1 - k + i, res + i);
            fmpz_neg(res + k, res + k);
        }

        for (i = n - len + 1; i < n - 1; i++)
            fmpz_mul(res + i, res + i, poly + len - 1);
        fmpz_one(rden);

        for (i = n - len; i > 0; i--)
        {
            fmpz_mul(rden, rden, poly + len - 1);
            fmpz_mul(res + i, res + i, rden);
        }
        fmpz_pow_ui(rden, poly + len - 1, n - 1);
        fmpz_mul_ui(res, rden, len - 1);
    }
}

void
fmpq_poly_power_sums(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpq_poly_power_sums_naive): Zero polynomial.\n");
    }
    else if (n <= 0 || poly->length == 1)
    {
        fmpq_poly_zero(res);
        return;
    }
    else if (n == 1)
    {
        fmpq_poly_fit_length(res, 1);
        fmpz_set_ui(res->coeffs, fmpq_poly_degree(poly));
        fmpz_one(res->den);
        _fmpq_poly_set_length(res, 1);
        return;
    }
    else if (poly == res)
    {
        fmpq_poly_t t;
        fmpq_poly_init(t);
        fmpq_poly_fit_length(t, n);
        _fmpq_poly_power_sums(t->coeffs, t->den, poly->coeffs,
                              poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }
    else
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_power_sums(res->coeffs, res->den, poly->coeffs,
                              poly->length, n);
    }
    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
    fmpq_poly_canonicalise(res);
}
