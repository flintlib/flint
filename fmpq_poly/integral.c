/*
    Copyright (C) 2011, 2021 Fredrik Johansson

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

static ulong _fmpz_gcd_big_small(const fmpz_t g, ulong h)
{
    __mpz_struct * z = COEFF_TO_PTR(*g);

    return n_gcd(mpn_mod_1(z->_mp_d, FLINT_ABS(z->_mp_size), h), h);
}

static ulong _fmpz_gcd_small(const fmpz_t g, ulong h)
{
    if (!COEFF_IS_MPZ(*g))
        return n_gcd(FLINT_ABS(*g), h);
    else
        return _fmpz_gcd_big_small(g, h);
}

void _fmpq_poly_integral(fmpz * rpoly, fmpz_t rden, 
                           const fmpz * poly, const fmpz_t den, slong len)
{
    slong k;
    ulong v, c, d;
    mp_ptr divisors;
    fmpz_t t, u;
    TMP_INIT;

    if (len <= 2)
    {
        if (len == 2)
            fmpz_set(rpoly + 1, poly);
        fmpz_zero(rpoly);
        fmpz_set(rden, den);
        return;
    }

    TMP_START;
    divisors = TMP_ALLOC(sizeof(ulong) * len);

    fmpz_init(t);
    fmpz_one(t);

    for (k = len - 1; k >= 2; k--)
    {
        if (fmpz_is_zero(poly + k - 1))
        {
            fmpz_zero(rpoly + k);
        }
        else
        {
            c = _fmpz_gcd_small(poly + k - 1, k);

            if (c == k)
            {
                fmpz_divexact_ui(rpoly + k, poly + k - 1, k);
                divisors[k] = 1;
            }
            else
            {
                if (c == 1)
                {
                    fmpz_set(rpoly + k, poly + k - 1);
                    divisors[k] = k;
                }
                else
                {
                    fmpz_divexact_ui(rpoly + k, poly + k - 1, c);
                    divisors[k] = k / c;
                }

                c = divisors[k];
                d = _fmpz_gcd_small(t, c);
                if (d != c)
                    fmpz_mul_ui(t, t, c / d);
            }
        }
    }

    fmpz_mul(rden, den, t);

    if (!fmpz_is_one(t))
    {
        if (!COEFF_IS_MPZ(*t))
        {
            v = *t;
            for (k = len - 1; k >= 2; k--)
            {
                if (!fmpz_is_zero(rpoly + k) && v != divisors[k])
                    fmpz_mul_ui(rpoly + k, rpoly + k, divisors[k] == 1 ? v : v / divisors[k]);
            }
        }
        else
        {
            fmpz_init(u);

            for (k = len - 1; k >= 2; k--)
            {
                if (!fmpz_is_zero(rpoly + k))
                {
                    if (divisors[k] == 1)
                    {
                        fmpz_mul(rpoly + k, rpoly + k, t);
                    }
                    else
                    {
                        fmpz_divexact_ui(u, t, divisors[k]);
                        fmpz_mul(rpoly + k, rpoly + k, u);
                    }
                }
            }

            fmpz_clear(u);
        }
    }

    fmpz_mul(rpoly + 1, poly + 0, t);
    fmpz_zero(rpoly);

    fmpz_clear(t);
    TMP_END;
}

void fmpq_poly_integral(fmpq_poly_t res, const fmpq_poly_t poly)
{
    slong len = poly->length;

    if (len == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    fmpq_poly_fit_length(res, len + 1);
    _fmpq_poly_integral(res->coeffs, res->den, poly->coeffs, poly->den, len + 1);
    _fmpq_poly_set_length(res, len + 1);
}
