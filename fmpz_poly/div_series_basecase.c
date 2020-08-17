/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
_fmpz_poly_div_series_basecase(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, slong n)
{
    slong i;
    fmpz_t r;

    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (Blen == 1)
    {
        if (fmpz_is_pm1(B))
        {
            if (fmpz_is_one(B))
                _fmpz_vec_set(Q, A, Alen);
            else
                _fmpz_vec_neg(Q, A, Alen);
        } else
        {
            fmpz_init(r);

            for (i = 0; i < Alen; i++)
            {   
                fmpz_fdiv_qr(Q + i, r, A + i, B + 0);

                if (!fmpz_is_zero(r))
                {
                    fmpz_clear(r);

                    flint_printf("Not an exact division\n");
                    flint_abort();
                }
            }

            fmpz_clear(r);
        }

        _fmpz_vec_zero(Q + Alen, n - Alen);
    }
    else
    {
        slong i, j;

        if (fmpz_is_pm1(B))
        {
            if (fmpz_is_one(B))
                fmpz_set(Q, A);
            else
                fmpz_neg(Q, A);
        } else
        {
            fmpz_init(r);

            fmpz_fdiv_qr(Q + 0, r, A + 0, B + 0);

            if (!fmpz_is_zero(r))
            {
                fmpz_clear(r);

                flint_printf("Not an exact division\n");
                flint_abort();
            }
        }

        for (i = 1; i < n; i++)
        {
            fmpz_mul(Q + i, B + 1, Q + i - 1);

            for (j = 2; j < FLINT_MIN(i + 1, Blen); j++)
                fmpz_addmul(Q + i, B + j, Q + i - j);

            if (i < Alen)
            {
                if (fmpz_is_pm1(B))
                {
                    if (fmpz_is_one(B))
                        fmpz_sub(Q + i, A + i, Q + i);
                    else
                        fmpz_sub(Q + i, Q + i, A + i);
                } else
                {
                    fmpz_sub(Q + i, A + i, Q + i);

                    fmpz_fdiv_qr(Q + i, r, Q + i, B + 0);

                    if (!fmpz_is_zero(r))
                    {
                        fmpz_clear(r);

                        flint_printf("Not an exact division\n");
                        flint_abort();
                    }
                }
            } else
            {
                if (fmpz_is_pm1(B))
                {
                    if (fmpz_is_one(B))
                        fmpz_neg(Q + i, Q + i);
                } else
                {
                    fmpz_neg(Q + i, Q + i);

                    fmpz_fdiv_qr(Q + i, r, Q + i, B + 0);

                    if (!fmpz_is_zero(r))
                    {
                        fmpz_clear(r);

                        flint_printf("Not an exact division\n");
                        flint_abort();
                    }
                }
            }
        }

        if (!fmpz_is_pm1(B))
            fmpz_clear(r);
    }
}

void fmpz_poly_div_series_basecase(fmpz_poly_t Q, const fmpz_poly_t A,
                                         const fmpz_poly_t B, slong n)
{
    slong Alen = FLINT_MIN(A->length, n);
    slong Blen = FLINT_MIN(B->length, n);

    if (Blen == 0)
    {
        flint_printf("Exception (fmpz_poly_div_series_basecase). Division by zero.\n");
        flint_abort();
    }

    if (Alen == 0)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_div_series_basecase(t->coeffs, A->coeffs, Alen, B->coeffs, Blen, n);
        fmpz_poly_swap(Q, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(Q, n);
        _fmpz_poly_div_series_basecase(Q->coeffs, A->coeffs, Alen, B->coeffs, Blen, n);
    }

    _fmpz_poly_set_length(Q, n);
    _fmpz_poly_normalise(Q);
}

