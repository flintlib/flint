/*
    Copyright (C) 2014, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_inv_series_basecase(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
{
    int neg;
    Qlen = FLINT_MIN(Qlen, n);

    neg = fmpz_is_one(Q + 0);
    fmpz_set(Qinv + 0, Q + 0);

    if (Qlen == 1)
    {
        _fmpz_vec_zero(Qinv + 1, n - 1);
    }
    else if (Qlen == 2 || _fmpz_vec_is_zero(Q + 1, Qlen - 2))
    {
        /* Special-case binomials */
        slong i, j, step;

        step = Qlen - 1;

        if (neg)
        {
            fmpz_neg(Qinv + step, Q + step);
            for (i = 2 * step; i < n; i += step)
                fmpz_mul(Qinv + i, Qinv + i - step, Qinv + step);
        }
        else
        {
            fmpz_neg(Qinv + step, Q + step);
            for (i = 2 * step; i < n; i += step)
                fmpz_mul(Qinv + i, Qinv + i - step, Q + step);
        }

        for (i = 0; i < n; i += step)
            for (j = i + 1; j < FLINT_MIN(n, i + step); j++)
                fmpz_zero(Qinv + j);
    }
    else
    {
        slong i, j, nsmall;
        char * Qbits;
        slong b, bits, Qinvbits;
        TMP_INIT;

        TMP_START;

        /* Qbits[i] = max(bits(Q[0]), ..., bits(Q[i])), as long as coeffs are small */
        Qbits = TMP_ALLOC(Qlen);
        Qbits[0] = 1;

        /* Maximum bits of all Qinv coefficients encountered so far */
        Qinvbits = 1;

        /* We have small coefficients for i < nsmall */
        for (nsmall = 1; nsmall < Qlen; nsmall++)
        {
            b = Q[nsmall];

            if (COEFF_IS_MPZ(b))
                break;

            b = FLINT_ABS(b);
            if ((b >> Qbits[nsmall - 1]) != 0)
                Qbits[nsmall] = FLINT_BIT_COUNT(b);
            else
                Qbits[nsmall] = Qbits[nsmall - 1];
        }

        for (i = 1; i < n; i++)
        {
            if (i >= nsmall || Qinvbits > SMALL_FMPZ_BITCOUNT_MAX || Qbits[i] > SMALL_FMPZ_BITCOUNT_MAX)
            {
                /* Can't use fast code. */
                bits = WORD_MAX;
            }
            else
            {
                /* Can maybe use fast code; bound bits. */
                b = FLINT_MIN(i, Qlen - 1);
                bits = FLINT_BIT_COUNT(b);

                /* Bit size of product. */
                bits += Qbits[i] + Qinvbits;

                /* Sign. */
                bits += 1;
            }

            if (bits <= 3 * FLINT_BITS - 1)
            {
                if (bits <= FLINT_BITS - 1)
                {
                    slong s, x, y;

                    s = 0;

                    for (j = 1; j < FLINT_MIN(i + 1, Qlen); j++)
                    {
                        x = Q[j];
                        y = Qinv[i - j];
                        s += x * y;
                    }

                    if (neg)
                        s = -s;

                    fmpz_set_si(Qinv + i, s);
                }
                else if (bits <= 2 * FLINT_BITS - 1)
                {
                    mp_limb_t hi, lo, shi, slo;
                    slong x, y;

                    shi = slo = 0;

                    for (j = 1; j < FLINT_MIN(i + 1, Qlen); j++)
                    {
                        x = Q[j];
                        y = Qinv[i - j];

                        smul_ppmm(hi, lo, x, y);
                        add_ssaaaa(shi, slo, shi, slo, hi, lo);
                    }

                    if (neg)
                        sub_ddmmss(shi, slo, 0, 0, shi, slo);

                    fmpz_set_signed_uiui(Qinv + i, shi, slo);
                }
                else
                {
                    mp_limb_t hi, lo, cy, shh, shi, slo;
                    slong x, y;

                    shh = shi = slo = 0;

                    for (j = 1; j < FLINT_MIN(i + 1, Qlen); j++)
                    {
                        x = Q[j];
                        y = Qinv[i - j];

                        smul_ppmm(hi, lo, x, y);
                        add_sssaaaaaa(cy, shi, slo, 0, shi, slo, 0, hi, lo);
                        shh += (0 <= (slong) hi) ? cy : cy - 1;
                    }

                    if (neg)
                        sub_dddmmmsss(shh, shi, slo, 0, 0, 0, shh, shi, slo);

                    fmpz_set_signed_uiuiui(Qinv + i, shh, shi, slo);
                }

                if (COEFF_IS_MPZ(*(Qinv + i)))
                {
                    /* Will no longer use fast code */
                    nsmall = i;
                }
                else
                {
                    /* Update Qinvbits */
                    b = FLINT_ABS(*(Qinv + i));
                    b = FLINT_BIT_COUNT(b);
                    Qinvbits = FLINT_MAX(Qinvbits, b);
                }
            }
            else
            {
                fmpz_mul(Qinv + i, Q + 1, Qinv + i - 1);

                for (j = 2; j < FLINT_MIN(i + 1, Qlen); j++)
                    fmpz_addmul(Qinv + i, Q + j, Qinv + i - j);

                if (neg)
                    fmpz_neg(Qinv + i, Qinv + i);
            }
        }

        TMP_END;
    }
}

void
fmpz_poly_inv_series_basecase(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_inv_series_basecase). Division by zero.\n");
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_inv_series_basecase(Qinv->coeffs, Q->coeffs, Qlen, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_inv_series_basecase(t->coeffs, Q->coeffs, Qlen, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }

    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
}

