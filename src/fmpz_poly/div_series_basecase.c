/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014, 2021 Fredrik Johansson
    Copyright (C) 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

static void
fmpz_divexact_checked(fmpz_t Q, const fmpz_t A, const fmpz_t B)
{
    fmpz_t r;
    fmpz_init(r);
    fmpz_fdiv_qr(Q, r, A, B);
    if (!fmpz_is_zero(r))
    {
        fmpz_clear(r); /* flint_throw */
        flint_throw(FLINT_ERROR, "Not an exact division\n");
    }
    /* no need to clear r */
}

void
_fmpz_poly_div_series_basecase(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, slong n)
{
    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (Blen == 1)
    {
        if (fmpz_is_pm1(B + 0))
        {
            if (fmpz_is_one(B + 0))
                _fmpz_vec_set(Q, A, Alen);
            else
                _fmpz_vec_neg(Q, A, Alen);
        }
        else
        {
            slong i;
            for (i = 0; i < Alen; i++)
                fmpz_divexact_checked(Q + i, A + i, B);
        }

        _fmpz_vec_zero(Q + Alen, n - Alen);
    }
    else if (Alen == 1 && fmpz_is_pm1(B + 0))
    {
        _fmpz_poly_inv_series_basecase(Q, B, Blen, n);
        if (!fmpz_is_one(A + 0))
            _fmpz_vec_scalar_mul_fmpz(Q, Q, n, A + 0);
    }
    else
    {
        slong i, j, nsmall;
        char *Bbits;
        slong b, bits, Qbits;
        TMP_INIT;

        TMP_START;

        if (fmpz_is_pm1(B + 0))
        {
            if (fmpz_is_one(B + 0))
                fmpz_set(Q + 0, A + 0);
            else
                fmpz_neg(Q + 0, A + 0);
        }
        else
        {
            fmpz_divexact_checked(Q + 0, A + 0, B + 0);
        }

        /* Bbits[i] = max(bits(B[0]), ..., bits(B[i])), as long as coeffs are small */
        Bbits = TMP_ALLOC(Blen);
        Bbits[0] = fmpz_bits(B + 0);

        /* Maximum bits of all Q coefficients encountered so far */
        Qbits = fmpz_bits(Q + 0);

        /* We have small coefficients for i < nsmall */
        for (nsmall = 0; nsmall < Blen; nsmall++)
        {
            b = B[nsmall];

            if (COEFF_IS_MPZ(b))
                break;

            b = FLINT_ABS(b);
            if (nsmall == 0 || (b >> Bbits[nsmall - 1]) != 0)
                Bbits[nsmall] = FLINT_BIT_COUNT(b);
            else
                Bbits[nsmall] = Bbits[nsmall - 1];
        }

        for (i = 1; i < n; i++)
        {
            if (i >= nsmall || Qbits > SMALL_FMPZ_BITCOUNT_MAX || Bbits[i] > SMALL_FMPZ_BITCOUNT_MAX)
            {
                /* Can't use fast code. */
                bits = WORD_MAX;
            }
            else
            {
                /* Can maybe use fast code; bound bits. */
                b = FLINT_MIN(i, Blen - 1);
                bits = FLINT_BIT_COUNT(b);

                /* Bit size of product. */
                bits += Bbits[i] + Qbits;

                /* Sign. */
                bits += 1;
            }

            if (bits <= 3 * FLINT_BITS - 1)
            {
                if (bits <= FLINT_BITS - 1)
                {
                    slong s, x, y;

                    s = 0;

                    for (j = 1; j < FLINT_MIN(i + 1, Blen); j++)
                    {
                        x = B[j];
                        y = Q[i - j];
                        s += x * y;
                    }

                    fmpz_set_si(Q + i, s);
                }
                else if (bits <= 2 * FLINT_BITS - 1)
                {
                    mp_limb_t hi, lo, shi, slo;
                    slong x, y;

                    shi = slo = 0;

                    for (j = 1; j < FLINT_MIN(i + 1, Blen); j++)
                    {
                        x = B[j];
                        y = Q[i - j];

                        smul_ppmm(hi, lo, x, y);
                        add_ssaaaa(shi, slo, shi, slo, hi, lo);
                    }

                    fmpz_set_signed_uiui(Q + i, shi, slo);
                }
                else
                {
                    mp_limb_t hi, lo, cy, shh, shi, slo;
                    slong x, y;

                    shh = shi = slo = 0;

                    for (j = 1; j < FLINT_MIN(i + 1, Blen); j++)
                    {
                        x = B[j];
                        y = Q[i - j];

                        smul_ppmm(hi, lo, x, y);
                        add_sssaaaaaa(cy, shi, slo, 0, shi, slo, 0, hi, lo);
                        shh += (0 <= (slong) hi) ? cy : cy - 1;
                    }

                    fmpz_set_signed_uiuiui(Q + i, shh, shi, slo);
                }
            }
            else
            {
                fmpz_mul(Q + i, B + 1, Q + i - 1);

                for (j = 2; j < FLINT_MIN(i + 1, Blen); j++)
                    fmpz_addmul(Q + i, B + j, Q + i - j);
            }

            if (i < Alen)
            {
                if (fmpz_is_pm1(B + 0))
                {
                    if (fmpz_is_one(B + 0))
                        fmpz_sub(Q + i, A + i, Q + i);
                    else
                        fmpz_sub(Q + i, Q + i, A + i);
                }
                else
                {
                    fmpz_sub(Q + i, A + i, Q + i);
                    fmpz_divexact_checked(Q + i, Q + i, B + 0);
                }
            }
            else
            {
                if (fmpz_is_pm1(B + 0))
                {
                    if (fmpz_is_one(B + 0))
                        fmpz_neg(Q + i, Q + i);
                }
                else
                {
                    fmpz_neg(Q + i, Q + i);
                    fmpz_divexact_checked(Q + i, Q + i, B + 0);
                }
            }

            if (COEFF_IS_MPZ(*(Q + i)))
            {
                /* Will no longer use fast code */
                nsmall = i;
            }
            else
            {
                /* Update Qbits */
                b = FLINT_ABS(*(Q + i));
                b = FLINT_BIT_COUNT(b);
                Qbits = FLINT_MAX(Qbits, b);
            }
        }

        TMP_END;
    }
}

void fmpz_poly_div_series_basecase(fmpz_poly_t Q, const fmpz_poly_t A,
                                         const fmpz_poly_t B, slong n)
{
    slong Alen = FLINT_MIN(A->length, n);
    slong Blen = FLINT_MIN(B->length, n);

    if (Blen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_div_series_basecase). Division by zero.\n");
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

