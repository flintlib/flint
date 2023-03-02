/*
    Copyright (C) 2012, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

/* Improve cache locality with huge coefficients */
#define BLOCK_SIZE 32

void
_fmpz_poly_taylor_shift_horner(fmpz * poly, const fmpz_t c, slong n)
{
    slong i, j, bits, out_bits;

    if (n <= 1 || fmpz_is_zero(c))
        return;

    if (!fmpz_is_one(c))
    {
        if (n <= 4 || (*c != WORD(-1) && n <= 10))
        {
            if (*c == WORD(-1))
            {
                for (i = n - 2; i >= 0; i--)
                    for (j = i; j < n - 1; j++)
                        fmpz_sub(poly + j, poly + j, poly + j + 1);
            }
            else
            {
                for (i = n - 2; i >= 0; i--)
                    for (j = i; j < n - 1; j++)
                        fmpz_addmul(poly + j, poly + j + 1, c);
            }
        }
        else
        {
            slong i;
            fmpz_t d, one;
            *one = 1;

            if (*c == WORD(-1))
            {
                for (i = 1; i < n; i += 2)
                    fmpz_neg(poly + i, poly + i);

                _fmpz_poly_taylor_shift_horner(poly, one, n);

                for (i = 1; i < n; i += 2)
                    fmpz_neg(poly + i, poly + i);
            }
            else
            {
                fmpz_init(d);

                for (i = 1; i < n; i++)
                {
                    if (i == 1)
                        fmpz_set(d, c);
                    else
                        fmpz_mul(d, d, c);

                    fmpz_mul(poly + i, poly + i, d);
                }

                _fmpz_poly_taylor_shift_horner(poly, one, n);

                for (i = n - 1; i >= 1; i--)
                {
                    fmpz_divexact(poly + i, poly + i, d);

                    if (i >= 3)
                        fmpz_divexact(d, d, c);
                    else if (i == 2)
                        fmpz_set(d, c);
                }

                fmpz_clear(d);
            }
        }

        return;
    }

    bits = _fmpz_vec_max_bits(poly, n);
    bits = FLINT_ABS(bits);
    out_bits = bits + n + 1;

    if (out_bits <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                poly[j] += poly[j + 1];
    }
    else if (n >= 5 && out_bits < 2 * FLINT_BITS)
    {
        mp_ptr t;
        TMP_INIT;
        TMP_START;

        t = TMP_ALLOC(2 * n * sizeof(mp_limb_t));

        for (i = 0; i < n; i++)
            fmpz_get_signed_uiui(t + 2*i + 1, t + 2*i + 0, poly + i);

        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                add_ssaaaa(t[2*j + 1], t[2*j], t[2*j + 1], t[2*j],
                                               t[2*(j + 1) + 1], t[2*(j + 1)]);

        for (i = 0; i < n; i++)
            fmpz_set_signed_uiui(poly + i, t[2*i + 1], t[2*i]);

        TMP_END;
    }
    else if (n >= 3 + FLINT_BIT_COUNT(bits) && bits <= 100 * FLINT_BITS)
    {
        slong B = BLOCK_SIZE, ii, jj, d;
        mp_ptr t;
        TMP_INIT;
        TMP_START;

        d = (out_bits + FLINT_BITS - 1) / FLINT_BITS;
        t = TMP_ALLOC(d * n * sizeof(mp_limb_t));

        for (i = 0; i < n; i++)
            fmpz_get_signed_ui_array(t + d * i, d, poly + i);

        if (d == 3)
        {
            for (i = n - 2; i >= 0; i--)
                for (j = i; j < n - 1; j++)
                    add_sssaaaaaa(t[3 * j + 2], t[3 * j + 1], t[3 * j],
                                  t[3 * j + 2], t[3 * j + 1], t[3 * j],
                                  t[3 * (j + 1) + 2], t[3 * (j + 1) + 1], t[3 * (j + 1)]);
        }
        else if (n < 2 * B)
        {
            for (i = n - 2; i >= 0; i--)
                for (j = i; j < n - 1; j++)
                    mpn_add_n(t + d * j, t + d * j, t + d * (j + 1), d);
        }
        else
        {
            for (ii = n - 2; ii >= 0; ii -= B)
                for (jj = 0; jj < n - ii - 1; jj += B)
                    for (i = ii; i >= FLINT_MAX(ii - B + 1, 0); i--)
                        for (j = jj; j < FLINT_MIN(jj + B, n - i - 1); j++)
                            mpn_add_n(t + d * (i + j), t + d * (i + j), t + d * (i + j + 1), d);
        }

        for (i = 0; i < n; i++)
            fmpz_set_signed_ui_array(poly + i, t + d * i, d);

        TMP_END;
    }
    else
    {
        slong B = BLOCK_SIZE, ii, jj;

        if (n < 2 * B)
        {
            for (i = n - 2; i >= 0; i--)
                for (j = i; j < n - 1; j++)
                    fmpz_add(poly + j, poly + j, poly + j + 1);
        }
        else
        {
            for (ii = n - 2; ii >= 0; ii -= B)
                for (jj = 0; jj < n - ii - 1; jj += B)
                    for (i = ii; i >= FLINT_MAX(ii - B + 1, 0); i--)
                        for (j = jj; j < FLINT_MIN(jj + B, n - i - 1); j++)
                            fmpz_add(poly + i + j, poly + i + j, poly + i + j + 1);
        }
    }
}

void
fmpz_poly_taylor_shift_horner(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c)
{
    if (f != g)
        fmpz_poly_set(g, f);

    _fmpz_poly_taylor_shift_horner(g->coeffs, c, g->length);
}
