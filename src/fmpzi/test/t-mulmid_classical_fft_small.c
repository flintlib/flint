/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpzi.h"

TEST_FUNCTION_START(fmpzi_poly_mulmid_classical_fft_small, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        slong len1 = 1 + n_randint(state, 10);
        slong len2 = 1 + n_randint(state, 10);
        slong bits = 9000 + n_randint(state, 9000);
        slong nlo, nhi, i, j;
        int sq = (n_randint(state, 4) == 0);
        int ok;
        fmpzi_struct * a, * b, * r1, * r2;
        fmpzi_t t, u;

        if (sq)
            len2 = len1;

        nlo = n_randint(state, len1 + len2 - 1);
        nhi = nlo + 1 + n_randint(state, len1 + len2 - 1 - nlo);

        a = flint_malloc(len1 * sizeof(fmpzi_struct));
        b = sq ? a : flint_malloc(len2 * sizeof(fmpzi_struct));
        r1 = flint_malloc((nhi - nlo) * sizeof(fmpzi_struct));
        r2 = flint_malloc((nhi - nlo) * sizeof(fmpzi_struct));
        for (i = 0; i < len1; i++)
            fmpzi_init(a + i);
        if (!sq)
            for (i = 0; i < len2; i++)
                fmpzi_init(b + i);
        for (i = 0; i < nhi - nlo; i++)
        {
            fmpzi_init(r1 + i);
            fmpzi_init(r2 + i);
        }
        fmpzi_init(t);
        fmpzi_init(u);

        for (i = 0; i < len1; i++)
        {
            fmpz_randbits(fmpzi_realref(a + i), state, bits);
            fmpz_randbits(fmpzi_imagref(a + i), state, bits);
            if (n_randint(state, 2))
                fmpz_neg(fmpzi_realref(a + i), fmpzi_realref(a + i));
            if (n_randint(state, 2))
                fmpz_neg(fmpzi_imagref(a + i), fmpzi_imagref(a + i));
            if (n_randint(state, 8) == 0)
                fmpz_zero(fmpzi_imagref(a + i));
        }
        if (!sq)
            for (i = 0; i < len2; i++)
            {
                fmpz_randbits(fmpzi_realref(b + i), state, bits);
                fmpz_randbits(fmpzi_imagref(b + i), state, bits);
                if (n_randint(state, 2))
                    fmpz_neg(fmpzi_realref(b + i),
                             fmpzi_realref(b + i));
                if (n_randint(state, 2))
                    fmpz_neg(fmpzi_imagref(b + i),
                             fmpzi_imagref(b + i));
            }

        /* either length order */
        if (n_randint(state, 2))
            ok = _fmpzi_poly_mulmid_classical_fft_small(r1, a, len1,
                                                        b, len2,
                                                        nlo, nhi);
        else
            ok = _fmpzi_poly_mulmid_classical_fft_small(r1, b, len2,
                                                        a, len1,
                                                        nlo, nhi);
        if (!ok)
            goto next;   /* gate or representation refusal: fallback */

        /* schoolbook referee */
        for (i = nlo; i < nhi; i++)
        {
            fmpzi_zero(r2 + i - nlo);
            for (j = FLINT_MAX(0, i - len1 + 1);
                 j <= FLINT_MIN(i, len2 - 1); j++)
            {
                fmpzi_mul(t, a + (i - j), b + j);
                fmpzi_add(r2 + i - nlo, r2 + i - nlo, t);
            }
        }

        for (i = 0; i < nhi - nlo; i++)
            if (!fmpzi_equal(r1 + i, r2 + i))
                TEST_FUNCTION_FAIL("len1 = %wd, len2 = %wd, "
                    "bits = %wd, window [%wd, %wd), sq = %d, "
                    "coeff %wd\n", len1, len2, bits, nlo, nhi, sq, i);

next:
        for (i = 0; i < len1; i++)
            fmpzi_clear(a + i);
        if (!sq)
            for (i = 0; i < len2; i++)
                fmpzi_clear(b + i);
        for (i = 0; i < nhi - nlo; i++)
        {
            fmpzi_clear(r1 + i);
            fmpzi_clear(r2 + i);
        }
        fmpzi_clear(t);
        fmpzi_clear(u);
        flint_free(a);
        if (!sq)
            flint_free(b);
        flint_free(r1);
        flint_free(r2);
    }

    /* deterministic coverage: a zero operand yields a zero window;
       and a shared-pointer call with DIFFERENT lengths exercises the
       squaring loop's unpaired terms (partners falling outside the
       kept range) */
#if FLINT_HAVE_FFT_SMALL
    {
        fmpzi_struct a[6], r[5];
        slong i, j;
        int ok;

        for (i = 0; i < 6; i++)
            fmpzi_init(a + i);
        for (i = 0; i < 5; i++)
            fmpzi_init(r + i);

        for (i = 0; i < 6; i++)
        {
            fmpz_randbits(fmpzi_realref(a + i), state, 12000);
            fmpz_randbits(fmpzi_imagref(a + i), state, 12000);
        }

        {
            fmpzi_struct zb[2];
            fmpzi_init(zb + 0);
            fmpzi_init(zb + 1);
            fmpz_one(fmpzi_realref(r + 0));
            ok = _fmpzi_poly_mulmid_classical_fft_small(r, a, 6,
                                                        zb, 2, 1, 4);
            if (!ok || !fmpzi_is_zero(r + 0) || !fmpzi_is_zero(r + 1)
                || !fmpzi_is_zero(r + 2))
                TEST_FUNCTION_FAIL("fmpzi zero operand window\n");
            fmpzi_clear(zb + 0);
            fmpzi_clear(zb + 1);
        }

        ok = _fmpzi_poly_mulmid_classical_fft_small(r, a, 6, a, 3,
                                                    2, 7);
        if (ok)
        {
            fmpzi_t t, s;
            fmpzi_init(t);
            fmpzi_init(s);
            for (i = 2; i < 7; i++)
            {
                fmpzi_zero(s);
                for (j = FLINT_MAX(0, i - 5); j <= FLINT_MIN(i, 2);
                     j++)
                {
                    fmpzi_mul(t, a + (i - j), a + j);
                    fmpzi_add(s, s, t);
                }
                if (!fmpzi_equal(r + i - 2, s))
                    TEST_FUNCTION_FAIL("shared unequal lengths: "
                        "coeff %wd\n", i);
            }
            fmpzi_clear(t);
            fmpzi_clear(s);
        }

        for (i = 0; i < 6; i++)
            fmpzi_clear(a + i);
        for (i = 0; i < 5; i++)
            fmpzi_clear(r + i);
    }
#endif

    TEST_FUNCTION_END(state);
}
