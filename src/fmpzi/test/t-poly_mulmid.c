/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpzi.h"

/* reference: (ar + i ai)(br + i bi) from four real products */
static void
_fmpzi_poly_mulmid_ref(fmpzi_struct * res, const fmpzi_struct * a, slong len1,
                       const fmpzi_struct * b, slong len2, slong nlo, slong nhi)
{
    fmpz_poly_t ar, ai, br, bi, t, u;
    slong i;

    fmpz_poly_init(ar); fmpz_poly_init(ai);
    fmpz_poly_init(br); fmpz_poly_init(bi);
    fmpz_poly_init(t); fmpz_poly_init(u);

    for (i = 0; i < len1; i++)
    {
        fmpz_poly_set_coeff_fmpz(ar, i, fmpzi_realref(a + i));
        fmpz_poly_set_coeff_fmpz(ai, i, fmpzi_imagref(a + i));
    }
    for (i = 0; i < len2; i++)
    {
        fmpz_poly_set_coeff_fmpz(br, i, fmpzi_realref(b + i));
        fmpz_poly_set_coeff_fmpz(bi, i, fmpzi_imagref(b + i));
    }

    fmpz_poly_mul(t, ar, br);
    fmpz_poly_mul(u, ai, bi);
    fmpz_poly_sub(t, t, u);
    for (i = nlo; i < nhi; i++)
        fmpz_poly_get_coeff_fmpz(fmpzi_realref(res + i - nlo), t, i);

    fmpz_poly_mul(t, ar, bi);
    fmpz_poly_mul(u, ai, br);
    fmpz_poly_add(t, t, u);
    for (i = nlo; i < nhi; i++)
        fmpz_poly_get_coeff_fmpz(fmpzi_imagref(res + i - nlo), t, i);

    fmpz_poly_clear(ar); fmpz_poly_clear(ai);
    fmpz_poly_clear(br); fmpz_poly_clear(bi);
    fmpz_poly_clear(t); fmpz_poly_clear(u);
}

static fmpzi_struct *
_vec_init(slong len)
{
    slong i;
    fmpzi_struct * v = flint_malloc(len * sizeof(fmpzi_struct));
    for (i = 0; i < len; i++)
        fmpzi_init(v + i);
    return v;
}

static void
_vec_clear(fmpzi_struct * v, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        fmpzi_clear(v + i);
    flint_free(v);
}

static void
_fmpzi_vec_randtest_parts(fmpzi_struct * v, flint_rand_t state, slong len,
                          slong rbits, slong ibits)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        fmpz_randtest(fmpzi_realref(v + i), state, rbits);
        fmpz_randtest(fmpzi_imagref(v + i), state, ibits);
    }
}

/* coefficients of maximal size M = 2^bits - 1 with signs making both
   parts of every product term add up: (M + iM)(M - iM) = 2M^2, etc. */
static void
_fmpzi_vec_extreme(fmpzi_struct * v, slong len, slong bits, int conj, int neg)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        fmpz_one(fmpzi_realref(v + i));
        fmpz_mul_2exp(fmpzi_realref(v + i), fmpzi_realref(v + i), bits);
        fmpz_sub_ui(fmpzi_realref(v + i), fmpzi_realref(v + i), 1);
        if (neg)
            fmpz_neg(fmpzi_realref(v + i), fmpzi_realref(v + i));
        if (conj)
            fmpz_neg(fmpzi_imagref(v + i), fmpzi_realref(v + i));
        else
            fmpz_set(fmpzi_imagref(v + i), fmpzi_realref(v + i));
    }
}

typedef void (*fmpzi_mulmid_func)(fmpzi_struct *, const fmpzi_struct *, slong,
                                  const fmpzi_struct *, slong, slong, slong);

TEST_FUNCTION_START(fmpzi_poly_mulmid, state)
{
    slong iter;

    for (iter = 0; iter < 900 * flint_test_multiplier(); iter++)
    {
        int which = n_randint(state, 5);
        fmpzi_mulmid_func funcs[5] = { _fmpzi_poly_mulmid, _fmpzi_poly_mulmid_KS,
            _fmpzi_poly_mulmid_SS, _fmpzi_poly_mulmid_SS_pair, _fmpzi_poly_mulmid_SS_single };
        const char * names[5] = { "", "_KS", "_SS", "_SS_pair", "_SS_single" };
        fmpzi_mulmid_func func = funcs[which];
        const char * name = names[which];
        slong len1, len2, nlo, nhi, len, i, bits;
        int squaring, alias, ok;
        fmpzi_struct * a, * b, * r1, * r2;

        if (n_randint(state, 4) == 0)
            len1 = 1 + n_randint(state, 800);
        else
            len1 = 1 + n_randint(state, 40);

        /* up to 3000 bits, and short inputs from 8000 bits, where the
           dispatcher reaches the classical algorithm on transforms */
        if (n_randint(state, 8) == 0)
        {
            slong len_max = 1 + n_randint(state, 30);

            bits = 8000 + n_randint(state, 12000);
            len1 = FLINT_MIN(len1, len_max);
        }
        else
            bits = 1 + n_randint(state, n_randint(state, 4) == 0 ? 3000 : 200);

        len2 = 1 + n_randint(state, n_randint(state, 2) ? len1 : 40);

        squaring = (n_randint(state, 5) == 0);
        if (squaring)
            len2 = len1;

        a = _vec_init(len1);
        b = squaring ? a : _vec_init(len2);

        /* sometimes unbalanced or zero real/imaginary parts, or
           coefficients of extreme size */
        switch (n_randint(state, 6))
        {
            case 5:
                _fmpzi_vec_extreme(a, len1, bits, 0, n_randint(state, 2));
                break;
            case 0:
                _fmpzi_vec_randtest_parts(a, state, len1, bits, 0);
                break;
            case 1:
                _fmpzi_vec_randtest_parts(a, state, len1, 0, bits);
                break;
            case 2:
                _fmpzi_vec_randtest_parts(a, state, len1, bits, 1 + bits / 8);
                break;
            default:
                _fmpzi_vec_randtest_parts(a, state, len1, bits, bits);
        }

        if (!squaring)
        {
            slong bits2 = 1 + n_randint(state, bits + 10);

            if (n_randint(state, 6) == 0)
                _fmpzi_vec_extreme(b, len2, bits2, n_randint(state, 2), n_randint(state, 2));
            else
            {
                /* sometimes real or purely imaginary */
                int k = n_randint(state, 6);
                _fmpzi_vec_randtest_parts(b, state, len2, (k == 0) ? 0 : bits2,
                                          (k == 1) ? 0 : bits2);
            }
        }

        len = len1 + len2 - 1;

        switch (n_randint(state, 3))
        {
            case 0:
                /* the classical middle product */
                nlo = FLINT_MIN(len1, len2) - 1;
                nhi = FLINT_MAX(len1, len2);
                break;
            case 1:
                nlo = 0;
                nhi = 1 + n_randint(state, len);
                break;
            default:
                nlo = n_randint(state, len);
                nhi = nlo + 1 + n_randint(state, len - nlo);
                break;
        }

        r1 = _vec_init(FLINT_MAX(nhi - nlo, len1));
        r2 = _vec_init(nhi - nlo);

        _fmpzi_poly_mulmid_ref(r2, a, len1, b, len2, nlo, nhi);

        /* the dispatcher, like the gr functions, need not support res
           aliasing an input */
        alias = which != 0 && !squaring && (nhi - nlo <= len1) && n_randint(state, 4) == 0;

        if (alias)
        {
            /* res aliasing the first input */
            for (i = 0; i < len1; i++)
                fmpzi_set(r1 + i, a + i);
            func(r1, r1, len1, b, len2, nlo, nhi);
        }
        else if (!squaring && n_randint(state, 2))
            func(r1, b, len2, a, len1, nlo, nhi);
        else
            func(r1, a, len1, b, len2, nlo, nhi);

        ok = 1;
        for (i = 0; i < nhi - nlo; i++)
            ok = ok && fmpzi_equal(r1 + i, r2 + i);

        if (!ok)
        {
            flint_printf("FAIL: fmpzi_poly_mulmid%s\n", name);
            flint_printf("len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd, bits = %wd, squaring = %d, alias = %d\n",
                len1, len2, nlo, nhi, bits, squaring, alias);
            for (i = 0; i < nhi - nlo; i++)
            {
                if (!fmpzi_equal(r1 + i, r2 + i))
                {
                    flint_printf("first difference at %wd\n", i + nlo);
                    break;
                }
            }
            fflush(stdout);
            flint_abort();
        }

        _vec_clear(a, len1);
        if (!squaring)
            _vec_clear(b, len2);
        _vec_clear(r1, FLINT_MAX(nhi - nlo, len1));
        _vec_clear(r2, nhi - nlo);
    }

    TEST_FUNCTION_END(state);
}
