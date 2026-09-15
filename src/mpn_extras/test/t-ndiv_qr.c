/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(flint_mpn_ndiv_qr, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, b, q, r, q2, r2;
        mp_size_t an, bn, n;
        fmpz_t fa, fb, fq, fr;
        int sgn, sgn2;

        bn = 1 + n_randint(state, 30);
        an = bn + n_randint(state, 40);
        if (n_randint(state, 30) == 0)
        {
            bn = 1 + n_randint(state, 600);
            an = bn + n_randint(state, 1000);
        }
        n = an - bn + 1;

        a = flint_malloc(an * sizeof(mp_limb_t));
        b = flint_malloc(bn * sizeof(mp_limb_t));
        q = flint_malloc((n + 1) * sizeof(mp_limb_t));
        r = flint_malloc(bn * sizeof(mp_limb_t));
        q2 = flint_malloc((n + 1) * sizeof(mp_limb_t));
        r2 = flint_malloc(bn * sizeof(mp_limb_t));

        flint_mpn_rrandom(a, state, an);
        flint_mpn_rrandom(b, state, bn);
        if (b[bn - 1] == 0)
            b[bn - 1] = 1;
        if (n_randint(state, 3) == 0)
        {
            /* exact or half-way cases */
            mp_ptr t = flint_malloc((n + bn) * sizeof(mp_limb_t));
            if (n >= bn)
                flint_mpn_mul(t, a, n, b, bn);
            else
                flint_mpn_mul(t, b, bn, a, n);
            flint_mpn_copyi(a, t, an);
            if (n_randint(state, 2))
            {
                b[0] &= ~UWORD(1);
                if (bn == 1 && b[0] == 0)
                    b[0] = 2;
                mpn_rshift(t, b, bn, 1);
                mpn_add(a, a, an, t, bn);
            }
            flint_free(t);
        }

        fmpz_init(fa); fmpz_init(fb); fmpz_init(fq); fmpz_init(fr);
        fmpz_set_ui_array(fa, a, an);
        fmpz_set_ui_array(fb, b, bn);
        /* round to nearest, ties to even (fmpz_ndiv_qr breaks ties
           towards zero) */
        fmpz_tdiv_qr(fq, fr, fa, fb);
        {
            int c = fmpz_cmp2abs(fb, fr);   /* compare |b| with 2|r| */
            if (c < 0 || (c == 0 && fmpz_is_odd(fq)))
            {
                fmpz_add_ui(fq, fq, 1);
                fmpz_sub(fr, fr, fb);
            }
        }
        sgn2 = fmpz_sgn(fr);
        fmpz_abs(fr, fr);
        fmpz_get_ui_array(q2, n + 1, fq);
        fmpz_get_ui_array(r2, bn, fr);

        sgn = flint_mpn_ndiv_qr(q, r, a, an, b, bn);

        if (sgn != sgn2 || mpn_cmp(q, q2, n + 1) != 0 || mpn_cmp(r, r2, bn) != 0)
            TEST_FUNCTION_FAIL("ndiv_qr: an = %wd, bn = %wd, sgn = %d, sgn2 = %d\n", an, bn, sgn, sgn2);

        flint_mpn_ndiv_q(q, a, an, b, bn);
        sgn = flint_mpn_ndiv_r(r, a, an, b, bn);

        if (sgn != sgn2 || mpn_cmp(q, q2, n + 1) != 0 || mpn_cmp(r, r2, bn) != 0)
            TEST_FUNCTION_FAIL("ndiv_q/ndiv_r: an = %wd, bn = %wd\n", an, bn);

        fmpz_clear(fa); fmpz_clear(fb); fmpz_clear(fq); fmpz_clear(fr);
        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(r);
        flint_free(q2);
        flint_free(r2);
    }

    TEST_FUNCTION_END(state);
}
