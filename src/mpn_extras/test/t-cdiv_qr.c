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

TEST_FUNCTION_START(flint_mpn_cdiv_qr, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, b, q, r, q2, r2;
        mp_size_t an, bn, n;
        fmpz_t fa, fb, fq, fr;

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
        if (n_randint(state, 5) == 0)
        {
            /* exact case */
            mp_ptr t = flint_malloc((n + bn) * sizeof(mp_limb_t));
            if (n >= bn)
                flint_mpn_mul(t, a, n, b, bn);
            else
                flint_mpn_mul(t, b, bn, a, n);
            flint_mpn_copyi(a, t, an);
            flint_free(t);
        }

        fmpz_init(fa); fmpz_init(fb); fmpz_init(fq); fmpz_init(fr);
        fmpz_set_ui_array(fa, a, an);
        fmpz_set_ui_array(fb, b, bn);
        fmpz_cdiv_qr(fq, fr, fa, fb);
        fmpz_neg(fr, fr);
        fmpz_get_ui_array(q2, n + 1, fq);
        fmpz_get_ui_array(r2, bn, fr);

        flint_mpn_cdiv_qr(q, r, a, an, b, bn);

        if (mpn_cmp(q, q2, n + 1) != 0 || mpn_cmp(r, r2, bn) != 0)
            TEST_FUNCTION_FAIL("cdiv_qr: an = %wd, bn = %wd\n", an, bn);

        flint_mpn_cdiv_q(q, a, an, b, bn);
        flint_mpn_cdiv_r(r, a, an, b, bn);

        if (mpn_cmp(q, q2, n + 1) != 0 || mpn_cmp(r, r2, bn) != 0)
            TEST_FUNCTION_FAIL("cdiv_q/cdiv_r: an = %wd, bn = %wd\n", an, bn);

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
