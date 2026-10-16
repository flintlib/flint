/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "acb_dft.h"

/* The module is a thin wrapper around gr_dft, whose test suite
   stress-tests the algorithms; here we verify that every public
   function behaves as documented (definition, inverse, plan reuse,
   aliasing, trivial lengths). */

/* naive X_k = sum_j v_j e(-jk/n) directly from the definition */
static void
_dft_naive(acb_ptr w, acb_srcptr v, slong n, slong prec)
{
    slong j, k;
    acb_ptr roots;
    acb_t t;

    roots = _acb_vec_init(n);
    acb_init(t);
    for (j = 0; j < n; j++)
    {
        arb_set_si(acb_realref(t), -2 * j);
        arb_div_ui(acb_realref(t), acb_realref(t), n, prec + 10);
        arb_zero(acb_imagref(t));
        acb_exp_pi_i(roots + j, t, prec + 10);
    }

    for (k = 0; k < n; k++)
    {
        acb_zero(w + k);
        for (j = 0; j < n; j++)
            acb_addmul(w + k, roots + (j * k) % n, v + j, prec);
    }

    _acb_vec_clear(roots, n);
    acb_clear(t);
}

TEST_FUNCTION_START(acb_dft, state)
{
    slong k;
    slong lens[9] = { 1, 2, 3, 4, 6, 8, 12, 30, 59 };

    for (k = 0; k < 9; k++)
    {
        slong len = lens[k], i;
        slong prec = 64 + n_randint(state, 128);
        acb_ptr v, w1, w2;
        acb_dft_pre_t pre;

        v = _acb_vec_init(len);
        w1 = _acb_vec_init(len);
        w2 = _acb_vec_init(len);

        for (i = 0; i < len; i++)
            acb_randtest_precise(v + i, state, prec, 0);

        /* the one-shot transform computes the DFT */
        _dft_naive(w1, v, len, prec);
        acb_dft(w2, v, len, prec);
        for (i = 0; i < len; i++)
            if (!acb_overlaps(w1 + i, w2 + i))
                TEST_FUNCTION_FAIL("dft vs definition, len = %wd i = %wd\n",
                        len, i);

        /* inverse roundtrip, aliased */
        _acb_vec_set(w2, v, len);
        acb_dft(w2, w2, len, prec);
        acb_dft_inverse(w2, w2, len, prec);
        for (i = 0; i < len; i++)
            if (!acb_overlaps(w2 + i, v + i))
                TEST_FUNCTION_FAIL("roundtrip, len = %wd i = %wd\n", len, i);

        /* the plan interface agrees with the one-shots and is
           reusable */
        acb_dft_precomp_init(pre, len, prec);
        acb_dft_precomp(w2, v, pre, prec);
        for (i = 0; i < len; i++)
            if (!acb_overlaps(w1 + i, w2 + i))
                TEST_FUNCTION_FAIL("precomp, len = %wd i = %wd\n", len, i);
        acb_dft_inverse_precomp(w2, w2, pre, prec);
        for (i = 0; i < len; i++)
            if (!acb_overlaps(w2 + i, v + i))
                TEST_FUNCTION_FAIL("precomp roundtrip, len = %wd i = %wd\n",
                        len, i);
        acb_dft_precomp_clear(pre);

        _acb_vec_clear(v, len);
        _acb_vec_clear(w1, len);
        _acb_vec_clear(w2, len);
    }

    /* len = 0 is a no-op for every entry point */
    {
        acb_dft_pre_t pre;
        acb_dft(NULL, NULL, 0, 64);
        acb_dft_inverse(NULL, NULL, 0, 64);
        acb_dft_precomp_init(pre, 0, 64);
        acb_dft_precomp(NULL, NULL, pre, 64);
        acb_dft_inverse_precomp(NULL, NULL, pre, 64);
        acb_dft_precomp_clear(pre);
    }

    TEST_FUNCTION_END(state);
}
