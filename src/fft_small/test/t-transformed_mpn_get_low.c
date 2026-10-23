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
#include "mpn_extras.h"
#include "gr.h"
#include "fft_small.h"

/* the two's-complement residue of (-1)^sign {mag, len} mod B^zn */
static void
_residue(nn_ptr r, slong zn, nn_srcptr mag, slong len, int sign)
{
    flint_mpn_zero(r, zn);
    if (len > 0)
        flint_mpn_copyi(r, mag, FLINT_MIN(len, zn));
    if (sign)
        mpn_neg(r, r, zn);
}

TEST_FUNCTION_START(gr_transformed_mpn_get_low, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong an, bn, i, need, fn, zn, sz, terms;
        int sa, sb, sc, sd, fs, mode, destroy;
        nn_ptr a, b, c, d, full, r1, r2;
        gr_ctx_t T;
        gr_ptr E;

        an = 1 + n_randint(state, n_randint(state, 8) ? 300 : 3000);
        bn = 1 + n_randint(state, an);
        sa = n_randint(state, 2);
        sb = n_randint(state, 2);
        sc = n_randint(state, 2);
        sd = n_randint(state, 2);
        /* 0: a b, 1: a b - c d (mixed signs), 2: a b + c d */
        mode = n_randint(state, 3);
        /* terms_bound 2 gives the two-prime plans below 40000 limbs,
           larger ones the wide plans */
        terms = n_randint(state, 2) ? 2 : 2 + n_randint(state, 6);

        a = flint_malloc(4 * an * sizeof(ulong));
        b = a + an;
        c = b + an;
        d = c + an;
        for (i = 0; i < 4 * an; i++)
            a[i] = n_randtest(state);

        /* values held as small side values (0, 1) or one limb */
        if (n_randint(state, 4) == 0)
        {
            a[0] = n_randint(state, 3);
            for (i = 1; i < an; i++)
                a[i] = 0;
        }

        if (gr_ctx_init_transformed_mpn(T, FLINT_BITS * (an + bn), terms, 1, 5,
                GR_TRANSFORMED_MPN_ALLOC_MALLOC) != GR_SUCCESS)
        {
            flint_free(a);
            continue;
        }

        sz = T->sizeof_elem;
        GR_TMP_INIT_VEC(E, 5, T);

        GR_MUST_SUCCEED(gr_transformed_mpn_set(GR_ENTRY(E, 0, sz), a, an, sa, T));
        GR_MUST_SUCCEED(gr_transformed_mpn_set(GR_ENTRY(E, 1, sz), b, bn, sb, T));
        GR_MUST_SUCCEED(gr_transformed_mpn_set(GR_ENTRY(E, 2, sz), c, an, sc, T));
        GR_MUST_SUCCEED(gr_transformed_mpn_set(GR_ENTRY(E, 3, sz), d, bn, sd, T));
        GR_MUST_SUCCEED(gr_mul(GR_ENTRY(E, 4, sz), GR_ENTRY(E, 0, sz), GR_ENTRY(E, 1, sz), T));
        if (mode == 1)
            GR_MUST_SUCCEED(gr_submul(GR_ENTRY(E, 4, sz), GR_ENTRY(E, 2, sz), GR_ENTRY(E, 3, sz), T));
        else if (mode == 2)
            GR_MUST_SUCCEED(gr_addmul(GR_ENTRY(E, 4, sz), GR_ENTRY(E, 2, sz), GR_ENTRY(E, 3, sz), T));

        need = gr_transformed_mpn_get_limbs(T, GR_ENTRY(E, 4, sz)) + 2;
        full = flint_malloc(need * sizeof(ulong));
        GR_MUST_SUCCEED(gr_transformed_mpn_get(full, need, &fn, &fs, GR_ENTRY(E, 4, sz), T));

        /* windows shorter and longer than the value */
        zn = 1 + n_randint(state, an + bn + 3);
        r1 = flint_malloc(zn * sizeof(ulong));
        r2 = flint_malloc(zn * sizeof(ulong));
        _residue(r1, zn, full, fn, fs);

        destroy = n_randint(state, 2);
        if (destroy)
            GR_MUST_SUCCEED(gr_transformed_mpn_get_low_destructive(r2, zn, GR_ENTRY(E, 4, sz), T));
        else
            GR_MUST_SUCCEED(gr_transformed_mpn_get_low(r2, zn, GR_ENTRY(E, 4, sz), T));

        if (mpn_cmp(r1, r2, zn) != 0)
        {
            flint_printf("FAIL: gr_transformed_mpn_get_low\n");
            flint_printf("an = %wd, bn = %wd, zn = %wd, mode = %d, terms = %wd, destroy = %d\n",
                         an, bn, zn, mode, terms, destroy);
            fflush(stdout);
            flint_abort();
        }

        flint_free(r1);
        flint_free(r2);
        flint_free(full);
        GR_TMP_CLEAR_VEC(E, 5, T);
        gr_ctx_clear(T);
        flint_free(a);
    }

    TEST_FUNCTION_END(state);
}
