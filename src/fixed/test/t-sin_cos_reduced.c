/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "fixed.h"

/* every algorithm of fixed_sin_cos_reduced agrees with arb within
   the documented budget on both outputs, for random (wn, r) with
   t < 2^-r */

static int
_check(nn_srcptr ysin, nn_srcptr yg, nn_srcptr t, slong wn,
    flint_rand_t state)
{
    arb_t ta, refs, refc, got, bnd;
    fmpz_t f;
    int ok = 1;

    arb_init(ta); arb_init(refs); arb_init(refc); arb_init(got);
    arb_init(bnd); fmpz_init(f);

    fmpz_set_ui_array(f, t, wn);
    arb_set_fmpz(ta, f);
    arb_mul_2exp_si(ta, ta, -FLINT_BITS * wn);
    arb_sin_cos(refs, refc, ta, FLINT_BITS * wn + 128);
    /* g = 1 - cos */
    arb_sub_ui(refc, refc, 1, FLINT_BITS * wn + 128);
    arb_neg(refc, refc);

    arb_set_ui(bnd, FIXED_SIN_COS_REDUCED_MAX_ERR);

    fmpz_set_ui_array(f, ysin, wn);
    arb_set_fmpz(got, f);
    arb_mul_2exp_si(got, got, -FLINT_BITS * wn);
    arb_sub(got, got, refs, FLINT_BITS * wn + 128);
    arb_mul_2exp_si(got, got, FLINT_BITS * wn);
    arb_abs(got, got);
    if (!arb_le(got, bnd))
        ok = 0;

    fmpz_set_ui_array(f, yg, wn);
    arb_set_fmpz(got, f);
    arb_mul_2exp_si(got, got, -FLINT_BITS * wn);
    arb_sub(got, got, refc, FLINT_BITS * wn + 128);
    arb_mul_2exp_si(got, got, FLINT_BITS * wn);
    arb_abs(got, got);
    if (!arb_le(got, bnd))
        ok = (ok == 0) ? 0 : -1;

    arb_clear(ta); arb_clear(refs); arb_clear(refc);
    arb_clear(got); arb_clear(bnd); fmpz_clear(f);
    return ok;
}

TEST_FUNCTION_START(fixed_sin_cos_reduced, state)
{
    slong iter;

    for (iter = 0; iter < 20 + 20 * flint_test_multiplier(); iter++)
    {
        slong wn = 1 + n_randint(state, (iter % 17 == 0) ? 700 : 60);
        flint_bitcnt_t r = 16 + n_randint(state,
            FLINT_MIN(FLINT_BITS * (ulong) wn, 3000));
        slong i, alg;
        nn_ptr t, ysin, yg;

        t = flint_malloc((wn + 2) * sizeof(ulong));
        ysin = flint_malloc((wn + 2) * sizeof(ulong));
        yg = flint_malloc((wn + 2) * sizeof(ulong));

        flint_mpn_urandomb(t, state, FLINT_BITS * wn);
        /* clear the top r bits so t < 2^-r */
        for (i = 0; i < (slong) (r / FLINT_BITS) && i < wn; i++)
            t[wn - 1 - i] = 0;
        if ((r % FLINT_BITS) && (slong) (r / FLINT_BITS) < wn)
            t[wn - 1 - r / FLINT_BITS] >>= (r % FLINT_BITS);

        for (alg = 0; alg <= 4; alg++)
        {
            if (alg == 1 || alg == 2)
                if (r < 32)
                    continue;
            fixed_sin_cos_reduced(ysin, yg, t, wn, r, (int) alg);
            if (_check(ysin, yg, t, wn, state) != 1)
                TEST_FUNCTION_FAIL("alg = %wd, wn = %wd, "
                    "r = %wd\n", alg, wn, (slong) r);
        }

        flint_free(t); flint_free(ysin); flint_free(yg);
    }

    /* targeted regimes that random sampling misses */
    {
        struct { slong wn; flint_bitcnt_t r; int alg; int fill; }
        cases[] = {
            { 64,   16, 0, 0 },   /* r < 32 forces a burst in auto */
            { 48,   16, 4, 1 },   /* whole-limb masks, all-ones */
            { 48,   64, 4, 1 },
            { 48,  128, 3, 1 },
            { 512,  16, 4, 0 },   /* wbot up-shift, sqrt slices */
            { 700,  16, 3, 0 },   /* one burst level + series */
            { 1100, 16, 0, 0 },   /* FULLBURST gate in auto */
            { 2100, 16, 3, 0 },   /* slice sqrt above the Newton
                                     cutoff */
            { 300,  16, 4, 3 },   /* sparse: dead low limbs fold
                                     into the slice frames */
            { 40,  256, 0, 2 },   /* t = 0: sin = g = 0 everywhere */
            { 40,   16, 3, 2 },   /* t = 0 through the burst: every
                                     slice skipped, zero numerator
                                     in the sine division */
            { 40,   16, 4, 2 },
            { 2170, 16, 4, 0 },   /* a second Newton-sqrt geometry
                                     (both parities of the V length
                                     across the two shapes) */
            { 5100, 16, 3, 0 },   /* N > 10000: the high-2-adic
                                     re-padding of the term count */
            { 60, 3820, 3, 0 },   /* r within a limb of 64 wn: the
                                     degenerate single-slice
                                     ladder */
        };
        slong c, i;

        for (c = 0; c < (slong) (sizeof(cases) / sizeof(cases[0]));
            c++)
        {
            slong wn = cases[c].wn;
            flint_bitcnt_t r = cases[c].r;
            nn_ptr t, ysin, yg;

            t = flint_malloc((wn + 2) * sizeof(ulong));
            ysin = flint_malloc((wn + 2) * sizeof(ulong));
            yg = flint_malloc((wn + 2) * sizeof(ulong));

            if (cases[c].fill == 2)
                flint_mpn_zero(t, wn);
            else if (cases[c].fill == 1)
                flint_mpn_store(t, wn, ~UWORD(0));
            else if (cases[c].fill == 3)
            {
                flint_mpn_zero(t, wn);
                for (i = 0; i < wn; i += 37)
                    t[i] = n_randtest(state);
                t[wn - 1] = n_randtest(state);
            }
            else
                flint_mpn_urandomb(t, state, FLINT_BITS * wn);
            for (i = 0; i < (slong) (r / FLINT_BITS) && i < wn; i++)
                t[wn - 1 - i] = 0;
            if ((r % FLINT_BITS) && (slong) (r / FLINT_BITS) < wn)
                t[wn - 1 - r / FLINT_BITS] >>= (r % FLINT_BITS);

            fixed_sin_cos_reduced(ysin, yg, t, wn, r, cases[c].alg);
            if (_check(ysin, yg, t, wn, state) != 1)
                TEST_FUNCTION_FAIL("targeted case %wd: wn = %wd, "
                    "r = %wd, alg = %d\n", c, wn, (slong) r,
                    cases[c].alg);

            flint_free(t); flint_free(ysin); flint_free(yg);
        }
    }

    TEST_FUNCTION_END(state);
}
