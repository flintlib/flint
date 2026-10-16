/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* shared helpers of the mp_real ball tests */

#ifndef MP_REAL_TEST_BALL_HELPERS_H
#define MP_REAL_TEST_BALL_HELPERS_H

/* random normalized mp_real; if allow_err, a random radius (an integer
   count of ulps of the bottom limb, which zero padding may place well
   below the significant limbs) */
static void
mp_real_randtest(mp_real_t x, flint_rand_t state, slong maxsize, int allow_err)
{
    slong size = n_randint(state, maxsize + 1);

    if (size == 0 && n_randint(state, 2))
    {
        mp_real_zero(x);
        x->exp = (slong) n_randint(state, 13) - 6;
        if (allow_err && n_randint(state, 2))
            x->err = 1 + n_randint(state, 1000);
        return;
    }

    size = FLINT_MAX(size, 1);
    mp_real_fit_length(x, size);
    flint_mpn_rrandom(x->d, state, size);
    x->d[size - 1] |= (UWORD(1) << (FLINT_BITS - 1 -
        n_randint(state, FLINT_BITS - 1)));
    x->size = size;
    x->negative = (int) n_randint(state, 2);
    x->exp = size + (slong) n_randint(state, 13) - 6;
    x->err = 0;
    if (allow_err && n_randint(state, 2))
    {
        x->err = 1 + n_randint(state, 1000);
        if (n_randint(state, 2))
            x->err = n_randtest(state) | 1;
        if (n_randint(state, 3) == 0)
        {
            /* zero padding: a radius far below the significant limbs */
            slong pad = 1 + n_randint(state, 4);
            mp_real_fit_length(x, x->size + pad);
            memmove(x->d + pad, x->d, x->size * sizeof(ulong));
            flint_mpn_zero(x->d, pad);
            x->size += pad;
        }
    }
    if (x->err == 0 && n_randint(state, 2))
    {
        /* exercise the low-zero-limb stripping */
        x->d[0] &= ~(ulong) n_randint(state, 2);
    }
    /* keep the normalization invariants */
    if (x->err == 0)
    {
        slong t = 0;
        while (t < x->size && x->d[t] == 0)
            t++;
        if (t == x->size)
        {
            mp_real_zero(x);
            return;
        }
        if (t > 0)
        {
            flint_mpn_copyi(x->d, x->d + t, x->size - t);
            x->size -= t;
        }
    }
}

/* a random true value inside the ball, as an exact arf */
static void
mp_real_random_point(arf_t t, const mp_real_t x, flint_rand_t state)
{
    arb_t b;
    arf_t u;

    arb_init(b);
    arf_init(u);
    mp_real_get_arb(b, x);

    arf_set(t, arb_midref(b));
    if (x->err != 0)
    {
        /* mid + (r / 2^30) * err * ulp, r in [-2^30, 2^30] */
        slong r = (slong) n_randint(state, UWORD(1) << 31)
                    - (slong) (UWORD(1) << 30);
        arf_set_ui(u, x->err);
        arf_mul_si(u, u, r, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul_2exp_si(u, u,
            FLINT_BITS * ((x->size == 0) ? x->exp : x->exp - x->size) - 30);
        arf_add(t, t, u, ARF_PREC_EXACT, ARF_RND_NEAR);
    }

    arb_clear(b);
    arf_clear(u);
}

static void
check_contains(const mp_real_t res, const arf_t truth, const char * op,
    slong iter)
{
    arb_t r;
    arb_init(r);
    mp_real_get_arb(r, res);
    if (!arb_contains_arf(r, truth))
    {
        flint_printf("FAIL: %s (iter %wd)\nres = ", op, iter);
        mp_real_print(res);
        flint_printf("truth = "); arf_printd(truth, 30);
        flint_printf("\n");
        flint_abort();
    }
    arb_clear(r);
}

#endif
