/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"
#include "ball_helpers.h"

/* the complex product against the exact arf products, with the parts
   of an operand at unrelated magnitudes and every aliasing */
static void
test_mul_complex(flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        mp_real_t ar, ai, br, bi, rr, ri, t1, t2;
        arf_t xr, xi, yr, yi, zr, zi, u;
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? 400 : 12);
        slong maxsize = (iter % 10 == 0) ? 400 : 10;
        int alias = (int) n_randint(state, 5);

        mp_real_init(ar); mp_real_init(ai); mp_real_init(br); mp_real_init(bi);
        mp_real_init(rr); mp_real_init(ri); mp_real_init(t1); mp_real_init(t2);
        arf_init(xr); arf_init(xi); arf_init(yr); arf_init(yi);
        arf_init(zr); arf_init(zi); arf_init(u);

        mp_real_randtest(ar, state, maxsize, 1);
        mp_real_randtest(ai, state, maxsize, 1);
        mp_real_randtest(br, state, maxsize, 1);
        mp_real_randtest(bi, state, maxsize, 1);
        /* often a tiny or an exact-zero imaginary part, as in the
           accumulation of exp(i x) for small x */
        if (n_randint(state, 3) == 0)
            mp_real_mul_2exp_si(ai, ai, -(slong) FLINT_BITS * (slong) n_randint(state, 2 * maxsize + 2));
        if (n_randint(state, 5) == 0)
            mp_real_zero(ai);
        if (n_randint(state, 3) == 0)
            mp_real_mul_2exp_si(bi, bi, -(slong) FLINT_BITS * (slong) n_randint(state, 2 * maxsize + 2));
        if (n_randint(state, 5) == 0)
            mp_real_zero(bi);
        mp_real_random_point(xr, ar, state);
        mp_real_random_point(xi, ai, state);
        mp_real_random_point(yr, br, state);
        mp_real_random_point(yi, bi, state);

        arf_mul(zr, xr, yr, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul(u, xi, yi, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_sub(zr, zr, u, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul(zi, xr, yi, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul(u, xi, yr, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_add(zi, zi, u, ARF_PREC_EXACT, ARF_RND_NEAR);

        if (alias == 0)
            mp_real_mul_complex(rr, ri, ar, ai, br, bi, n);
        else if (alias == 1)
        {
            mp_real_set(rr, ar); mp_real_set(ri, ai);
            mp_real_mul_complex(rr, ri, rr, ri, br, bi, n);
        }
        else if (alias == 2)
        {
            mp_real_set(rr, br); mp_real_set(ri, bi);
            mp_real_mul_complex(rr, ri, ar, ai, rr, ri, n);
        }
        else if (alias == 3)
        {
            /* a square */
            mp_real_set(br, ar); mp_real_set(bi, ai);
            arf_set(yr, xr); arf_set(yi, xi);
            arf_mul(zr, xr, yr, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_mul(u, xi, yi, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_sub(zr, zr, u, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_mul(zi, xr, yi, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_mul(u, xi, yr, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_add(zi, zi, u, ARF_PREC_EXACT, ARF_RND_NEAR);
            mp_real_mul_complex(rr, ri, ar, ai, ar, ai, n);
        }
        else
        {
            mp_real_set(rr, ai); mp_real_set(ri, ar);
            mp_real_mul_complex(rr, ri, ri, rr, br, bi, n);
        }

        check_contains(rr, zr, "mul_complex re", iter);
        check_contains(ri, zi, "mul_complex im", iter);

        /* accuracy: about n limbs relative to the larger part of the
           result (one frame for both parts), or that of the operands:
           against the four products, the radius of each part is
           within a few limbs of the reference's or of the frame */
        {
            mp_real_t t3, t4;
            double lr, li, ref_r, ref_i, top, fr;
            mp_real_init(t3); mp_real_init(t4);
            mp_real_mul(t1, ar, br, n);
            mp_real_mul(t2, ai, bi, n);
            mp_real_sub(t1, t1, t2, n);
            mp_real_mul(t3, ar, bi, n);
            mp_real_mul(t4, ai, br, n);
            mp_real_add(t3, t3, t4, n);
            /* log2 of a radius, and of the frame's ulp */
#define RADBITS(x) ((x)->err == 0 ? -1e300 : (double) FLINT_BITS * ((x)->exp - (x)->size) + FLINT_BIT_COUNT((x)->err))
#define TOPBITS(x) ((x)->size == 0 ? -1e300 : (double) FLINT_BITS * (x)->exp)
            lr = RADBITS(rr); li = RADBITS(ri);
            ref_r = RADBITS(t1); ref_i = RADBITS(t3);
            top = FLINT_MAX(TOPBITS(t1), TOPBITS(t3));
            fr = top - FLINT_BITS * (n - 2);
            if (lr > FLINT_MAX(FLINT_MAX(ref_r, ref_i), fr) + 3 * FLINT_BITS
                || li > FLINT_MAX(FLINT_MAX(ref_r, ref_i), fr) + 3 * FLINT_BITS)
            {
                flint_printf("FAIL: mul_complex accuracy (iter %wd): rad %g %g ref %g %g frame %g\n", iter, lr, li, ref_r, ref_i, fr);
                mp_real_print(ar); mp_real_print(ai); mp_real_print(br); mp_real_print(bi);
                mp_real_print(rr); mp_real_print(ri); mp_real_print(t1); mp_real_print(t3);
                flint_abort();
            }
            mp_real_clear(t3); mp_real_clear(t4);
        }

        mp_real_clear(ar); mp_real_clear(ai); mp_real_clear(br); mp_real_clear(bi);
        mp_real_clear(rr); mp_real_clear(ri); mp_real_clear(t1); mp_real_clear(t2);
        arf_clear(xr); arf_clear(xi); arf_clear(yr); arf_clear(yi);
        arf_clear(zr); arf_clear(zi); arf_clear(u);
    }
}

TEST_FUNCTION_START(mp_real_mul_complex, state)
{
    test_mul_complex(state, 2000 + 2000 * flint_test_multiplier());

    TEST_FUNCTION_END(state);
}
