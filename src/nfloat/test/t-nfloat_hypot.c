/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "gr.h"
#include "nfloat.h"

/* Checks that nfloat_hypot and nfloat_complex_abs are correctly rounded in
   every rounding mode (truncated in the default mode), against the square
   root of the exact sum of squares in arf. Inputs include short mantissas
   (exact results such as hypot(3, 4) = 5 and near-ties), exponent
   differences around the threshold of the tiny case, and aliasing. */

static void
randtest_mant(nfloat_ptr x, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);

    for (i = 0; i < n; i++)
        NFLOAT_D(x)[i] = n_randtest(state);

    if (n_randint(state, 4) == 0)
    {
        /* few bits */
        flint_mpn_zero(NFLOAT_D(x), n - 1);
        NFLOAT_D(x)[n - 1] &= ~((UWORD(1) << n_randint(state, FLINT_BITS)) - 1);
    }

    NFLOAT_D(x)[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    NFLOAT_SGNBIT(x) = n_randint(state, 2);
}

TEST_FUNCTION_START(nfloat_hypot, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, cctx;
        slong n, prec, d;
        int rnd, flags, aliasing;
        nfloat_ptr x, y, z;
        nfloat_complex_ptr w;
        arf_t a, b, s, r, zf;
        arf_rnd_t arnd;

        if (n_randint(state, 8) == 0)
            n = 1 + n_randint(state, NFLOAT_MAX_LIMBS);
        else
            n = 1 + n_randint(state, 4);

        prec = n * FLINT_BITS;
        rnd = n_randint(state, 3);
        flags = (rnd == 0) ? 0 : (rnd == 1) ? NFLOAT_RND_FLOOR : NFLOAT_RND_CEIL;
        arnd = (rnd == 0) ? ARF_RND_DOWN : (rnd == 1) ? ARF_RND_FLOOR : ARF_RND_CEIL;
        nfloat_ctx_init(ctx, prec, flags);
        nfloat_complex_ctx_init(cctx, prec, flags);

        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);
        w = gr_heap_init(cctx);
        arf_init(a);
        arf_init(b);
        arf_init(s);
        arf_init(r);
        arf_init(zf);

        randtest_mant(x, state, ctx);
        randtest_mant(y, state, ctx);
        NFLOAT_EXP(x) = (slong) n_randint(state, 200) - 100;
        /* exponent differences up to past the threshold of the tiny case */
        d = n_randint(state, FLINT_BITS / 2 * n + 8);
        if (n_randint(state, 4) == 0)
            d = n_randint(state, 4);
        NFLOAT_EXP(y) = NFLOAT_EXP(x) + (n_randint(state, 2) ? d : -d);
        if (n_randint(state, 16) == 0)
            GR_MUST_SUCCEED(nfloat_zero(y, ctx));

        aliasing = n_randint(state, 3);
        if (aliasing == 1)
        {
            GR_MUST_SUCCEED(nfloat_set(z, x, ctx));
            GR_MUST_SUCCEED(nfloat_hypot(z, z, y, ctx));
        }
        else if (aliasing == 2)
        {
            GR_MUST_SUCCEED(nfloat_set(z, y, ctx));
            GR_MUST_SUCCEED(nfloat_hypot(z, x, z, ctx));
        }
        else
            GR_MUST_SUCCEED(nfloat_hypot(z, x, y, ctx));

        GR_MUST_SUCCEED(nfloat_get_arf(a, x, ctx));
        GR_MUST_SUCCEED(nfloat_get_arf(b, y, ctx));
        arf_mul(s, a, a, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_addmul(s, b, b, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_sqrt(r, s, prec, arnd);
        GR_MUST_SUCCEED(nfloat_get_arf(zf, z, ctx));

        if (!arf_equal(zf, r))
        {
            flint_printf("FAIL (hypot): not correctly rounded\n");
            flint_printf("n = %wd, rnd = %d, aliasing = %d\n", n, rnd, aliasing);
            flint_printf("x = "); arf_printd(a, prec / 3 + 5); flint_printf("\n");
            flint_printf("y = "); arf_printd(b, prec / 3 + 5); flint_printf("\n");
            flint_printf("z = "); arf_printd(zf, prec / 3 + 5); flint_printf("\n");
            flint_printf("r = "); arf_printd(r, prec / 3 + 5); flint_printf("\n");
            flint_abort();
        }

        /* complex abs */
        GR_MUST_SUCCEED(nfloat_set(NFLOAT_COMPLEX_RE(w, cctx), x, ctx));
        GR_MUST_SUCCEED(nfloat_set(NFLOAT_COMPLEX_IM(w, cctx), y, ctx));
        GR_MUST_SUCCEED(nfloat_complex_abs(w, w, cctx));
        GR_MUST_SUCCEED(nfloat_get_arf(zf, NFLOAT_COMPLEX_RE(w, cctx), ctx));

        if (!arf_equal(zf, r) || !NFLOAT_IS_ZERO(NFLOAT_COMPLEX_IM(w, cctx)))
        {
            flint_printf("FAIL (complex_abs)\n");
            flint_printf("n = %wd, rnd = %d\n", n, rnd);
            flint_abort();
        }

        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);
        gr_heap_clear(w, cctx);
        arf_clear(a);
        arf_clear(b);
        arf_clear(s);
        arf_clear(r);
        arf_clear(zf);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
