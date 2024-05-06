/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_vec.h"
#include "arf.h"
#include "nfloat.h"

void
nfloat_test_dot(flint_rand_t state, slong iters, gr_ctx_t ctx)
{
    slong iter, i, len, prec;
    gr_ptr s0, vec1, vec2, res;
    int subtract, initial, reverse;
    arf_ptr avec1, avec2;
    arf_t as0, ares, ares2, amag, err, t;

    prec = NFLOAT_CTX_PREC(ctx);

    for (iter = 0; iter < iters; iter++)
    {
        len = n_randint(state, 30);

        initial = n_randint(state, 2);
        subtract = n_randint(state, 2);
        reverse = n_randint(state, 2);

        vec1 = gr_heap_init_vec(len, ctx);
        vec2 = gr_heap_init_vec(len, ctx);
        s0 = gr_heap_init(ctx);
        res = gr_heap_init(ctx);

        avec1 = _arf_vec_init(len);
        avec2 = _arf_vec_init(len);
        arf_init(as0);
        arf_init(ares);
        arf_init(ares2);
        arf_init(amag);
        arf_init(t);
        arf_init(err);

        if (initial)
        {
            arf_randtest(as0, state, prec, 10);
            GR_MUST_SUCCEED(nfloat_set_arf(s0, as0, ctx));

            arf_set(ares, as0);
            arf_abs(t, as0);
            arf_add(amag, amag, t, prec, ARF_RND_DOWN);
        }

        for (i = 0; i < len; i++)
        {
            if (i > 0 && n_randint(state, 2))
            {
                arf_set(avec1 + i, avec1 + len - 1 - i);
                arf_neg(avec2 + i, avec2 + len - 1 - i);
            }
            else
            {
                arf_randtest(avec1 + i, state, prec, 10);
                arf_randtest(avec2 + i, state, prec, 10);
            }
        }

        for (i = 0; i < len; i++)
        {
            arf_mul(t, avec1 + i, avec2 + (reverse ? len - 1 - i : i), 2 * prec, ARF_RND_DOWN);

            if (subtract)
                arf_sub(ares, ares, t, 2 * prec, ARF_RND_DOWN);
            else
                arf_add(ares, ares, t, 2 * prec, ARF_RND_DOWN);

            arf_abs(t, t);
            arf_add(amag, amag, t, prec, ARF_RND_DOWN);
        }

        /* tolerance */
        arf_mul_2exp_si(t, amag, -prec + 3);

        for (i = 0; i < len; i++)
        {
            GR_MUST_SUCCEED(nfloat_set_arf(GR_ENTRY(vec1, i, ctx->sizeof_elem), avec1 + i, ctx));
            GR_MUST_SUCCEED(nfloat_set_arf(GR_ENTRY(vec2, i, ctx->sizeof_elem), avec2 + i, ctx));
        }

        if (reverse)
            GR_MUST_SUCCEED(_nfloat_vec_dot_rev(res, initial ? s0 : NULL, subtract, vec1, vec2, len, ctx));
        else
            GR_MUST_SUCCEED(_nfloat_vec_dot(res, initial ? s0 : NULL, subtract, vec1, vec2, len, ctx));

        GR_MUST_SUCCEED(nfloat_get_arf(ares2, res, ctx));

        arf_sub(err, ares, ares2, prec, ARF_RND_DOWN);
        arf_abs(err, err);

        if (arf_cmpabs(err, t) > 0)
        {
            flint_printf("FAIL: dot\n");
            gr_ctx_println(ctx);

            flint_printf("reverse = %d, subtract = %d\n", reverse, subtract);

            if (initial)
            {
                flint_printf("\n\ninitial = ");
                arf_printd(as0, 2 + prec / 3.33);
            }

            flint_printf("\n\nvec1 = ");
            _gr_vec_print(vec1, len, ctx);
            flint_printf("\n\nvec2 = ");
            _gr_vec_print(vec2, len, ctx);
            flint_printf("\n\nares = \n");
            arf_printd(ares, 2 + prec / 3.33);
            flint_printf("\n\nares2 = \n");
            arf_printd(ares2, 2 + prec / 3.33);
            flint_printf("\n\ntol = \n");
            arf_printd(t, 10);
            flint_printf("\n\nerr = \n");
            arf_printd(err, 10);
            flint_printf("\n\n");

            flint_abort();
        }

        gr_heap_clear_vec(vec1, len, ctx);
        gr_heap_clear_vec(vec2, len, ctx);
        gr_heap_clear(s0, ctx);
        gr_heap_clear(res, ctx);

        _arf_vec_clear(avec1, len);
        _arf_vec_clear(avec2, len);
        arf_clear(as0);
        arf_clear(ares);
        arf_clear(ares2);
        arf_clear(amag);
        arf_clear(t);
        arf_clear(err);
    }
}

TEST_FUNCTION_START(nfloat, state)
{
    gr_ctx_t ctx;
    slong prec;

    for (prec = NFLOAT_MIN_LIMBS * FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        nfloat_ctx_init(ctx, prec, 0);
        gr_test_floating_point(ctx, 100 * flint_test_multiplier(), 0);
        nfloat_test_dot(state, (prec <= 128 ? 10000 : 100) * flint_test_multiplier(), ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
