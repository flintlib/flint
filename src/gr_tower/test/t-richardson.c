/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    Random positive real numbers in the lazy field, built from rationals,
    square roots, exp, log, sums and products, and identities between
    them which require the relation search.
*/

/* Numbers of moderate size are used as arguments of exp, so that
   enclosures remain representable (exp(exp(exp(3))) is not). */
static void
random_positive(gr_ptr res, flint_rand_t state, slong depth, int small, gr_ctx_t K)
{
    int which = n_randint(state, depth == 0 ? 2 : 6);

    if (small && which == 3)
        which = 4;

    switch (which)
    {
        case 0:
            GR_MUST_SUCCEED(gr_set_ui(res, 1 + n_randint(state, 5), K));
            break;
        case 1:
            {
                fmpq_t q;
                fmpq_init(q);
                fmpq_set_si(q, 1 + n_randint(state, 5), 1 + n_randint(state, 4));
                GR_MUST_SUCCEED(gr_set_fmpq(res, q, K));
                fmpq_clear(q);
            }
            break;
        case 2:
            random_positive(res, state, depth - 1, small, K);
            GR_MUST_SUCCEED(gr_sqrt(res, res, K));
            break;
        case 3:
            random_positive(res, state, depth - 1, 1, K);
            GR_MUST_SUCCEED(gr_exp(res, res, K));
            break;
        case 4:
            {
                /* log of something > 1 */
                gr_ptr t;
                GR_TMP_INIT(t, K);
                random_positive(t, state, depth - 1, small, K);
                GR_MUST_SUCCEED(gr_add_ui(t, t, 1, K));
                GR_MUST_SUCCEED(gr_log(res, t, K));
                GR_TMP_CLEAR(t, K);
            }
            break;
        default:
            {
                gr_ptr t;
                GR_TMP_INIT(t, K);
                random_positive(res, state, depth - 1, small, K);
                random_positive(t, state, depth - 1, small, K);
                if (n_randint(state, 2))
                    GR_MUST_SUCCEED(gr_add(res, res, t, K));
                else
                    GR_MUST_SUCCEED(gr_mul(res, res, t, K));
                GR_TMP_CLEAR(t, K);
            }
            break;
    }
}

TEST_FUNCTION_START(gr_tower_richardson, state)
{
    gr_ctx_t QQ, K;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ptr a, b, x, y, z;
        int which;
        truth_t t;

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(a, b, x, y, z, K);

        random_positive(a, state, 1 + n_randint(state, 3), 1, K);
        random_positive(b, state, 1 + n_randint(state, 3), 1, K);

        which = n_randint(state, 5);

        switch (which)
        {
            case 0:
                /* exp(a) exp(b) = exp(a + b) */
                GR_MUST_SUCCEED(gr_exp(x, a, K));
                GR_MUST_SUCCEED(gr_exp(y, b, K));
                GR_MUST_SUCCEED(gr_mul(x, x, y, K));
                GR_MUST_SUCCEED(gr_add(y, a, b, K));
                GR_MUST_SUCCEED(gr_exp(y, y, K));
                break;
            case 1:
                /* log(a) + log(b) = log(a b) (positive reals) */
                GR_MUST_SUCCEED(gr_log(x, a, K));
                GR_MUST_SUCCEED(gr_log(y, b, K));
                GR_MUST_SUCCEED(gr_add(x, x, y, K));
                GR_MUST_SUCCEED(gr_mul(y, a, b, K));
                GR_MUST_SUCCEED(gr_log(y, y, K));
                break;
            case 2:
                /* exp(log(a)) = a */
                GR_MUST_SUCCEED(gr_log(x, a, K));
                GR_MUST_SUCCEED(gr_exp(x, x, K));
                GR_MUST_SUCCEED(gr_set(y, a, K));
                break;
            case 3:
                /* log(exp(a)) = a (real a) */
                GR_MUST_SUCCEED(gr_exp(x, a, K));
                GR_MUST_SUCCEED(gr_log(x, x, K));
                GR_MUST_SUCCEED(gr_set(y, a, K));
                break;
            default:
                /* exp(a)^n = exp(n a), sqrt(exp(a)) = exp(a/2) */
                {
                    ulong n = 2 + n_randint(state, 3);
                    GR_MUST_SUCCEED(gr_exp(x, a, K));
                    GR_MUST_SUCCEED(gr_pow_ui(x, x, n, K));
                    GR_MUST_SUCCEED(gr_mul_ui(y, a, n, K));
                    GR_MUST_SUCCEED(gr_exp(y, y, K));
                    if (n_randint(state, 2))
                    {
                        GR_MUST_SUCCEED(gr_sqrt(x, x, K));
                        GR_MUST_SUCCEED(gr_div_ui(y, a, 2, K));
                        GR_MUST_SUCCEED(gr_mul_ui(y, y, n, K));
                        GR_MUST_SUCCEED(gr_exp(y, y, K));
                    }
                }
                break;
        }

        t = gr_equal(x, y, K);
        if (t != T_TRUE)
        {
            flint_printf("FAIL: identity %d not proved (%s)\n", which, t == T_FALSE ? "claimed false" : "unknown");
            gr_println(a, K);
            gr_println(b, K);
            gr_println(x, K);
            gr_println(y, K);
            gr_tower_lazy_ctx_stats(K);
            {
                slong level;
                gr_tower_print(gr_tower_lazy_get_tower(&level, x, K));
                gr_tower_print(gr_tower_lazy_get_tower(&level, y, K));
            }
            flint_abort();
        }

        /* a perturbed version is not an identity */
        GR_MUST_SUCCEED(gr_sub(z, x, y, K));
        {
            fmpq_t q;
            fmpq_init(q);
            fmpq_set_si(q, 1 + n_randint(state, 10), 1 + n_randint(state, 10));
            if (n_randint(state, 2))
                fmpq_neg(q, q);
            GR_MUST_SUCCEED(gr_add_fmpq(z, z, q, K));
            fmpq_clear(q);
        }
        if (n_randint(state, 2))
            GR_MUST_SUCCEED(gr_mul(z, z, a, K));

        t = gr_is_zero(z, K);
        if (t != T_FALSE)
        {
            flint_printf("FAIL: perturbed identity %d (%s)\n", which, t == T_TRUE ? "claimed zero" : "unknown");
            gr_println(z, K);
            flint_abort();
        }

        GR_TMP_CLEAR5(a, b, x, y, z, K);
        gr_ctx_clear(K);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
