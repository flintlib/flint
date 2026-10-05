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
#include "qqbar.h"
#include "acb.h"
#include "gr_vec.h"
#include "gr_tower.h"

/* checks that the enclosure of the flat element x agrees with z */
static void
check_acb(gr_tower_t T, gr_srcptr x, gr_ctx_t F, const acb_t z, const char * what)
{
    acb_t w;
    acb_init(w);
    GR_MUST_SUCCEED(gr_tower_field_flat_get_acb(w, x, 64, F));
    if (!acb_overlaps(w, z) || acb_rel_accuracy_bits(w) < 40)
    {
        flint_printf("FAIL: %s\n", what);
        gr_println(x, F);
        acb_printn(w, 20, 0); flint_printf("\n");
        acb_printn(z, 20, 0); flint_printf("\n");
        flint_abort();
    }
    acb_clear(w);
}

TEST_FUNCTION_START(gr_tower_transcendental, state)
{
    gr_ctx_t QQ;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    /* Q(pi), Q(e), and the growth of the flat context beyond its capacity */
    {
        gr_tower_t T;
        gr_ctx_t F;
        gr_ptr x, y;
        acb_t z, w;
        slong j;

        acb_init(z);
        acb_init(w);
        gr_tower_init(T, QQ);

        GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));

        /* e = exp(1), then exp(k) and log(k) for k = 2, 3, 4: eight
           generators in total (pi, e, exp(2), log(2), ...), exceeding
           the initial capacity of the flat context */
        for (j = 1; j <= 4; j++)
        {
            fmpz_mpoly_q_t c;
            fmpz_mpoly_ctx_struct * mctx = gr_tower_flat(T)->mctx;

            fmpz_mpoly_q_init(c, mctx);
            fmpz_mpoly_q_set_si(c, j, mctx);
            GR_MUST_SUCCEED(gr_tower_adjoin_exp_flat(T, c, mctx, j == 1 ? "e" : NULL));
            if (j > 1)
                GR_MUST_SUCCEED(gr_tower_adjoin_log_flat(T, c, mctx, NULL));
            fmpz_mpoly_q_clear(c, mctx);
        }

        /* log(0) is rejected */
        {
            gr_ptr c = gr_heap_init(gr_tower_field(T));
            GR_MUST_SUCCEED(gr_zero(c, gr_tower_field(T)));
            if (gr_tower_adjoin_log(T, c, NULL) != GR_DOMAIN)
            {
                flint_printf("FAIL: log(0)\n");
                flint_abort();
            }
            gr_heap_clear(c, gr_tower_field(T));
        }

        if (T->num_trans != 8 || T->length != 0)
        {
            flint_printf("FAIL: number of generators\n");
            flint_abort();
        }

        /* enclosures of the generators */
        GR_MUST_SUCCEED(gr_tower_trans_get_acb(z, T, 1, 100));
        acb_const_pi(w, 100);
        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < 90)
        {
            flint_printf("FAIL: pi\n");
            flint_abort();
        }
        for (j = 1; j <= 4; j++)
        {
            GR_MUST_SUCCEED(gr_tower_trans_get_acb(z, T, (j == 1) ? 2 : 2 * j - 1, 100));
            acb_set_ui(w, j); acb_exp(w, w, 100);
            if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < 90)
            {
                flint_printf("FAIL: exp(%wd)\n", j);
                gr_tower_print(T);
                acb_printn(z, 30, 0); flint_printf("\n");
                acb_printn(w, 30, 0); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
            if (j == 1)
                continue;
            GR_MUST_SUCCEED(gr_tower_trans_get_acb(z, T, 2 * j, 100));
            acb_set_ui(w, j); acb_log(w, w, 100);
            if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < 90)
            {
                flint_printf("FAIL: log(%wd)\n", j);
                flint_abort();
            }
        }

        /* arithmetic in Q(pi, e, ...) as a fixed flat field */
        gr_ctx_init_tower_field_flat(F, T);
        gr_test_ring(F, 5 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);

        x = gr_heap_init(F);
        y = gr_heap_init(F);
        {
            gr_vec_t gens;
            gr_vec_init(gens, 0, F);
            GR_MUST_SUCCEED(gr_gens(gens, F));
            if (gens->length != 8)
            {
                flint_printf("FAIL: number of flat generators\n");
                flint_abort();
            }
            /* x = pi + e, y = pi + 2e */
            GR_MUST_SUCCEED(gr_add(x, gr_vec_entry_ptr(gens, 0, F), gr_vec_entry_ptr(gens, 1, F), F));
            GR_MUST_SUCCEED(gr_add(y, x, gr_vec_entry_ptr(gens, 1, F), F));
            gr_vec_clear(gens, F);
        }
        GR_MUST_SUCCEED(gr_div(x, x, y, F));   /* (pi + e) / (pi + 2e) */
        {
            acb_t e;
            acb_init(e);
            acb_const_pi(z, 64);
            acb_one(e); acb_exp(e, e, 64);
            acb_add(w, z, e, 64);
            acb_addmul_ui(w, e, 1, 64);        /* pi + 2e */
            acb_add(z, z, e, 64);              /* pi + e */
            acb_div(z, z, w, 64);
            acb_clear(e);
        }
        check_acb(T, x, F, z, "(pi + e) / (pi + 2e)");

        gr_heap_clear(x, F);
        gr_heap_clear(y, F);
        gr_ctx_clear(F);
        gr_tower_clear(T);
        acb_clear(z);
        acb_clear(w);
    }

    /* mixed towers: rebasing an algebraic tower, then adjoining more
       algebraic steps over the new base */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T, T2;
        gr_ctx_t K, F;
        qqbar_t a;
        acb_t z, w, u;
        gr_ptr x;
        slong k, n1 = n_randint(state, 3), n2 = n_randint(state, 2);

        qqbar_init(a);
        acb_init(z);
        acb_init(w);
        acb_init(u);
        gr_tower_init(T, QQ);

        for (k = 0; k < n1; k++)
        {
            qqbar_randtest(a, state, 1 + n_randint(state, 3), 6);
            GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, a, NULL));
        }

        /* enclosure of the top generator before rebasing */
        if (n1 > 0)
            GR_MUST_SUCCEED(gr_tower_step_get_acb(u, T, n1, 64));

        /* adjoin exp(x) or log(x) for a random top element x, or pi */
        {
            gr_ctx_struct * top = gr_tower_field(T);
            gr_ptr e;
            fmpz_mpoly_q_t fe;
            fmpz_mpoly_ctx_struct * mctx;
            int which = n_randint(state, 3);

            /* nested elements do not survive the change of base field,
               so the argument is converted to flat form first */
            GR_TMP_INIT(e, top);
            GR_MUST_SUCCEED(gr_randtest(e, state, top));
            GR_MUST_SUCCEED(gr_tower_get_acb(w, e, 64, T));
            mctx = gr_tower_flat(T)->mctx;
            fmpz_mpoly_q_init(fe, mctx);
            GR_MUST_SUCCEED(gr_tower_flat_set_nested_at(fe, e, T->length, gr_tower_flat(T)));
            GR_TMP_CLEAR(e, top);

            if (which == 0)
            {
                GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
                acb_const_pi(w, 64);
            }
            else if (which == 1)
            {
                int status = gr_tower_adjoin_exp_flat(T, fe, mctx, NULL);
                if (status == GR_DOMAIN)     /* exp(0) */
                {
                    GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
                    acb_const_pi(w, 64);
                }
                else
                {
                    GR_MUST_SUCCEED(status);
                    acb_exp(w, w, 64);
                }
            }
            else
            {
                int status = gr_tower_adjoin_log_flat(T, fe, mctx, NULL);
                if (status == GR_DOMAIN)     /* log(0) or log(1) */
                {
                    GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
                    acb_const_pi(w, 64);
                }
                else
                {
                    GR_MUST_SUCCEED(status);
                    acb_log(w, w, 64);
                }
            }
            fmpz_mpoly_q_clear(fe, mctx);
        }

        /* the nested chain was rebuilt: same degrees, same enclosures */
        if (T->length != n1 || T->num_trans != 1)
        {
            flint_printf("FAIL: structure after rebase\n");
            flint_abort();
        }
        if (n1 > 0)
        {
            GR_MUST_SUCCEED(gr_tower_step_get_acb(z, T, n1, 64));
            if (!acb_overlaps(z, u))
            {
                flint_printf("FAIL: enclosure after rebase\n");
                flint_abort();
            }
        }
        GR_MUST_SUCCEED(gr_tower_trans_get_acb(z, T, 1, 64));
        if (!acb_overlaps(z, w))
        {
            flint_printf("FAIL: transcendental generator enclosure\n");
            gr_tower_print(T);
            acb_printn(z, 20, 0); flint_printf("\n");
            acb_printn(w, 20, 0); flint_printf("\n");
            flint_abort();
        }

        /* more algebraic steps over Q(t) */
        for (k = 0; k < n2; k++)
        {
            gr_ctx_struct * top = gr_tower_field(T);
            gr_ptr e;
            GR_TMP_INIT(e, top);
            GR_MUST_SUCCEED(gr_randtest(e, state, top));
            GR_IGNORE(gr_tower_adjoin_root_ui(T, e, 2 + n_randint(state, 2), NULL));
            GR_TMP_CLEAR(e, top);
        }

        gr_ctx_init_tower_field(K, T);
        gr_ctx_init_tower_field_flat(F, T);

        gr_test_ring(F, 3 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);

        /* nested vs flat: arithmetic and evaluation */
        {
            gr_ptr p, q, r, fp, fq, fr, back;

            GR_TMP_INIT3(p, q, r, K);
            GR_TMP_INIT3(fp, fq, fr, F);
            GR_TMP_INIT(back, K);

            GR_MUST_SUCCEED(gr_randtest(p, state, K));
            GR_MUST_SUCCEED(gr_randtest(q, state, K));
            GR_MUST_SUCCEED(gr_set_other(fp, p, K, F));
            GR_MUST_SUCCEED(gr_set_other(fq, q, K, F));

            GR_MUST_SUCCEED(gr_mul(r, p, q, K));
            GR_MUST_SUCCEED(gr_add(r, r, p, K));
            GR_MUST_SUCCEED(gr_mul(fr, fp, fq, F));
            GR_MUST_SUCCEED(gr_add(fr, fr, fp, F));
            GR_MUST_SUCCEED(gr_tower_flat_get_nested(back, fr, F));

            if (gr_equal(back, r, K) != T_TRUE)
            {
                flint_printf("FAIL: flat vs nested (mixed tower)\n");
                gr_tower_print(T);
                gr_println(r, K); gr_println(fr, F); gr_println(back, K);
                flint_abort();
            }

            GR_MUST_SUCCEED(gr_tower_get_acb(z, r, 64, T));
            GR_MUST_SUCCEED(gr_tower_field_flat_get_acb(w, fr, 64, F));
            if (!acb_overlaps(z, w))
            {
                flint_printf("FAIL: nested vs flat evaluation\n");
                flint_abort();
            }

            if (gr_is_zero(r, K) == T_FALSE)
            {
                GR_MUST_SUCCEED(gr_div(r, q, r, K));
                GR_MUST_SUCCEED(gr_div(fr, fq, fr, F));
                GR_MUST_SUCCEED(gr_tower_flat_get_nested(back, fr, F));
                if (gr_equal(back, r, K) != T_TRUE)
                {
                    flint_printf("FAIL: flat vs nested division (mixed tower)\n");
                    flint_abort();
                }
            }

            GR_TMP_CLEAR3(p, q, r, K);
            GR_TMP_CLEAR3(fp, fq, fr, F);
            GR_TMP_CLEAR(back, K);
        }

        /* copying preserves the transcendental generators */
        gr_tower_init(T2, QQ);
        gr_tower_set(T2, T);
        if (T2->num_trans != T->num_trans || T2->length != T->length)
        {
            flint_printf("FAIL: copy\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_tower_trans_get_acb(z, T2, 1, 64));
        GR_MUST_SUCCEED(gr_tower_trans_get_acb(w, T, 1, 64));
        if (!acb_overlaps(z, w))
        {
            flint_printf("FAIL: copy enclosure\n");
            flint_abort();
        }
        x = gr_heap_init(gr_tower_field(T2));
        GR_MUST_SUCCEED(gr_randtest(x, state, gr_tower_field(T2)));
        GR_MUST_SUCCEED(gr_tower_get_acb(z, x, 64, T2));
        gr_heap_clear(x, gr_tower_field(T2));
        gr_tower_clear(T2);

        gr_ctx_clear(F);
        gr_ctx_clear(K);
        gr_tower_clear(T);
        qqbar_clear(a);
        acb_clear(z);
        acb_clear(w);
        acb_clear(u);
    }

    /* one pi per tower */
    {
        gr_tower_t T;
        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
        if (gr_tower_adjoin_pi(T, "p2") != GR_DOMAIN || T->num_gens != 1)
        {
            flint_printf("FAIL: a second pi adjoined\n");
            flint_abort();
        }
        gr_tower_clear(T);
    }

    /* free generators: the rational function field Q(x, y) */
    {
        gr_tower_t T;
        gr_ctx_t K, KF;
        gr_ptr a, b, c, t, f;
        acb_t z;

        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_free(T, "x"));
        GR_MUST_SUCCEED(gr_tower_adjoin_free(T, "y"));
        gr_ctx_init_tower_field(K, T);
        a = gr_heap_init(K); b = gr_heap_init(K); c = gr_heap_init(K); t = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_tower_gen_get(a, T, 0));
        GR_MUST_SUCCEED(gr_tower_gen_get(b, T, 1));
        /* (x + y)^2 - x^2 - 2 x y - y^2 = 0, x - y != 0, with 1/(x - y) */
        GR_MUST_SUCCEED(gr_add(c, a, b, K));
        GR_MUST_SUCCEED(gr_sqr(c, c, K));
        GR_MUST_SUCCEED(gr_submul(c, a, a, K));
        GR_MUST_SUCCEED(gr_submul(c, b, b, K));
        GR_MUST_SUCCEED(gr_mul(t, a, b, K));
        GR_MUST_SUCCEED(gr_mul_ui(t, t, 2, K));
        GR_MUST_SUCCEED(gr_sub(c, c, t, K));
        if (gr_is_zero(c, K) != T_TRUE)
        {
            flint_printf("FAIL: free generators, zero test\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_sub(c, a, b, K));
        if (gr_is_zero(c, K) != T_FALSE || gr_inv(t, c, K) != GR_SUCCESS)
        {
            flint_printf("FAIL: free generators, nonzero test or inverse\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_mul(t, t, c, K));
        if (gr_is_one(t, K) != T_TRUE)
        {
            flint_printf("FAIL: free generators, (x - y) / (x - y)\n");
            flint_abort();
        }
        /* no numerical value */
        acb_init(z);
        if (gr_tower_get_acb(z, c, 64, T) == GR_SUCCESS)
        {
            flint_printf("FAIL: free generators have no enclosure\n");
            flint_abort();
        }
        acb_clear(z);
        /* the flat field */
        gr_ctx_init_tower_field_flat(KF, T);
        f = gr_heap_init(KF);
        GR_MUST_SUCCEED(gr_tower_flat_set_nested(f, c, KF));
        if (gr_is_zero(f, KF) != T_FALSE)
        {
            flint_printf("FAIL: free generators, flat zero test\n");
            flint_abort();
        }
        gr_heap_clear(f, KF);
        gr_ctx_clear(KF);
        gr_heap_clear(a, K); gr_heap_clear(b, K); gr_heap_clear(c, K); gr_heap_clear(t, K);
        gr_ctx_clear(K);
        gr_tower_clear(T);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
