/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include "gr_tower/impl.h"

/*
    The real trigonometric generators tan(u) and atan(u) and Richardson's
    algorithm on their angles: relations with pi (Machin's formula),
    between tangents of multiple angles (in both definition orders),
    between atan and tan, and, with i present, with the exponentials of
    the complex field; random identities of the real lazy field
    (addition formulas at random rational arguments) are decided, and
    perturbed ones refuted.
*/

#define CHECK_TRUTH(x, v, what) \
    do { truth_t _t = (x); if (_t != (v)) { flint_printf("FAIL: %s (%d)\n", what, _t); flint_abort(); } } while (0)

/* the rational c as a flat element in T's context */
static void
_flat_fmpq(fmpz_mpoly_q_t res, slong p, ulong q, gr_tower_t T)
{
    fmpq_t c;
    fmpq_init(c);
    fmpq_set_si(c, p, q);
    gr_tower_flat_ensure(&T->flat);
    fmpz_mpoly_q_set_fmpq(res, c, T->flat.mctx);
    fmpq_clear(c);
}

static int
_adjoin_const(gr_tower_t T, int kind, slong p, ulong q)
{
    fmpz_mpoly_q_t u;
    int status;
    gr_tower_flat_ensure(&T->flat);
    fmpz_mpoly_q_init(u, T->flat.mctx);
    _flat_fmpq(u, p, q, T);
    status = (kind == GR_TOWER_TAN) ? gr_tower_adjoin_tan_flat(T, u, T->flat.mctx, NULL)
                                    : gr_tower_adjoin_atan_flat(T, u, T->flat.mctx, NULL);
    fmpz_mpoly_q_clear(u, T->flat.mctx);
    return status;
}

TEST_FUNCTION_START(gr_tower_trig, state)
{
    gr_ctx_t QQ;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    /* tower level: Machin's formula 4 atan(1/5) - atan(1/239) = pi/4, and
       atan(1/2) + atan(1/3) = pi/4 */
    {
        gr_tower_t T;
        gr_ctx_t K;
        gr_ptr a, b, c, d, pi, x;

        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
        GR_MUST_SUCCEED(_adjoin_const(T, GR_TOWER_ATAN, 1, 5));
        GR_MUST_SUCCEED(_adjoin_const(T, GR_TOWER_ATAN, 1, 239));
        GR_MUST_SUCCEED(_adjoin_const(T, GR_TOWER_ATAN, 1, 2));
        GR_MUST_SUCCEED(_adjoin_const(T, GR_TOWER_ATAN, 1, 3));

        gr_ctx_init_tower_field(K, T);
        pi = gr_heap_init(K); a = gr_heap_init(K); b = gr_heap_init(K);
        c = gr_heap_init(K); d = gr_heap_init(K); x = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_tower_gen_get(pi, T, 0));
        GR_MUST_SUCCEED(gr_tower_gen_get(a, T, 1));
        GR_MUST_SUCCEED(gr_tower_gen_get(b, T, 2));
        GR_MUST_SUCCEED(gr_tower_gen_get(c, T, 3));
        GR_MUST_SUCCEED(gr_tower_gen_get(d, T, 4));

        GR_MUST_SUCCEED(gr_mul_si(x, a, 4, K));
        GR_MUST_SUCCEED(gr_sub(x, x, b, K));
        GR_MUST_SUCCEED(gr_div_si(pi, pi, 4, K));
        GR_MUST_SUCCEED(gr_sub(x, x, pi, K));
        CHECK_TRUTH(gr_is_zero(x, K), T_TRUE, "Machin");

        GR_MUST_SUCCEED(gr_add(x, c, d, K));
        GR_MUST_SUCCEED(gr_sub(x, x, pi, K));
        CHECK_TRUTH(gr_is_zero(x, K), T_TRUE, "atan(1/2) + atan(1/3)");
        GR_MUST_SUCCEED(gr_add_si(x, x, 1, K));
        CHECK_TRUTH(gr_is_zero(x, K), T_FALSE, "atan(1/2) + atan(1/3) + 1");

        gr_heap_clear(pi, K); gr_heap_clear(a, K); gr_heap_clear(b, K);
        gr_heap_clear(c, K); gr_heap_clear(d, K); gr_heap_clear(x, K);
        gr_ctx_clear(K);
        gr_tower_clear(T);
    }

    /* tan(2) = 2 tan(1) / (1 - tan(1)^2), with the generators in both
       orders; atan(tan(1)) = 1 */
    {
        slong order;
        for (order = 0; order < 2; order++)
        {
            gr_tower_t T;
            gr_ctx_t K;
            gr_ptr t1, t2, x, y;

            gr_tower_init(T, QQ);
            GR_MUST_SUCCEED(_adjoin_const(T, GR_TOWER_TAN, order ? 2 : 1, 1));
            GR_MUST_SUCCEED(_adjoin_const(T, GR_TOWER_TAN, order ? 1 : 2, 1));

            gr_ctx_init_tower_field(K, T);
            t1 = gr_heap_init(K); t2 = gr_heap_init(K); x = gr_heap_init(K); y = gr_heap_init(K);
            GR_MUST_SUCCEED(gr_tower_gen_get(t1, T, order ? 1 : 0));
            GR_MUST_SUCCEED(gr_tower_gen_get(t2, T, order ? 0 : 1));

            GR_MUST_SUCCEED(gr_sqr(x, t1, K));
            GR_MUST_SUCCEED(gr_neg(x, x, K));
            GR_MUST_SUCCEED(gr_add_si(x, x, 1, K));
            GR_MUST_SUCCEED(gr_mul_si(y, t1, 2, K));
            GR_MUST_SUCCEED(gr_div(x, y, x, K));
            CHECK_TRUTH(gr_equal(x, t2, K), T_TRUE, "tan(2) double angle");
            CHECK_TRUTH(gr_equal(x, t1, K), T_FALSE, "tan(2) != tan(1)");

            gr_heap_clear(t1, K); gr_heap_clear(t2, K); gr_heap_clear(x, K); gr_heap_clear(y, K);
            gr_ctx_clear(K);
            gr_tower_clear(T);
        }
    }

    /* the real lazy field */
    {
        gr_ctx_t C, R;
        gr_ptr x, y, z, w, pi;
        char * s;

        gr_ctx_init_tower_lazy(C, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_ctx_init_tower_lazy_view(R, C, GR_TOWER_LAZY_REAL);
        x = gr_heap_init(R); y = gr_heap_init(R); z = gr_heap_init(R);
        w = gr_heap_init(R); pi = gr_heap_init(R);

        /* cos(1) in terms of tan(1/2), without complex generators */
        GR_MUST_SUCCEED(gr_set_si(x, 1, R));
        GR_MUST_SUCCEED(gr_cos(y, x, R));
        GR_MUST_SUCCEED(gr_get_str(&s, y, R));
        if (strstr(s, "tan(1/2)") == NULL || strstr(s, "i") != NULL)
        {
            flint_printf("FAIL: cos(1) = %s\n", s);
            flint_abort();
        }
        flint_free(s);

        /* atan(tan(1)) = 1, atan(tan(2)) = 2 - pi, tan(atan(2)) = 2 */
        GR_MUST_SUCCEED(gr_tan(y, x, R));
        GR_MUST_SUCCEED(gr_atan(y, y, R));
        CHECK_TRUTH(gr_is_one(y, R), T_TRUE, "atan(tan(1))");
        GR_MUST_SUCCEED(gr_set_si(x, 2, R));
        GR_MUST_SUCCEED(gr_tan(y, x, R));
        GR_MUST_SUCCEED(gr_atan(y, y, R));
        GR_MUST_SUCCEED(gr_pi(pi, R));
        GR_MUST_SUCCEED(gr_sub(z, x, pi, R));
        CHECK_TRUTH(gr_equal(y, z, R), T_TRUE, "atan(tan(2))");
        GR_MUST_SUCCEED(gr_atan(y, x, R));
        GR_MUST_SUCCEED(gr_tan(y, y, R));
        CHECK_TRUTH(gr_equal(y, x, R), T_TRUE, "tan(atan(2))");

        /* the same value in the complex field (exp form) */
        GR_MUST_SUCCEED(gr_set_si(x, 1, R));
        GR_MUST_SUCCEED(gr_cos(y, x, R));
        GR_MUST_SUCCEED(gr_set_si(z, 1, C));
        GR_MUST_SUCCEED(gr_cos(z, z, C));
        CHECK_TRUTH(gr_equal(y, z, C), T_TRUE, "cos(1) real vs complex");

        /* a basis change of tangents in which one new argument is an old
           one: with t1 = tan(3 - s), t2 = tan(2 - s) present (s = sqrt 2),
           tan(s/2) gets the basis tan(3 - s), tan((2 - s)/2) with t1, t2
           eliminated (all degrees 1), rather than a modulus of degree 4
           over Q(t1, t2) whose arithmetic is very slow */
        {
            gr_ptr s, a, b;
            gr_tower_struct * T;
            slong lev, d;
            GR_TMP_INIT3(s, a, b, R);
            GR_MUST_SUCCEED(gr_set_si(s, 2, R));
            GR_MUST_SUCCEED(gr_sqrt(s, s, R));
            GR_MUST_SUCCEED(gr_sub_si(a, s, 3, R));
            GR_MUST_SUCCEED(gr_neg(a, a, R));
            GR_MUST_SUCCEED(gr_tan(a, a, R));
            GR_MUST_SUCCEED(gr_sub_si(b, s, 2, R));
            GR_MUST_SUCCEED(gr_neg(b, b, R));
            GR_MUST_SUCCEED(gr_tan(b, b, R));
            GR_MUST_SUCCEED(gr_div_si(x, s, 2, R));
            GR_MUST_SUCCEED(gr_tan(x, x, R));
            T = gr_tower_lazy_get_tower(&lev, x, R);
            for (d = 0; d < T->num_gens; d++)
            {
                if (GR_TOWER_GEN(T, d)->kind == GR_TOWER_ALGEBRAIC && gr_tower_step_degree(T, GR_TOWER_GEN(T, d)->index) > 2)
                {
                    flint_printf("FAIL: tan(sqrt(2)/2) of degree > 2 over the other tangents\n");
                    gr_tower_print(T);
                    flint_abort();
                }
            }
            GR_TMP_CLEAR3(s, a, b, R);
        }

        /* random addition formulas: sin(a + b) = sin a cos b + cos a sin b
           and tan(a + b) = (tan a + tan b)/(1 - tan a tan b) at random
           rationals (plus sqrt(2) sometimes) */
        for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
        {
            gr_ptr a, b, t, u;
            fmpq_t q;
            int which = n_randint(state, 2);

            /* (a fresh context now and then: the towers of one context
               accumulate all the tangents) */
            if (iter % 5 == 4)
            {
                gr_heap_clear(x, R); gr_heap_clear(y, R); gr_heap_clear(z, R);
                gr_heap_clear(w, R); gr_heap_clear(pi, R);
                gr_ctx_clear(R);
                gr_ctx_clear(C);
                gr_ctx_init_tower_lazy(C, QQ, GR_TOWER_MERGE_EXPRESS);
                gr_ctx_init_tower_lazy_view(R, C, GR_TOWER_LAZY_REAL);
                x = gr_heap_init(R); y = gr_heap_init(R); z = gr_heap_init(R);
                w = gr_heap_init(R); pi = gr_heap_init(R);
            }

            a = gr_heap_init(R); b = gr_heap_init(R); t = gr_heap_init(R); u = gr_heap_init(R);
            fmpq_init(q);

            /* (small numerators over the denominators 1 and 2: the
               tangents of commensurable angles are algebraic over each
               other, and a common primitive angle with a large
               multiplicity would make the test slow) */
            fmpq_set_si(q, (slong) n_randint(state, 7) - 3, 1 + n_randint(state, 2));
            GR_MUST_SUCCEED(gr_set_fmpq(a, q, R));
            fmpq_set_si(q, (slong) n_randint(state, 7) - 3, 1 + n_randint(state, 2));
            GR_MUST_SUCCEED(gr_set_fmpq(b, q, R));
            if (n_randint(state, 2))
            {
                GR_MUST_SUCCEED(gr_set_ui(t, 2, R));
                GR_MUST_SUCCEED(gr_sqrt(t, t, R));
                GR_MUST_SUCCEED(gr_add(a, a, t, R));
            }

            if (gr_is_zero(a, R) != T_FALSE || gr_is_zero(b, R) != T_FALSE)
                goto next;

            GR_MUST_SUCCEED(gr_add(t, a, b, R));
            if (which == 0)
            {
                GR_MUST_SUCCEED(gr_sin(x, t, R));
                GR_MUST_SUCCEED(gr_sin(y, a, R));
                GR_MUST_SUCCEED(gr_cos(z, b, R));
                GR_MUST_SUCCEED(gr_mul(y, y, z, R));
                GR_MUST_SUCCEED(gr_cos(z, a, R));
                GR_MUST_SUCCEED(gr_sin(u, b, R));
                GR_MUST_SUCCEED(gr_mul(z, z, u, R));
                GR_MUST_SUCCEED(gr_add(y, y, z, R));
            }
            else
            {
                GR_MUST_SUCCEED(gr_tan(x, t, R));
                GR_MUST_SUCCEED(gr_tan(y, a, R));
                GR_MUST_SUCCEED(gr_tan(z, b, R));
                GR_MUST_SUCCEED(gr_mul(u, y, z, R));
                GR_MUST_SUCCEED(gr_neg(u, u, R));
                GR_MUST_SUCCEED(gr_add_si(u, u, 1, R));
                GR_MUST_SUCCEED(gr_add(y, y, z, R));
                if (gr_is_zero(u, R) != T_FALSE)
                    goto next;
                GR_MUST_SUCCEED(gr_div(y, y, u, R));
            }

            CHECK_TRUTH(gr_equal(x, y, R), T_TRUE, "addition formula");
            GR_MUST_SUCCEED(gr_set_si(w, 1, R));
            GR_MUST_SUCCEED(gr_div_si(w, w, 1000, R));
            GR_MUST_SUCCEED(gr_add(y, y, w, R));
            CHECK_TRUTH(gr_equal(x, y, R), T_FALSE, "perturbed addition formula");

        next:
            fmpq_clear(q);
            gr_heap_clear(a, R); gr_heap_clear(b, R); gr_heap_clear(t, R); gr_heap_clear(u, R);
        }

        gr_heap_clear(x, R); gr_heap_clear(y, R); gr_heap_clear(z, R);
        gr_heap_clear(w, R); gr_heap_clear(pi, R);
        gr_ctx_clear(R);
        gr_ctx_clear(C);
    }

    /* the tangent normal form of the real trigonometric constants: at
       r pi (r = p/q), sin, cos and tan are polynomials in one generator
       tan(pi/M) of the level (no roots of unity), with the exact
       identities sin^2 + cos^2 = 1, tan cos = sin, and the values of
       acb; sin(2 pi/5) through tan(pi/5), tan(pi/12) = 2 - sqrt(3) */
    {
        gr_ctx_t C, R;
        gr_ptr x, s, c, t, u;
        char * str;
        acb_t z, v;

        gr_ctx_init_tower_lazy(C, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_ctx_init_tower_lazy_view(R, C, GR_TOWER_LAZY_REAL);
        x = gr_heap_init(R); s = gr_heap_init(R); c = gr_heap_init(R);
        t = gr_heap_init(R); u = gr_heap_init(R);
        acb_init(z);
        acb_init(v);

        GR_MUST_SUCCEED(gr_pi(x, R));
        GR_MUST_SUCCEED(gr_mul_ui(x, x, 2, R));
        GR_MUST_SUCCEED(gr_div_ui(x, x, 5, R));
        GR_MUST_SUCCEED(gr_sin(s, x, R));
        GR_MUST_SUCCEED(gr_get_str(&str, s, R));
        if (strstr(str, "tan(pi/5)") == NULL || strstr(str, "exp(") != NULL)
        {
            flint_printf("FAIL: sin(2 pi/5) = %s\n", str);
            flint_abort();
        }
        flint_free(str);

        GR_MUST_SUCCEED(gr_pi(x, R));
        GR_MUST_SUCCEED(gr_div_ui(x, x, 12, R));
        GR_MUST_SUCCEED(gr_tan(t, x, R));
        GR_MUST_SUCCEED(gr_set_ui(u, 3, R));
        GR_MUST_SUCCEED(gr_sqrt(u, u, R));
        GR_MUST_SUCCEED(gr_neg(u, u, R));
        GR_MUST_SUCCEED(gr_add_ui(u, u, 2, R));
        CHECK_TRUTH(gr_equal(t, u, R), T_TRUE, "tan(pi/12) = 2 - sqrt(3)");
        GR_MUST_SUCCEED(gr_get_str(&str, t, R));
        if (strstr(str, "tan(pi/24)") != NULL)
        {
            flint_printf("FAIL: tan(pi/12) = %s\n", str);
            flint_abort();
        }
        flint_free(str);

        for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
        {
            slong q = 1 + n_randint(state, 30), p = (slong) n_randint(state, 4 * q + 1) - 2 * q;
            fmpq_t r;
            int which;

            fmpq_init(r);
            fmpq_set_si(r, p, q);
            GR_MUST_SUCCEED(gr_pi(x, R));
            GR_MUST_SUCCEED(gr_mul_fmpq(x, x, r, R));
            GR_MUST_SUCCEED(gr_sin(s, x, R));
            GR_MUST_SUCCEED(gr_cos(c, x, R));

            for (which = 0; which < 2; which++)
            {
                char * q2;
                GR_MUST_SUCCEED(gr_get_str(&str, which ? c : s, R));
                for (q2 = str; *q2; q2++)
                    if (q2[0] == 'p' && q2[1] == 'i')
                        q2[0] = q2[1] = 'P';
                if (strstr(str, "exp(") != NULL || strstr(str, "i") != NULL)
                {
                    flint_printf("FAIL: tangent normal form at %wd/%wd: %s\n", p, q, str);
                    flint_abort();
                }
                flint_free(str);

                GR_MUST_SUCCEED(gr_tower_lazy_get_acb(z, which ? c : s, 128, R));
                acb_zero(v);
                if (which)
                    arb_cos_pi_fmpq(acb_realref(v), r, 128);
                else
                    arb_sin_pi_fmpq(acb_realref(v), r, 128);
                if (!acb_overlaps(z, v))
                {
                    flint_printf("FAIL: value at %wd/%wd (%d)\n", p, q, which);
                    flint_abort();
                }
            }

            GR_MUST_SUCCEED(gr_sqr(t, s, R));
            GR_MUST_SUCCEED(gr_sqr(u, c, R));
            GR_MUST_SUCCEED(gr_add(t, t, u, R));
            CHECK_TRUTH(gr_is_one(t, R), T_TRUE, "sin^2 + cos^2 = 1");

            if (gr_is_zero(c, R) == T_FALSE)
            {
                GR_MUST_SUCCEED(gr_tan(t, x, R));
                GR_MUST_SUCCEED(gr_mul(t, t, c, R));
                CHECK_TRUTH(gr_equal(t, s, R), T_TRUE, "tan cos = sin");
            }
            else if (gr_tan(t, x, R) != GR_DOMAIN)
            {
                flint_printf("FAIL: tan at a pole\n");
                flint_abort();
            }

            fmpq_clear(r);
        }

        acb_clear(z);
        acb_clear(v);
        gr_heap_clear(x, R); gr_heap_clear(s, R); gr_heap_clear(c, R);
        gr_heap_clear(t, R); gr_heap_clear(u, R);
        gr_ctx_clear(R);
        gr_ctx_clear(C);
    }

    /* tan(atan(2)/3): the angle relation 3 A_g = A_a gives a basis whose
       tangent is the generator itself (the relation search must not
       restart on it); the degree-6 modulus is refined to the cubic
       tan(3u) = 2 by the zero test */
    {
        gr_tower_t T;
        gr_ctx_t K;
        gr_ptr u, t, n, d;

        gr_tower_init(T, QQ);
        gr_ctx_init_tower_field(K, T);
        u = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_set_ui(u, 2, K));
        GR_MUST_SUCCEED(gr_tower_adjoin_atan(T, u, "A"));
        gr_heap_clear(u, K);
        gr_ctx_clear(K);
        gr_ctx_init_tower_field(K, T);
        u = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_tower_gen_get(u, T, 0));
        GR_MUST_SUCCEED(gr_div_ui(u, u, 3, K));
        GR_MUST_SUCCEED(gr_tower_adjoin_tan(T, u, "B"));
        gr_heap_clear(u, K);
        gr_ctx_clear(K);
        /* the relation search (which used to loop forever here) finds
           the degree-6 relation tan(3 B) = 2 */
        gr_tower_flat_ensure(&T->flat);
        if (!_gr_tower_search_relations(&T->flat, 64) || T->length != 1 || gr_tower_step_degree(T, 1) != 6)
        {
            flint_printf("FAIL: tan(atan(2)/3) not found algebraic of degree 6\n");
            gr_tower_print(T);
            flint_abort();
        }
        gr_ctx_init_tower_field(K, T);
        t = gr_heap_init(K); n = gr_heap_init(K); d = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_tower_gen_get(t, T, 1));
        /* tan(3u) = (3t - t^3) / (1 - 3t^2) = 2 */
        GR_MUST_SUCCEED(gr_pow_ui(n, t, 3, K));
        GR_MUST_SUCCEED(gr_neg(n, n, K));
        GR_MUST_SUCCEED(gr_addmul_ui(n, t, 3, K));
        GR_MUST_SUCCEED(gr_sqr(d, t, K));
        GR_MUST_SUCCEED(gr_mul_ui(d, d, 3, K));
        GR_MUST_SUCCEED(gr_neg(d, d, K));
        GR_MUST_SUCCEED(gr_add_ui(d, d, 1, K));
        GR_MUST_SUCCEED(gr_mul_ui(d, d, 2, K));
        GR_MUST_SUCCEED(gr_sub(n, n, d, K));
        CHECK_TRUTH(gr_is_zero(n, K), T_TRUE, "tan(3 atan(2)/3) = 2");
        if (gr_tower_step_degree(T, 1) != 3)
        {
            flint_printf("FAIL: tan(atan(2)/3) should have degree 3 (%wd)\n", gr_tower_step_degree(T, 1));
            flint_abort();
        }
        gr_heap_clear(t, K); gr_heap_clear(n, K); gr_heap_clear(d, K);
        gr_ctx_clear(K);
        gr_tower_clear(T);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
