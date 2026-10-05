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
#include "fmpq_poly.h"
#include "acb.h"
#include "qqbar.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_tower.h"

static void
_adjoin_sqrt_ui(gr_tower_t T, ulong n, const char * name)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_ptr x;
    GR_TMP_INIT(x, top);
    GR_MUST_SUCCEED(gr_set_ui(x, n, top));
    GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, name));
    GR_TMP_CLEAR(x, top);
}

TEST_FUNCTION_START(gr_tower_map, state)
{
    gr_ctx_t QQ;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    /* Express: sqrt(6), sqrt(2) + sqrt(3) and (sqrt(2)+sqrt(3))^3 / 7 in Q(sqrt2, sqrt3);
       sqrt(5) is not in the field. */
    {
        gr_tower_t T;
        gr_ctx_struct * top;
        qqbar_t a, b, c;
        gr_ptr x, y;
        int status;

        gr_tower_init(T, QQ);
        _adjoin_sqrt_ui(T, 2, "s2");
        _adjoin_sqrt_ui(T, 3, "s3");
        top = gr_tower_field(T);

        qqbar_init(a);
        qqbar_init(b);
        qqbar_init(c);
        GR_TMP_INIT2(x, y, top);

        qqbar_set_ui(a, 2); qqbar_sqrt(a, a);
        qqbar_set_ui(b, 3); qqbar_sqrt(b, b);
        qqbar_add(c, a, b);
        qqbar_pow_ui(c, c, 3);
        qqbar_div_ui(c, c, 7);

        status = gr_tower_express_qqbar(x, c, T);
        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: express (s2+s3)^3/7, status %d\n", status);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_tower_get_qqbar(a, x, T));
        if (!qqbar_equal(a, c))
        {
            flint_printf("FAIL: express round trip\n");
            qqbar_print(a); flint_printf("\n");
            qqbar_print(c); flint_printf("\n");
            flint_abort();
        }

        qqbar_set_ui(c, 6); qqbar_sqrt(c, c);
        GR_MUST_SUCCEED(gr_tower_express_qqbar(x, c, T));
        /* x must equal s2 * s3 */
        {
            gr_vec_t gens;
            gr_vec_init(gens, 0, top);
            GR_MUST_SUCCEED(gr_gens_recursive(gens, top));
            GR_MUST_SUCCEED(gr_mul(y, gr_vec_entry_ptr(gens, 0, top), gr_vec_entry_ptr(gens, 1, top), top));
            if (gr_tower_equal(x, y, T) != T_TRUE)
            {
                flint_printf("FAIL: sqrt(6) != s2 * s3\n");
                gr_println(x, top);
                flint_abort();
            }
            /* and -sqrt(6) must map to the negative */
            qqbar_neg(c, c);
            GR_MUST_SUCCEED(gr_tower_express_qqbar(x, c, T));
            GR_MUST_SUCCEED(gr_neg(y, y, top));
            if (gr_tower_equal(x, y, T) != T_TRUE)
            {
                flint_printf("FAIL: -sqrt(6)\n");
                flint_abort();
            }
            gr_vec_clear(gens, top);
        }

        qqbar_set_ui(c, 5); qqbar_sqrt(c, c);
        status = gr_tower_express_qqbar(x, c, T);
        if (status == GR_SUCCESS)
        {
            flint_printf("FAIL: sqrt(5) expressed in Q(sqrt2, sqrt3)\n");
            flint_abort();
        }

        qqbar_set_ui(c, 2); qqbar_root_ui(c, c, 3);
        status = gr_tower_express_qqbar(x, c, T);
        if (status != GR_DOMAIN)
        {
            flint_printf("FAIL: cbrt(2) should be rejected by the degree test\n");
            flint_abort();
        }

        qqbar_clear(a);
        qqbar_clear(b);
        qqbar_clear(c);
        GR_TMP_CLEAR2(x, y, top);
        gr_tower_clear(T);
    }

    /* Merge: A = Q(s2, s3), B = Q(s6, s5). U = Q(s2, s3, s5) with s6 expressed. */
    {
        gr_tower_t A, B, U;
        gr_tower_map_t mapA, mapB;
        gr_ctx_struct * top;
        gr_vec_t gA, gB;
        gr_ptr x, y;

        gr_tower_init(A, QQ);
        _adjoin_sqrt_ui(A, 2, "s2");
        _adjoin_sqrt_ui(A, 3, "s3");
        gr_tower_init(B, QQ);
        _adjoin_sqrt_ui(B, 6, "s6");
        _adjoin_sqrt_ui(B, 5, "s5");

        gr_tower_init(U, QQ);
        GR_MUST_SUCCEED(gr_tower_merge(U, mapA, mapB, A, B, GR_TOWER_MERGE_EXPRESS));

        if (gr_tower_length(U) != 3 || gr_tower_degree(U) != 8)
        {
            flint_printf("FAIL: merged tower\n");
            gr_tower_print(U);
            flint_abort();
        }

        top = gr_tower_field(U);
        gr_vec_init(gA, 0, gr_tower_field(A));
        gr_vec_init(gB, 0, gr_tower_field(B));
        GR_MUST_SUCCEED(gr_gens_recursive(gA, gr_tower_field(A)));
        GR_MUST_SUCCEED(gr_gens_recursive(gB, gr_tower_field(B)));
        GR_TMP_INIT2(x, y, top);

        /* mapA(s2) * mapA(s3) == mapB(s6) */
        GR_MUST_SUCCEED(gr_tower_map_apply(x, gr_vec_entry_ptr(gA, 0, gr_tower_field(A)), mapA));
        GR_MUST_SUCCEED(gr_tower_map_apply(y, gr_vec_entry_ptr(gA, 1, gr_tower_field(A)), mapA));
        GR_MUST_SUCCEED(gr_mul(x, x, y, top));
        GR_MUST_SUCCEED(gr_tower_map_apply(y, gr_vec_entry_ptr(gB, 0, gr_tower_field(B)), mapB));
        if (gr_tower_equal(x, y, U) != T_TRUE)
        {
            flint_printf("FAIL: s2 s3 != s6 in the merged tower\n");
            gr_println(x, top); gr_println(y, top);
            flint_abort();
        }

        /* mapB(s5)^2 == 5 and mapB(s5) is the new generator */
        GR_MUST_SUCCEED(gr_tower_map_apply(y, gr_vec_entry_ptr(gB, 1, gr_tower_field(B)), mapB));
        GR_MUST_SUCCEED(gr_sqr(x, y, top));
        if (gr_sub_ui(x, x, 5, top) != GR_SUCCESS || gr_tower_is_zero(x, U) != T_TRUE)
        {
            flint_printf("FAIL: s5^2 != 5\n");
            flint_abort();
        }

        /* map applied to an expression: (s6 + s5)^2 = 11 + 2 s30 */
        {
            gr_ctx_struct * Btop = gr_tower_field(B);
            gr_ptr e;
            GR_TMP_INIT(e, Btop);
            GR_MUST_SUCCEED(gr_add(e, gr_vec_entry_ptr(gB, 0, Btop), gr_vec_entry_ptr(gB, 1, Btop), Btop));
            GR_MUST_SUCCEED(gr_sqr(e, e, Btop));
            GR_MUST_SUCCEED(gr_tower_map_apply(x, e, mapB));
            /* compare with the same computation in U */
            GR_MUST_SUCCEED(gr_tower_map_apply(y, gr_vec_entry_ptr(gB, 0, Btop), mapB));
            {
                gr_ptr w;
                GR_TMP_INIT(w, top);
                GR_MUST_SUCCEED(gr_tower_map_apply(w, gr_vec_entry_ptr(gB, 1, Btop), mapB));
                GR_MUST_SUCCEED(gr_add(y, y, w, top));
                GR_TMP_CLEAR(w, top);
            }
            GR_MUST_SUCCEED(gr_sqr(y, y, top));
            if (gr_tower_equal(x, y, U) != T_TRUE)
            {
                flint_printf("FAIL: map is not a homomorphism\n");
                flint_abort();
            }
            GR_TMP_CLEAR(e, Btop);
        }

        GR_TMP_CLEAR2(x, y, top);
        gr_vec_clear(gA, gr_tower_field(A));
        gr_vec_clear(gB, gr_tower_field(B));
        gr_tower_map_clear(mapA);
        gr_tower_map_clear(mapB);
        gr_tower_clear(U);
        gr_tower_clear(A);
        gr_tower_clear(B);
    }

    /* Merge with transcendental generators. A = Q(s2)(t = exp(s2)),
       B = Q(u = exp(1), s2', s8) where s2' has the same definition as s2.
       U = A + B: exp(1) is adjoined, s2' is identified with s2 (by
       definition id), s8 = 2 s2 is expressed. */
    {
        gr_tower_t A, B, U;
        gr_tower_map_t mapA, mapB;
        gr_ctx_struct * top;
        gr_ptr x, y;
        acb_t z, w;

        acb_init(z);
        acb_init(w);

        gr_tower_init(A, QQ);
        _adjoin_sqrt_ui(A, 2, "s2");
        gr_tower_gen_set_def_id(A, GR_TOWER_STEP(A, 0)->def_order, 1001);
        {
            /* nested elements do not survive adjoining a transcendental
               generator, so the argument is passed in flat form */
            fmpz_mpoly_q_t g;
            fmpz_mpoly_ctx_struct * mctx = gr_tower_flat(A)->mctx;
            fmpz_mpoly_q_init(g, mctx);
            fmpz_mpoly_q_gen(g, GR_TOWER_FLAT_VAR(gr_tower_flat(A), 1), mctx);
            GR_MUST_SUCCEED(gr_tower_adjoin_exp_flat(A, g, mctx, "t"));
            fmpz_mpoly_q_clear(g, mctx);
            gr_tower_gen_set_def_id(A, GR_TOWER_TRANS(A, 0)->def_order, 1002);
        }

        gr_tower_init(B, QQ);
        {
            gr_ptr g = gr_heap_init(QQ);
            GR_MUST_SUCCEED(gr_one(g, QQ));
            GR_MUST_SUCCEED(gr_tower_adjoin_exp(B, g, "u"));
            gr_heap_clear(g, QQ);
            gr_tower_gen_set_def_id(B, GR_TOWER_TRANS(B, 0)->def_order, 1003);
        }
        _adjoin_sqrt_ui(B, 2, "s2");
        gr_tower_gen_set_def_id(B, GR_TOWER_STEP(B, 0)->def_order, 1001);
        _adjoin_sqrt_ui(B, 8, "s8");
        gr_tower_gen_set_def_id(B, GR_TOWER_STEP(B, 1)->def_order, 1004);

        gr_tower_init(U, QQ);
        GR_MUST_SUCCEED(gr_tower_merge(U, mapA, mapB, A, B, GR_TOWER_MERGE_EXPRESS));

        if (gr_tower_length(U) != 1 || U->num_trans != 2 || gr_tower_num_gens(U) != 3)
        {
            flint_printf("FAIL: merged mixed tower\n");
            gr_tower_print(U);
            flint_abort();
        }

        top = gr_tower_field(U);
        GR_TMP_INIT2(x, y, top);

        /* mapB(s8) == 2 mapA(s2) */
        {
            gr_ptr g = gr_heap_init(gr_tower_field(B));
            GR_MUST_SUCCEED(gr_gen(g, gr_tower_field(B)));
            GR_MUST_SUCCEED(gr_tower_map_apply(x, g, mapB));
            gr_heap_clear(g, gr_tower_field(B));
        }
        {
            gr_ptr g = gr_heap_init(gr_tower_field(A));
            GR_MUST_SUCCEED(gr_gen(g, gr_tower_field(A)));
            GR_MUST_SUCCEED(gr_tower_map_apply(y, g, mapA));
            gr_heap_clear(g, gr_tower_field(A));
        }
        GR_MUST_SUCCEED(gr_mul_ui(y, y, 2, top));
        if (gr_tower_equal(x, y, U) != T_TRUE)
        {
            flint_printf("FAIL: s8 != 2 s2 in the merged mixed tower\n");
            gr_println(x, top); gr_println(y, top);
            flint_abort();
        }

        /* the images of the transcendental generators evaluate correctly:
           mapA(t) = exp(sqrt(2)), mapB(u) = e */
        {
            fmpz_mpoly_q_t img;
            gr_tower_map_sync(mapA);
            fmpz_mpoly_q_init(img, mapA->mctx);
            fmpz_mpoly_q_set(img, gr_tower_map_image(mapA, 1), mapA->mctx);
            GR_MUST_SUCCEED(gr_tower_flat_get_acb(z, img, 64, gr_tower_flat(U)));
            acb_set_ui(w, 2); acb_sqrt(w, w, 64); acb_exp(w, w, 64);
            if (!acb_overlaps(z, w))
            {
                flint_printf("FAIL: image of exp(s2)\n");
                flint_abort();
            }
            fmpz_mpoly_q_clear(img, mapA->mctx);

            gr_tower_map_sync(mapB);
            fmpz_mpoly_q_init(img, mapB->mctx);
            fmpz_mpoly_q_set(img, gr_tower_map_image(mapB, 0), mapB->mctx);
            GR_MUST_SUCCEED(gr_tower_flat_get_acb(z, img, 64, gr_tower_flat(U)));
            acb_one(w); acb_exp(w, w, 64);
            if (!acb_overlaps(z, w))
            {
                flint_printf("FAIL: image of exp(1)\n");
                flint_abort();
            }
            fmpz_mpoly_q_clear(img, mapB->mctx);
        }

        /* prefix copy in definition order: the first two generators of B */
        {
            gr_tower_t P;
            gr_tower_init(P, QQ);
            gr_tower_set_prefix(P, B, 2);
            if (P->num_trans != 1 || P->length != 1 || gr_tower_num_gens(P) != 2 ||
                gr_tower_prefix_length(B, 2) != 1 || gr_tower_prefix_num_trans(B, 2) != 1)
            {
                flint_printf("FAIL: prefix\n");
                gr_tower_print(P);
                flint_abort();
            }
            gr_tower_clear(P);
        }

        GR_TMP_CLEAR2(x, y, top);
        gr_tower_map_clear(mapA);
        gr_tower_map_clear(mapB);
        gr_tower_clear(U);
        gr_tower_clear(A);
        gr_tower_clear(B);
        acb_clear(z);
        acb_clear(w);
    }

    /* Eliminate: after refinement of Q(s2)(s8), the second step has degree 1 */
    {
        gr_tower_t T, U;
        gr_tower_map_t map;
        gr_ctx_t K;
        gr_vec_t gens;
        gr_ptr t, u;
        qqbar_t q;

        gr_tower_init(T, QQ);
        _adjoin_sqrt_ui(T, 2, "s2");
        qqbar_init(q);
        qqbar_set_ui(q, 8);
        qqbar_sqrt(q, q);
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, q, "s8"));

        gr_ctx_init_tower_field(K, T);
        gr_vec_init(gens, 0, K);
        GR_MUST_SUCCEED(gr_gens(gens, K));
        GR_TMP_INIT(t, K);
        GR_MUST_SUCCEED(gr_mul_si(t, gr_vec_entry_ptr(gens, 0, K), 2, K));
        GR_MUST_SUCCEED(gr_sub(t, gr_vec_entry_ptr(gens, 1, K), t, K));
        if (gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: refinement\n");
            flint_abort();
        }

        gr_tower_init(U, QQ);
        GR_MUST_SUCCEED(gr_tower_eliminate(U, map, T));

        if (gr_tower_length(U) != 1 || gr_tower_degree(U) != 2)
        {
            flint_printf("FAIL: eliminate\n");
            gr_tower_print(U);
            flint_abort();
        }

        /* image of s8 is 2 s2 */
        GR_TMP_INIT(u, gr_tower_field(U));
        GR_MUST_SUCCEED(gr_tower_map_apply(u, gr_vec_entry_ptr(gens, 1, K), map));
        {
            gr_ptr g;
            GR_TMP_INIT(g, gr_tower_field(U));
            GR_MUST_SUCCEED(gr_gen(g, gr_tower_field(U)));
            GR_MUST_SUCCEED(gr_mul_si(g, g, 2, gr_tower_field(U)));
            if (gr_tower_equal(u, g, U) != T_TRUE)
            {
                flint_printf("FAIL: image of s8\n");
                gr_println(u, gr_tower_field(U));
                flint_abort();
            }
            GR_TMP_CLEAR(g, gr_tower_field(U));
        }

        GR_TMP_CLEAR(u, gr_tower_field(U));
        GR_TMP_CLEAR(t, K);
        gr_vec_clear(gens, K);
        gr_ctx_clear(K);
        gr_tower_map_clear(map);
        gr_tower_clear(U);
        gr_tower_clear(T);
        qqbar_clear(q);
    }

    /* Random: two towers from random qqbars, merge, and check that
       arithmetic commutes with the maps and with qqbar arithmetic. */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_tower_t A, B, U;
        gr_tower_map_t mapA, mapB;
        qqbar_t a, b, c, d;
        gr_ctx_struct * Atop, * Btop, * Utop;
        gr_ptr x, y, z;
        int flags = n_randint(state, 2) ? GR_TOWER_MERGE_EXPRESS : 0;

        qqbar_init(a); qqbar_init(b); qqbar_init(c); qqbar_init(d);
        qqbar_randtest(a, state, 1 + n_randint(state, 3), 6);
        qqbar_randtest(b, state, 1 + n_randint(state, 3), 6);

        gr_tower_init(A, QQ);
        gr_tower_init(B, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(A, a, "a"));
        if (n_randint(state, 2))
            GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(A, b, "b"));
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(B, b, "b"));
        if (n_randint(state, 2))
            GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(B, a, "a"));

        gr_tower_init(U, QQ);
        GR_MUST_SUCCEED(gr_tower_merge(U, mapA, mapB, A, B, flags));

        Atop = gr_tower_field(A);
        Btop = gr_tower_field(B);
        Utop = gr_tower_field(U);
        GR_TMP_INIT3(x, y, z, Utop);

        /* x = image of a (first generator of A), y = image of b (first generator of B) */
        {
            gr_vec_t gv;
            gr_vec_init(gv, 0, Atop);
            GR_MUST_SUCCEED(gr_gens_recursive(gv, Atop));
            GR_MUST_SUCCEED(gr_tower_map_apply(x, gr_vec_entry_ptr(gv, 0, Atop), mapA));
            gr_vec_clear(gv, Atop);
            gr_vec_init(gv, 0, Btop);
            GR_MUST_SUCCEED(gr_gens_recursive(gv, Btop));
            GR_MUST_SUCCEED(gr_tower_map_apply(y, gr_vec_entry_ptr(gv, 0, Btop), mapB));
            gr_vec_clear(gv, Btop);
        }

        /* z = x*y - x + 2y, compare with qqbar */
        GR_MUST_SUCCEED(gr_mul(z, x, y, Utop));
        GR_MUST_SUCCEED(gr_sub(z, z, x, Utop));
        GR_MUST_SUCCEED(gr_mul_si(y, y, 2, Utop));
        GR_MUST_SUCCEED(gr_add(z, z, y, Utop));

        qqbar_mul(c, a, b);
        qqbar_sub(c, c, a);
        qqbar_mul_si(d, b, 2);
        qqbar_add(c, c, d);

        GR_MUST_SUCCEED(gr_tower_get_qqbar(d, z, U));
        if (!qqbar_equal(c, d))
        {
            flint_printf("FAIL: merge arithmetic\n");
            gr_tower_print(A); gr_tower_print(B); gr_tower_print(U);
            flint_abort();
        }

        GR_TMP_CLEAR3(x, y, z, Utop);
        gr_tower_map_clear(mapA);
        gr_tower_map_clear(mapB);
        gr_tower_clear(U);
        gr_tower_clear(A);
        gr_tower_clear(B);
        qqbar_clear(a); qqbar_clear(b); qqbar_clear(c); qqbar_clear(d);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
