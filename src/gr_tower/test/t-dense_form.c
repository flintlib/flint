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
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/* The dense forms of elements of number fields of one generator (and
   the rational fast paths): random sequences of operations in a lazy
   field with them and in one without them (the option
   GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT = 0) give the same printed
   values; dense elements passed to the other functions (square roots,
   polynomials, matrices, shallow copies) give the same results. */

#define DF_NUM 6

#define DF_NUM_SMALL 11
static const char * df_gens[] = {
    "sqrt(2)", "sqrt(-3)", "sqrt(-163)", "exp(2*pi*i/5)", "exp(2*pi*i/7)",
    "exp(2*pi*i/9)", "2^(1/3)", "exp(2*pi*i/16)", "sqrt(5)", "i", "(1+sqrt(5))/2",
    /* (degrees above 16, where sparse elements stay flat: the sparsity
       rule, and degrees above 64, where the products use fast division) */
    "exp(2*pi*i/37)", "2^(1/37)", "exp(2*pi*i/67)", "exp(2*pi*i/128)", "exp(2*pi*i/101)"
};

static int
df_random_op(gr_ptr * v, slong n, flint_rand_t state, gr_ctx_t K)
{
    slong i = n_randint(state, n), j = n_randint(state, n), k = n_randint(state, n);
    slong c = (slong) n_randint(state, 7) - 3;
    int status = GR_SUCCESS;
    fmpq_t q;
    fmpz_t z;

    fmpq_init(q);
    fmpz_init(z);
    fmpq_set_si(q, c, 1 + n_randint(state, 4));
    fmpz_set_si(z, c);

    switch (n_randint(state, 20))
    {
        case 0: status = gr_add(v[i], v[j], v[k], K); break;
        case 1: status = gr_sub(v[i], v[j], v[k], K); break;
        case 2: status = gr_mul(v[i], v[j], v[k], K); break;
        case 3: status = gr_div(v[i], v[j], v[k], K); break;
        case 4: status = gr_inv(v[i], v[j], K); break;
        case 5: status = gr_neg(v[i], v[j], K); break;
        case 6: status = gr_set(v[i], v[j], K); break;
        case 7: status = gr_add_si(v[i], v[j], c, K); break;
        case 8: status = gr_mul_si(v[i], v[j], c, K); break;
        case 9: status = gr_div_si(v[i], v[j], c, K); break;
        case 10: status = gr_add_fmpq(v[i], v[j], q, K); break;
        case 11: status = gr_mul_fmpq(v[i], v[j], q, K); break;
        case 12: status = gr_sub_fmpz(v[i], v[j], z, K); break;
        case 13: status = gr_set_si(v[i], c, K); break;
        case 14: status = gr_mul(v[i], v[i], v[i], K); break;
        case 15: status = gr_add(v[i], v[i], v[j], K); break;
        case 16: status = gr_sqr(v[i], v[j], K); break;
        case 17: status = gr_pow_ui(v[i], v[j], n_randint(state, 5), K); break;
        case 18: status = gr_div_fmpz(v[i], v[j], z, K); break;
        default: status = gr_sub_fmpq(v[i], v[j], q, K); break;
    }

    fmpq_clear(q);
    fmpz_clear(z);
    return status;
}

TEST_FUNCTION_START(gr_tower_dense_form, state)
{
    slong iter;

    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t QQ, K[2];
        gr_ptr v[2][DF_NUM];
        slong i, j, steps = 1 + n_randint(state, 40);
        flint_rand_t st[2];
        int status[2], w;
        const char * g1 = df_gens[n_randint(state, sizeof(df_gens) / sizeof(df_gens[0]))];
        /* (the second field among the small ones: composita of two
           large fields would make the test slow) */
        const char * g2 = df_gens[n_randint(state, DF_NUM_SMALL)];
        int two = (n_randint(state, 4) == 0);
        ulong seed1 = n_randlimb(state), seed2 = n_randlimb(state);

        gr_ctx_init_fmpq(QQ);
        for (w = 0; w < 2; w++)
        {
            gr_ctx_init_tower_lazy(K[w], QQ, GR_TOWER_MERGE_EXPRESS);
            if (w == 1)
                gr_tower_lazy_ctx_set_option(K[w], GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT, 0);
            for (i = 0; i < DF_NUM; i++)
                v[w][i] = gr_heap_init(K[w]);
            flint_rand_init(st[w]);
            flint_rand_set_seed(st[w], seed1, seed2);
        }

        for (w = 0; w < 2; w++)
        {
            status[w] = GR_SUCCESS;
            for (i = 0; i < DF_NUM; i++)
            {
                /* elements of one field (or two, sometimes) */
                status[w] |= gr_set_str(v[w][i], (two && i % 3 == 2) ? g2 : g1, K[w]);
                status[w] |= gr_add_si(v[w][i], v[w][i], (slong) n_randint(st[w], 5) - 2, K[w]);
                status[w] |= gr_div_ui(v[w][i], v[w][i], 1 + n_randint(st[w], 3), K[w]);
            }
            for (j = 0; j < steps; j++)
                status[w] |= df_random_op(v[w], DF_NUM, st[w], K[w]);
        }

        if (status[0] != status[1])
        {
            flint_printf("FAIL: status %d %d (%s, %s)\n", status[0], status[1], g1, g2);
            flint_abort();
        }

        if (status[0] == GR_SUCCESS)
        {
            for (i = 0; i < DF_NUM; i++)
            {
                char * s0, * s1;
                truth_t e;
                GR_MUST_SUCCEED(gr_get_str(&s0, v[0][i], K[0]));
                GR_MUST_SUCCEED(gr_get_str(&s1, v[1][i], K[1]));
                if (strcmp(s0, s1) != 0)
                {
                    flint_printf("FAIL: values (%s, %s)\n%s\n%s\n", g1, g2, s0, s1);
                    flint_abort();
                }
                flint_free(s0);
                flint_free(s1);

                /* consistency of the predicates with the dense forms */
                e = gr_equal(v[0][i], v[0][(i + 1) % DF_NUM], K[0]);
                if (e != gr_equal(v[1][i], v[1][(i + 1) % DF_NUM], K[1]) && e != T_UNKNOWN)
                {
                    flint_printf("FAIL: equal (%s, %s)\n", g1, g2);
                    flint_abort();
                }
                if (gr_is_zero(v[0][i], K[0]) != gr_is_zero(v[1][i], K[1]) ||
                    gr_is_one(v[0][i], K[0]) != gr_is_one(v[1][i], K[1]))
                {
                    flint_printf("FAIL: predicates (%s, %s)\n", g1, g2);
                    flint_abort();
                }
            }

            /* x * (1/x) = 1, and (x + y)^2 = x^2 + 2 x y + y^2 */
            {
                gr_ptr a, b, c;
                gr_ctx_struct * Kd = K[0];
                GR_TMP_INIT3(a, b, c, Kd);
                for (i = 0; i < DF_NUM; i++)
                {
                    gr_ptr x = v[0][i], y = v[0][(i + 2) % DF_NUM];
                    if (gr_is_zero(x, Kd) == T_FALSE)
                    {
                        GR_MUST_SUCCEED(gr_inv(a, x, Kd));
                        GR_MUST_SUCCEED(gr_mul(a, a, x, Kd));
                        if (gr_is_one(a, Kd) != T_TRUE)
                        {
                            flint_printf("FAIL: inverse (%s, %s)\n", g1, g2);
                            flint_abort();
                        }
                    }
                    GR_MUST_SUCCEED(gr_add(a, x, y, Kd));
                    GR_MUST_SUCCEED(gr_sqr(a, a, Kd));
                    GR_MUST_SUCCEED(gr_mul(b, x, y, Kd));
                    GR_MUST_SUCCEED(gr_mul_two(b, b, Kd));
                    GR_MUST_SUCCEED(gr_sqr(c, x, Kd));
                    GR_MUST_SUCCEED(gr_add(b, b, c, Kd));
                    GR_MUST_SUCCEED(gr_sqr(c, y, Kd));
                    GR_MUST_SUCCEED(gr_add(b, b, c, Kd));
                    if (gr_equal(a, b, Kd) != T_TRUE)
                    {
                        flint_printf("FAIL: square of a sum (%s, %s)\n", g1, g2);
                        flint_abort();
                    }
                    /* (a dense element through a locked function) */
                    if (gr_sqrt(c, a, Kd) == GR_SUCCESS)
                    {
                        GR_MUST_SUCCEED(gr_sqr(c, c, Kd));
                        if (gr_equal(c, a, Kd) != T_TRUE)
                        {
                            flint_printf("FAIL: sqrt (%s, %s)\n", g1, g2);
                            flint_abort();
                        }
                    }
                }
                GR_TMP_CLEAR3(a, b, c, Kd);
            }

            /* matrices and polynomials of the elements (shallow copies) */
            {
                gr_mat_t A, B, C, D;
                gr_poly_t P, Q, R, S;
                gr_mat_init(A, 2, 3, K[0]);
                gr_mat_init(B, 3, 2, K[0]);
                gr_mat_init(C, 2, 2, K[0]);
                gr_mat_init(D, 2, 2, K[1]);
                for (i = 0; i < 6; i++)
                {
                    GR_MUST_SUCCEED(gr_set(gr_mat_entry_ptr(A, i / 3, i % 3, K[0]), v[0][i], K[0]));
                    GR_MUST_SUCCEED(gr_set(gr_mat_entry_ptr(B, i / 2, i % 2, K[0]), v[0][(i + 1) % DF_NUM], K[0]));
                }
                GR_MUST_SUCCEED(gr_mat_mul(C, A, B, K[0]));
                gr_mat_clear(A, K[0]);
                gr_mat_clear(B, K[0]);
                gr_mat_init(A, 2, 3, K[1]);
                gr_mat_init(B, 3, 2, K[1]);
                for (i = 0; i < 6; i++)
                {
                    GR_MUST_SUCCEED(gr_set(gr_mat_entry_ptr(A, i / 3, i % 3, K[1]), v[1][i], K[1]));
                    GR_MUST_SUCCEED(gr_set(gr_mat_entry_ptr(B, i / 2, i % 2, K[1]), v[1][(i + 1) % DF_NUM], K[1]));
                }
                GR_MUST_SUCCEED(gr_mat_mul(D, A, B, K[1]));
                for (i = 0; i < 4; i++)
                {
                    char * s0, * s1;
                    GR_MUST_SUCCEED(gr_get_str(&s0, gr_mat_entry_ptr(C, i / 2, i % 2, K[0]), K[0]));
                    GR_MUST_SUCCEED(gr_get_str(&s1, gr_mat_entry_ptr(D, i / 2, i % 2, K[1]), K[1]));
                    if (strcmp(s0, s1) != 0)
                    {
                        /* (the generators may be ordered differently in the
                           two contexts, depending on which towers were
                           merged: then the values must agree) */
                        acb_t z0, z1;
                        acb_init(z0);
                        acb_init(z1);
                        if (gr_tower_lazy_get_acb(z0, gr_mat_entry_ptr(C, i / 2, i % 2, K[0]), 256, K[0]) != GR_SUCCESS ||
                            gr_tower_lazy_get_acb(z1, gr_mat_entry_ptr(D, i / 2, i % 2, K[1]), 256, K[1]) != GR_SUCCESS ||
                            !acb_overlaps(z0, z1) || acb_rel_accuracy_bits(z0) < 200)
                        {
                            flint_printf("FAIL: matrix product (%s, %s)\n%s\n%s\n", g1, g2, s0, s1);
                            flint_abort();
                        }
                        acb_clear(z0);
                        acb_clear(z1);
                    }
                    flint_free(s0);
                    flint_free(s1);
                }
                gr_mat_clear(A, K[1]);
                gr_mat_clear(B, K[1]);
                gr_mat_clear(C, K[0]);
                gr_mat_clear(D, K[1]);

                gr_poly_init(P, K[0]);
                gr_poly_init(Q, K[0]);
                gr_poly_init(R, K[1]);
                gr_poly_init(S, K[1]);
                for (i = 0; i < 3; i++)
                {
                    GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(P, i, v[0][i], K[0]));
                    GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(R, i, v[1][i], K[1]));
                }
                GR_MUST_SUCCEED(gr_poly_mul(Q, P, P, K[0]));
                GR_MUST_SUCCEED(gr_poly_mul(S, R, R, K[1]));
                for (i = 0; i < Q->length || i < S->length; i++)
                {
                    char * s0, * s1;
                    gr_ptr c0 = gr_heap_init(K[0]), c1 = gr_heap_init(K[1]);
                    GR_MUST_SUCCEED(gr_poly_get_coeff_scalar(c0, Q, i, K[0]));
                    GR_MUST_SUCCEED(gr_poly_get_coeff_scalar(c1, S, i, K[1]));
                    GR_MUST_SUCCEED(gr_get_str(&s0, c0, K[0]));
                    GR_MUST_SUCCEED(gr_get_str(&s1, c1, K[1]));
                    if (strcmp(s0, s1) != 0)
                    {
                        /* (the generators may be ordered differently in
                           the two contexts: compare the values exactly) */
                        gr_ptr c2 = gr_heap_init(K[0]);
                        if (gr_set_other(c2, c1, K[1], K[0]) != GR_SUCCESS || gr_equal(c0, c2, K[0]) != T_TRUE)
                        {
                            flint_printf("FAIL: polynomial product (%s, %s)\n%s\n%s\n", g1, g2, s0, s1);
                            flint_abort();
                        }
                        gr_heap_clear(c2, K[0]);
                    }
                    flint_free(s0);
                    flint_free(s1);
                    gr_heap_clear(c0, K[0]);
                    gr_heap_clear(c1, K[1]);
                }
                gr_poly_clear(P, K[0]);
                gr_poly_clear(Q, K[0]);
                gr_poly_clear(R, K[1]);
                gr_poly_clear(S, K[1]);
            }
        }

        for (w = 0; w < 2; w++)
        {
            for (i = 0; i < DF_NUM; i++)
                gr_heap_clear(v[w][i], K[w]);
            flint_rand_clear(st[w]);
            gr_ctx_clear(K[w]);
        }
        gr_ctx_clear(QQ);
    }

    TEST_FUNCTION_END(state);
}
