/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include "test_helpers.h"
#include "fmpq.h"
#include "qqbar.h"
#include "acb.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    Random algebraic expressions (rationals, i, roots of unity, radicals
    of small integers, sums, products, quotients, square and cube roots,
    conjugates, real parts) evaluated in the lazy tower field and in
    qqbar, which serves as an independent oracle: the two values must
    agree, both after converting the tower value to qqbar and by the
    zero test in the tower field (which merges the expression's tower
    with the tower built from the qqbar value).
*/

#define DEGREE_LIMIT 256

static char trace[10000];

/* a random expression, evaluated in both fields; depth bounds the tree */
static int
_random_expr(gr_ptr x, gr_ptr q, slong depth, flint_rand_t state, gr_ctx_t K, gr_ctx_t QQbar)
{
    int status = GR_SUCCESS;
    int op = n_randint(state, depth == 0 ? 5 : 12);

    if (strlen(trace) < 9000)
        flint_sprintf(trace + strlen(trace), "%d ", op);

    switch (op)
    {
        case 0:   /* small rational */
        {
            fmpq_t c;
            fmpq_init(c);
            fmpq_set_si(c, (slong) n_randint(state, 11) - 5, 1 + n_randint(state, 4));
            status |= gr_set_fmpq(x, c, K);
            status |= gr_set_fmpq(q, c, QQbar);
            fmpq_clear(c);
            break;
        }
        case 1:   /* i */
            status |= gr_i(x, K);
            status |= gr_i(q, QQbar);
            break;
        case 2:   /* root of unity exp(2 pi i k / n) */
        {
            ulong n = 1 + n_randint(state, 12);
            ulong k = n_randint(state, n);
            qqbar_t z;
            if (strlen(trace) < 9000)
                flint_sprintf(trace + strlen(trace), "(z %wu/%wu) ", k, n);
            qqbar_init(z);
            qqbar_root_of_unity(z, k, n);
            status |= gr_set_other(q, z, QQbar, QQbar);
            status |= gr_set_other(x, z, QQbar, K);
            qqbar_clear(z);
            break;
        }
        case 3:   /* n-th root of a small positive integer */
        case 4:
        {
            ulong n = 2 + n_randint(state, 3);
            ulong a = 1 + n_randint(state, 12);
            fmpq_t e;
            if (strlen(trace) < 9000)
                flint_sprintf(trace + strlen(trace), "(%wu^(1/%wu)) ", a, n);
            fmpq_init(e);
            fmpq_set_si(e, 1, n);
            status |= gr_set_ui(x, a, K);
            status |= gr_pow_fmpq(x, x, e, K);
            status |= gr_set_ui(q, a, QQbar);
            status |= gr_pow_fmpq(q, q, e, QQbar);
            fmpq_clear(e);
            break;
        }
        case 5: case 6:   /* sum, difference */
        {
            gr_ptr x2, q2;
            x2 = gr_heap_init(K);
            q2 = gr_heap_init(QQbar);
            status |= _random_expr(x, q, depth - 1, state, K, QQbar);
            status |= _random_expr(x2, q2, depth - 1, state, K, QQbar);
            if (qqbar_degree(q) * qqbar_degree(q2) > DEGREE_LIMIT)
            {
                /* (beyond the reach of the qqbar oracle: keep one operand) */
            }
            else if (op == 5)
            {
                status |= gr_add(x, x, x2, K);
                status |= gr_add(q, q, q2, QQbar);
            }
            else
            {
                status |= gr_sub(x, x, x2, K);
                status |= gr_sub(q, q, q2, QQbar);
            }
            gr_heap_clear(x2, K);
            gr_heap_clear(q2, QQbar);
            break;
        }
        case 7: case 8:   /* product, quotient */
        {
            gr_ptr x2, q2;
            x2 = gr_heap_init(K);
            q2 = gr_heap_init(QQbar);
            status |= _random_expr(x, q, depth - 1, state, K, QQbar);
            status |= _random_expr(x2, q2, depth - 1, state, K, QQbar);
            if (qqbar_degree(q) * qqbar_degree(q2) > DEGREE_LIMIT)
            {
            }
            else if (op == 7 || gr_is_zero(q2, QQbar) == T_TRUE)
            {
                status |= gr_mul(x, x, x2, K);
                status |= gr_mul(q, q, q2, QQbar);
            }
            else
            {
                status |= gr_div(x, x, x2, K);
                status |= gr_div(q, q, q2, QQbar);
            }
            gr_heap_clear(x2, K);
            gr_heap_clear(q2, QQbar);
            break;
        }
        case 9:   /* principal square or cube root */
        {
            ulong n = 2 + n_randint(state, 2);
            fmpq_t e;
            if (strlen(trace) < 9000)
                flint_sprintf(trace + strlen(trace), "(root %wu) ", n);
            fmpq_init(e);
            fmpq_set_si(e, 1, n);
            status |= _random_expr(x, q, depth - 1, state, K, QQbar);
            if (qqbar_degree(q) * n <= DEGREE_LIMIT)
            {
                status |= gr_pow_fmpq(x, x, e, K);
                status |= gr_pow_fmpq(q, q, e, QQbar);
            }
            fmpq_clear(e);
            break;
        }
        case 10:  /* conjugate */
            status |= _random_expr(x, q, depth - 1, state, K, QQbar);
            status |= gr_conj(x, x, K);
            status |= gr_conj(q, q, QQbar);
            break;
        case 11:  /* real part (of degree up to the square of the degree) */
            status |= _random_expr(x, q, depth - 1, state, K, QQbar);
            if (qqbar_degree(q) * qqbar_degree(q) <= DEGREE_LIMIT)
            {
                status |= gr_re(x, x, K);
                status |= gr_re(q, q, QQbar);
            }
            break;
    }

    return status;
}

TEST_FUNCTION_START(gr_tower_lazy_qqbar, state)
{
    gr_ctx_t QQ, K, QQbar;
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_complex_qqbar(QQbar);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ptr x, y, d, q, r;
        slong depth = 1 + n_randint(state, 4);
        int status;

        /* a fresh context every few iterations, otherwise the registry
           of towers accumulates (which is also worth testing) */
        if (iter % 10 == 0)
        {
            if (iter > 0)
                gr_ctx_clear(K);
            gr_ctx_init_tower_lazy(K, QQ, n_randint(state, 2) ? GR_TOWER_MERGE_EXPRESS : 0);
            /* (sometimes with the generator policies: composite roots of
               unity, composite square roots) */
            if (n_randint(state, 3) == 0)
                gr_tower_lazy_ctx_set_gen_flags(K, 2 * (1 + n_randint(state, 3)));
            /* (sometimes with a random cyclotomic degree cap) */
            if (n_randint(state, 4) == 0)
                gr_tower_lazy_ctx_set_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT, 1 << n_randint(state, 9));
            /* (sometimes with primitive elements for small number fields) */
            if (n_randint(state, 4) == 0)
                gr_tower_lazy_ctx_set_option(K, GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT, 4 << n_randint(state, 3));
        }

        x = gr_heap_init(K);
        y = gr_heap_init(K);
        d = gr_heap_init(K);
        q = gr_heap_init(QQbar);
        r = gr_heap_init(QQbar);

        trace[0] = 0;
        status = _random_expr(x, q, depth, state, K, QQbar);

        if (status == GR_SUCCESS)
        {
            /* the tower value as a qqbar (for towers of moderate degree:
               the coordinates of the powers of an element of a field of
               degree in the hundreds are large) */
            {
                slong level;
                status = (gr_tower_degree(gr_tower_lazy_get_tower(&level, x, K)) <= 64) ?
                    gr_set_other(r, x, K, QQbar) : GR_UNABLE;
            }
            if (status == GR_SUCCESS && gr_equal(r, q, QQbar) != T_TRUE)
            {
                flint_printf("FAIL: values differ (qqbar conversion)\n%s\n", trace);
                flint_printf("tower: "); gr_println(x, K);
                flint_printf("       "); gr_println(r, QQbar);
                flint_printf("qqbar: "); gr_println(q, QQbar);
                flint_abort();
            }

            /* the qqbar value in the tower field, compared exactly there
               (for moderate degrees: a qqbar of large degree is adjoined
               as a generator of its own, and the merged tower is large) */
            {
                slong level;
                status = (qqbar_degree(q) <= 24 && gr_tower_degree(gr_tower_lazy_get_tower(&level, x, K)) <= 64) ?
                    gr_set_other(y, q, QQbar, K) : GR_UNABLE;
            }
            if (status == GR_SUCCESS)
            {
                status = gr_sub(d, x, y, K);
                if (status == GR_SUCCESS && gr_is_zero(d, K) != T_TRUE)
                {
                    flint_printf("FAIL: values differ (zero test in the tower field)\n%s\n", trace);
                    flint_printf("tower: "); gr_println(x, K);
                    flint_printf("qqbar: "); gr_println(q, QQbar);
                    flint_printf("diff:  "); gr_println(d, K);
                    flint_abort();
                }

                /* and a nonzero difference is found nonzero */
                if (status == GR_SUCCESS && gr_is_zero(q, QQbar) != T_TRUE)
                {
                    status = gr_add(d, d, y, K);
                    if (status == GR_SUCCESS && gr_is_zero(d, K) != T_FALSE)
                    {
                        flint_printf("FAIL: nonzero value not found nonzero\n%s\n", trace);
                        flint_printf("tower: "); gr_println(x, K);
                        flint_printf("qqbar: "); gr_println(q, QQbar);
                        flint_abort();
                    }
                }
            }

            /* enclosures agree */
            {
                acb_t z, w;
                acb_init(z);
                acb_init(w);
                if (gr_tower_lazy_get_acb(z, x, 64, K) == GR_SUCCESS)
                {
                    qqbar_get_acb(w, q, 64);
                    if (!acb_overlaps(z, w))
                    {
                        flint_printf("FAIL: enclosures do not overlap\n%s\n", trace);
                        flint_printf("tower: "); gr_println(x, K); acb_printd(z, 15); flint_printf("\n");
                        flint_printf("qqbar: "); gr_println(q, QQbar); acb_printd(w, 15); flint_printf("\n");
                        flint_abort();
                    }
                }
                acb_clear(z);
                acb_clear(w);
            }
        }

        gr_heap_clear(x, K);
        gr_heap_clear(y, K);
        gr_heap_clear(d, K);
        gr_heap_clear(q, QQbar);
        gr_heap_clear(r, QQbar);
    }

    gr_ctx_clear(K);

    /* square roots of integers next to roots of unity whose field
       contains them (Gauss sums): sqrt(3) with zeta_3 and i, sqrt(2)
       with zeta_8, sqrt(5) with zeta_5, sqrt(-7) with zeta_7 */
    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        static const slong cs[] = { 2, 3, 5, 6, 7, 10, 11, 13, -1, -2, -3, -5, -7, 12, 18 };
        gr_ptr s, z, w, a, b;
        slong c = cs[n_randint(state, 15)], c2, N, k;
        char str[100];

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        s = gr_heap_init(K);
        z = gr_heap_init(K);
        w = gr_heap_init(K);
        a = gr_heap_init(K);
        b = gr_heap_init(K);

        /* zeta_N with 4 * 8 * 3 * 5 * 7 * 11 * 13 / N small enough,
           times zeta_37 */
        {
            static const slong Ns[] = { 24, 40, 56, 88, 104, 120, 168, 12, 20, 28 };
            N = Ns[n_randint(state, 10)];
        }
        flint_sprintf(str, "exp(2*pi*i/%wd)*exp(2*pi*i/37)", N);
        GR_MUST_SUCCEED(gr_set_str(z, str, K));
        GR_MUST_SUCCEED(gr_set_si(s, c, K));
        GR_MUST_SUCCEED(gr_sqrt(s, s, K));

        /* (s + z)^2 = c + 2 s z + z^2, the factors combined in different
           orders */
        GR_MUST_SUCCEED(gr_add(a, s, z, K));
        GR_MUST_SUCCEED(gr_sqr(a, a, K));
        GR_MUST_SUCCEED(gr_mul(b, s, z, K));
        GR_MUST_SUCCEED(gr_mul_ui(b, b, 2, K));
        GR_MUST_SUCCEED(gr_sqr(w, z, K));
        GR_MUST_SUCCEED(gr_add(b, b, w, K));
        GR_MUST_SUCCEED(gr_add_si(b, b, c, K));
        if (gr_equal(a, b, K) != T_TRUE)
        {
            flint_printf("FAIL: Gauss sums, c = %wd, N = %wd\n", c, N);
            gr_println(a, K);
            gr_println(b, K);
            flint_abort();
        }

        /* sqrt(c) z != sqrt(c2) z */
        c2 = cs[n_randint(state, 15)];
        k = n_randint(state, 4);
        if (c2 != c)
        {
            GR_MUST_SUCCEED(gr_set_si(w, c2, K));
            GR_MUST_SUCCEED(gr_sqrt(w, w, K));
            GR_MUST_SUCCEED(gr_mul(a, s, z, K));
            GR_MUST_SUCCEED(gr_mul(b, w, z, K));
            GR_MUST_SUCCEED(gr_pow_ui(a, a, k + 1, K));
            GR_MUST_SUCCEED(gr_pow_ui(b, b, k + 1, K));
            /* (equal only for (k + 1) = 0 mod 4 and c = -c2) */
            if (gr_equal(a, b, K) != (((k + 1) % 4 == 0 && c == -c2) ? T_TRUE : T_FALSE))
            {
                flint_printf("FAIL: Gauss sums (nonzero), c = %wd, c2 = %wd, N = %wd\n", c, c2, N);
                flint_abort();
            }
        }

        gr_heap_clear(s, K);
        gr_heap_clear(z, K);
        gr_heap_clear(w, K);
        gr_heap_clear(a, K);
        gr_heap_clear(b, K);
        gr_ctx_clear(K);
    }

    /* the generator policies: zeta_N of a composite order as one
       generator (of the least common multiple of the orders requested),
       sqrt(A) of a composite A as one generator; identities across them */
    {
        static const char * ids[] = {
            "exp(2*pi*i/15)^5 - exp(2*pi*i/3)",
            "exp(2*pi*i/12) - exp(2*pi*i/60)^5",
            "exp(2*pi*i/7)*exp(2*pi*i/9) - exp(2*pi*i*16/63)",
            "sqrt(5) - 1 - 2*(exp(2*pi*i/5) + exp(-2*pi*i/5))",
            "sqrt(2) - exp(2*pi*i/8) - exp(-2*pi*i/8)",
            "sqrt(-3) - 2*exp(2*pi*i/3) - 1",
            "sqrt(10)*sqrt(15) - 5*sqrt(6)",
            "sqrt(6) - sqrt(2)*sqrt(3)",
            "sqrt(12)/sqrt(3) - 2",
            "conj(exp(2*pi*i/35)) - exp(-2*pi*i/35)",
        };
        int flags;
        slong j;
        char * str;

        for (flags = 2; flags <= 6; flags += 2)
        {
            gr_ptr x;
            gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
            gr_tower_lazy_ctx_set_gen_flags(K, flags);
            x = gr_heap_init(K);

            if (flags & GR_TOWER_GENS_COMPOSITE_ROOTS)
            {
                /* (one generator: exp(2 pi i / 60) after i and zeta_15, in the
                   fresh context) */
                GR_MUST_SUCCEED(gr_set_str(x, "exp(2*pi*i/15)*i", K));
                GR_MUST_SUCCEED(gr_get_str(&str, x, K));
                if (strstr(str, "exp(2*pi*i/60)") == NULL || strstr(str, "exp(2*pi*i/15)") != NULL)
                {
                    flint_printf("FAIL: composite root of unity: %s\n", str);
                    flint_abort();
                }
                flint_free(str);
            }
            if (flags & GR_TOWER_GENS_COMPOSITE_RADICALS)
            {
                GR_MUST_SUCCEED(gr_set_str(x, "sqrt(24)", K));
                GR_MUST_SUCCEED(gr_get_str(&str, x, K));
                if (strstr(str, "sqrt(6)") == NULL)
                {
                    flint_printf("FAIL: composite square root: %s\n", str);
                    flint_abort();
                }
                flint_free(str);
            }

            for (j = 0; j < (slong) (sizeof(ids) / sizeof(ids[0])); j++)
            {
                if (gr_set_str(x, ids[j], K) != GR_SUCCESS || gr_is_zero(x, K) != T_TRUE)
                {
                    flint_printf("FAIL: generator policy %d: %s\n", flags, ids[j]);
                    flint_abort();
                }
                GR_MUST_SUCCEED(gr_set_str(x, ids[j], K));
                GR_MUST_SUCCEED(gr_add_ui(x, x, 1, K));
                if (gr_is_zero(x, K) != T_FALSE)
                {
                    flint_printf("FAIL: generator policy %d (perturbed): %s\n", flags, ids[j]);
                    flint_abort();
                }
            }

            gr_heap_clear(x, K);
            gr_ctx_clear(K);
        }
    }

    /* the options: the cyclotomic degree cap (prime powers beyond it),
       the tangent form in the complex field */
    {
        gr_ptr x;
        char * str;
        slong cap;

        for (cap = 8; cap <= 16; cap += 8)
        {
            gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
            gr_tower_lazy_ctx_set_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT, cap);
            if (gr_tower_lazy_ctx_get_option(K, GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT) != cap)
            {
                flint_printf("FAIL: get_option\n");
                flint_abort();
            }
            x = gr_heap_init(K);
            GR_MUST_SUCCEED(gr_set_str(x, "exp(2*pi*i/15)*i", K));
            GR_MUST_SUCCEED(gr_get_str(&str, x, K));
            if ((strstr(str, "exp(2*pi*i/60)") != NULL) != (cap == 16))
            {
                flint_printf("FAIL: cyclotomic degree cap %wd: %s\n", cap, str);
                flint_abort();
            }
            flint_free(str);
            GR_MUST_SUCCEED(gr_set_str(x, "exp(2*pi*i/15)^15 - 1", K));
            if (gr_is_zero(x, K) != T_TRUE)
            {
                flint_printf("FAIL: cyclotomic degree cap %wd (zero)\n", cap);
                flint_abort();
            }
            gr_heap_clear(x, K);
            gr_ctx_clear(K);
        }

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_tower_lazy_ctx_set_option(K, GR_TOWER_OPT_TRIG_FORM, GR_TOWER_TRIG_TANGENT);
        x = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_set_str(x, "pi/5", K));
        GR_MUST_SUCCEED(gr_cos(x, x, K));
        GR_MUST_SUCCEED(gr_get_str(&str, x, K));
        if (strstr(str, "tan(pi/5)") == NULL)
        {
            flint_printf("FAIL: tangent form: %s\n", str);
            flint_abort();
        }
        flint_free(str);
        GR_MUST_SUCCEED(gr_set_str(x, "cos(pi/5) - (1 + sqrt(5))/4", K));
        if (gr_is_zero(x, K) != T_TRUE)
        {
            flint_printf("FAIL: tangent form (zero)\n");
            flint_abort();
        }
        gr_heap_clear(x, K);
        gr_ctx_clear(K);
    }

    /* inverses and quotients in towers with triangular moduli (the
       linear algebra inverse) */
    {
        const char * towers[][3] = {
            { "sqrt(2)", "sqrt(1+sqrt(2))", "3" },
            { "sqrt(2)", "sqrt(2+sqrt(2))", "sqrt(2+sqrt(2+sqrt(2)))" },
            { "sqrt(3)", "(1+sqrt(3))^(1/3)", "sqrt(5)" },
            { "2^(1/3)", "sqrt(1+2^(1/3))", "sqrt(3+sqrt(1+2^(1/3)))" },
        };
        slong w, rep, k;

        for (w = 0; w < 4; w++)
        {
            gr_ptr g[3], x, y, t, u;
            gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
            for (k = 0; k < 3; k++)
            {
                g[k] = gr_heap_init(K);
                GR_MUST_SUCCEED(gr_set_str(g[k], towers[w][k], K));
            }
            x = gr_heap_init(K); y = gr_heap_init(K); t = gr_heap_init(K); u = gr_heap_init(K);
            for (rep = 0; rep < 5; rep++)
            {
                GR_MUST_SUCCEED(gr_set_si(x, 1 + n_randint(state, 5), K));
                GR_MUST_SUCCEED(gr_set_si(y, (slong) n_randint(state, 5) - 2, K));
                for (k = 0; k < 3; k++)
                {
                    GR_MUST_SUCCEED(gr_mul_si(t, g[k], (slong) n_randint(state, 7) - 3, K));
                    GR_MUST_SUCCEED(gr_add(x, x, t, K));
                    GR_MUST_SUCCEED(gr_mul(t, t, g[(k + 1) % 3], K));
                    GR_MUST_SUCCEED(gr_add(y, y, t, K));
                }
                if (gr_is_zero(x, K) != T_FALSE)
                    continue;
                GR_MUST_SUCCEED(gr_inv(t, x, K));
                GR_MUST_SUCCEED(gr_mul(u, t, x, K));
                if (gr_is_one(u, K) != T_TRUE)
                {
                    flint_printf("FAIL: triangular inverse (%s, %s, %s)\n", towers[w][0], towers[w][1], towers[w][2]);
                    flint_abort();
                }
                GR_MUST_SUCCEED(gr_div(t, y, x, K));
                GR_MUST_SUCCEED(gr_mul(u, t, x, K));
                if (gr_equal(u, y, K) != T_TRUE)
                {
                    flint_printf("FAIL: triangular quotient (%s, %s, %s)\n", towers[w][0], towers[w][1], towers[w][2]);
                    flint_abort();
                }
            }
            for (k = 0; k < 3; k++)
                gr_heap_clear(g[k], K);
            gr_heap_clear(x, K); gr_heap_clear(y, K); gr_heap_clear(t, K); gr_heap_clear(u, K);
            gr_ctx_clear(K);
        }
    }

    gr_ctx_clear(QQbar);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
