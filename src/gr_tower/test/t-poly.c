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
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    Polynomials over the lazy field: gcds and resultants (by the
    subresultant PRS, the Euclidean sequence over Q(pi, e)(sqrt 2, ...)
    having exponential coefficient growth), and roots of quadratics,
    binomials and products of known factors.
*/

/* a random element: a small combination of sqrt 2, sqrt 3, pi, e */
static int
_random_elem(gr_ptr res, flint_rand_t state, gr_ptr * gens, gr_ctx_t K)
{
    int status = GR_SUCCESS;
    slong i;
    status |= gr_set_si(res, (slong) n_randint(state, 7) - 3, K);
    for (i = 0; i < 4; i++)
    {
        if (n_randint(state, 3) == 0)
        {
            gr_ptr t;
            GR_TMP_INIT(t, K);
            status |= gr_mul_si(t, gens[i], (slong) n_randint(state, 5) - 2, K);
            status |= gr_add(res, res, t, K);
            GR_TMP_CLEAR(t, K);
        }
    }
    return status;
}

static int
_random_poly(gr_poly_t f, slong len, flint_rand_t state, gr_ptr * gens, gr_ctx_t K)
{
    slong i;
    int status = GR_SUCCESS;
    gr_poly_fit_length(f, len, K);
    for (i = 0; i < len; i++)
        status |= _random_elem(gr_poly_coeff_ptr(f, i, K), state, gens, K);
    _gr_poly_set_length(f, len, K);
    _gr_poly_normalise(f, K);
    return status;
}

TEST_FUNCTION_START(gr_tower_poly, state)
{
    gr_ctx_t QQ, K;
    gr_ptr gens[4];
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);

    gens[0] = gr_heap_init(K);
    gens[1] = gr_heap_init(K);
    gens[2] = gr_heap_init(K);
    gens[3] = gr_heap_init(K);
    GR_MUST_SUCCEED(gr_set_ui(gens[0], 2, K));
    GR_MUST_SUCCEED(gr_sqrt(gens[0], gens[0], K));
    GR_MUST_SUCCEED(gr_set_ui(gens[1], 3, K));
    GR_MUST_SUCCEED(gr_sqrt(gens[1], gens[1], K));
    GR_MUST_SUCCEED(gr_pi(gens[2], K));
    GR_MUST_SUCCEED(gr_one(gens[3], K));
    GR_MUST_SUCCEED(gr_exp(gens[3], gens[3], K));

    /* gcd(g a, g b) = g (up to a unit) for random a, b, g with coprime
       a, b (generically); the resultant of g a and g b vanishes */
    for (iter = 0; iter < 6 * flint_test_multiplier(); iter++)
    {
        gr_poly_t a, b, g, ga, gb, d, r;
        gr_ptr res;
        int status = GR_SUCCESS;

        gr_poly_init(a, K);
        gr_poly_init(b, K);
        gr_poly_init(g, K);
        gr_poly_init(ga, K);
        gr_poly_init(gb, K);
        gr_poly_init(d, K);
        gr_poly_init(r, K);
        res = gr_heap_init(K);

        status |= _random_poly(a, 2 + n_randint(state, 3), state, gens, K);
        status |= _random_poly(b, 2 + n_randint(state, 3), state, gens, K);
        status |= _random_poly(g, 2 + n_randint(state, 3), state, gens, K);
        status |= gr_poly_mul(ga, g, a, K);
        status |= gr_poly_mul(gb, g, b, K);

        if (status == GR_SUCCESS && a->length >= 2 && b->length >= 2 && g->length >= 2)
        {
            status = gr_poly_gcd(d, ga, gb, K);
            if (status == GR_SUCCESS)
            {
                /* g divides d */
                gr_poly_t q;
                gr_poly_init(q, K);
                status = gr_poly_divrem(q, r, d, g, K);
                if (status == GR_SUCCESS && !(gr_poly_is_zero(r, K) == T_TRUE))
                {
                    flint_printf("FAIL: gcd not divisible by the common factor\n");
                    flint_printf("g = "); gr_poly_print(g, K); flint_printf("\n");
                    flint_printf("d = "); gr_poly_print(d, K); flint_printf("\n");
                    flint_abort();
                }
                /* d divides ga and gb */
                status |= gr_poly_divrem(q, r, ga, d, K);
                if (status == GR_SUCCESS && !(gr_poly_is_zero(r, K) == T_TRUE))
                {
                    flint_printf("FAIL: gcd does not divide\n");
                    flint_abort();
                }
                gr_poly_clear(q, K);
            }

            if (status == GR_SUCCESS)
            {
                status = gr_poly_resultant(res, ga, gb, K);
                if (status == GR_SUCCESS && gr_is_zero(res, K) != T_TRUE)
                {
                    flint_printf("FAIL: resultant of polynomials with a common factor\n");
                    flint_abort();
                }
            }
        }

        gr_poly_clear(a, K);
        gr_poly_clear(b, K);
        gr_poly_clear(g, K);
        gr_poly_clear(ga, K);
        gr_poly_clear(gb, K);
        gr_poly_clear(d, K);
        gr_poly_clear(r, K);
        gr_heap_clear(res, K);
    }

    /* roots: products of linear factors with known roots, quadratics,
       binomials; every root is verified exactly */
    for (iter = 0; iter < 6 * flint_test_multiplier(); iter++)
    {
        gr_poly_t f, lin;
        gr_vec_t roots;
        fmpz_vec_t mult;
        gr_ptr c, v;
        slong n = 1 + n_randint(state, 4), i, expected;
        int status = GR_SUCCESS;
        int kind = n_randint(state, 3);

        gr_poly_init(f, K);
        gr_poly_init(lin, K);
        gr_vec_init(roots, 0, K);
        fmpz_vec_init(mult, 0);
        c = gr_heap_init(K);
        v = gr_heap_init(K);

        if (kind == 0)
        {
            /* product of x - c_i, distinct c_i generically */
            status |= gr_poly_one(f, K);
            for (i = 0; i < n; i++)
            {
                status |= _random_elem(c, state, gens, K);
                status |= gr_poly_zero(lin, K);
                status |= gr_poly_set_coeff_si(lin, 1, 1, K);
                status |= gr_neg(c, c, K);
                status |= gr_poly_set_coeff_scalar(lin, 0, c, K);
                status |= gr_poly_mul(f, f, lin, K);
            }
            expected = -1;   /* (repeated roots possible) */
        }
        else if (kind == 1)
        {
            /* a quadratic */
            status |= _random_poly(f, 3, state, gens, K);
            expected = (f->length == 3) ? 2 : -1;
        }
        else
        {
            /* a x^n + b */
            n = 2 + n_randint(state, 4);
            status |= gr_poly_zero(f, K);
            status |= _random_elem(c, state, gens, K);
            status |= gr_poly_set_coeff_scalar(f, 0, c, K);
            status |= _random_elem(c, state, gens, K);
            if (gr_is_zero(c, K) == T_TRUE)
                status |= gr_one(c, K);
            status |= gr_poly_set_coeff_scalar(f, n, c, K);
            expected = (gr_is_zero(gr_poly_coeff_srcptr(f, 0, K), K) == T_TRUE) ? -1 : n;
        }

        if (status == GR_SUCCESS && f->length >= 2)
        {
            status = gr_poly_roots(roots, mult, f, 0, K);
            if (status == GR_SUCCESS)
            {
                slong total = 0;
                for (i = 0; i < roots->length; i++)
                {
                    status |= gr_poly_evaluate(v, f, gr_vec_entry_ptr(roots, i, K), K);
                    if (status == GR_SUCCESS && gr_is_zero(v, K) != T_TRUE)
                    {
                        flint_printf("FAIL: not a root\n");
                        flint_printf("f = "); gr_poly_print(f, K); flint_printf("\n");
                        flint_printf("r = "); gr_println(gr_vec_entry_ptr(roots, i, K), K);
                        flint_abort();
                    }
                    total += fmpz_get_si(mult->entries + i);
                }
                if (status == GR_SUCCESS && total != f->length - 1)
                {
                    flint_printf("FAIL: number of roots %wd for degree %wd\n", total, f->length - 1);
                    flint_printf("f = "); gr_poly_print(f, K); flint_printf("\n");
                    flint_abort();
                }
                if (status == GR_SUCCESS && expected >= 0 && roots->length != expected)
                {
                    flint_printf("FAIL: %wd distinct roots, expected %wd\n", roots->length, expected);
                    flint_printf("f = "); gr_poly_print(f, K); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        gr_poly_clear(f, K);
        gr_poly_clear(lin, K);
        gr_vec_clear(roots, K);
        fmpz_vec_clear(mult);
        gr_heap_clear(c, K);
        gr_heap_clear(v, K);
    }

    /* structured roots: sqrt(4 sqrt 2) = 2 root_4(2), the roots of
       x^4 - sqrt 2 are root_8(2) times the fourth roots of unity */
    {
        gr_ptr x, y;
        gr_poly_t f;
        gr_vec_t roots;
        fmpz_vec_t mult;
        fmpq_t q;

        x = gr_heap_init(K);
        y = gr_heap_init(K);
        fmpq_init(q);
        GR_MUST_SUCCEED(gr_mul_ui(x, gens[0], 4, K));
        GR_MUST_SUCCEED(gr_sqrt(x, x, K));
        fmpq_set_si(q, 1, 4);
        GR_MUST_SUCCEED(gr_set_ui(y, 2, K));
        GR_MUST_SUCCEED(gr_pow_fmpq(y, y, q, K));
        GR_MUST_SUCCEED(gr_mul_ui(y, y, 2, K));
        if (gr_equal(x, y, K) != T_TRUE)
        {
            flint_printf("FAIL: sqrt(4 sqrt 2)\n");
            flint_abort();
        }
        {
            char * s;
            GR_MUST_SUCCEED(gr_get_str(&s, x, K));
            if (strstr(s, "root(2, 4)") == NULL)
            {
                flint_printf("FAIL: sqrt(4 sqrt 2) not structured: %s\n", s);
                flint_abort();
            }
            flint_free(s);
        }

        gr_poly_init(f, K);
        gr_vec_init(roots, 0, K);
        fmpz_vec_init(mult, 0);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 4, 1, K));
        GR_MUST_SUCCEED(gr_neg(y, gens[0], K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, 0, y, K));
        GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K));
        if (roots->length != 4)
        {
            flint_printf("FAIL: roots of x^4 - sqrt 2\n");
            flint_abort();
        }
        {
            slong lv;
            gr_tower_struct * T = gr_tower_lazy_get_tower(&lv, gr_vec_entry_ptr(roots, 0, K), K);
            if (gr_tower_degree(T) > 16)
            {
                flint_printf("FAIL: roots of x^4 - sqrt 2 in a tower of degree %wd\n", gr_tower_degree(T));
                flint_abort();
            }
        }
        gr_poly_clear(f, K);
        gr_vec_clear(roots, K);
        fmpz_vec_clear(mult);
        gr_heap_clear(x, K);
        gr_heap_clear(y, K);
        fmpq_clear(q);
    }

    gr_heap_clear(gens[0], K);
    gr_heap_clear(gens[1], K);
    gr_heap_clear(gens[2], K);
    gr_heap_clear(gens[3], K);
    gr_ctx_clear(K);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
