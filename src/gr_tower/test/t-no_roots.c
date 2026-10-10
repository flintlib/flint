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
#include "fmpz_poly.h"
#include "gr_poly.h"
#include "acb.h"
#include "ulong_extras.h"
#include "gr_tower.h"

/*
    gr_tower_poly_no_roots_modular never excludes the roots of a
    polynomial having a root in the tower ((x - r) h for a random element
    r, also with transcendental generators), and excludes them for
    polynomials without roots most of the time. gr_tower_set_subset keeps
    the marked generators with their moduli and enclosures.
*/

/* a random tower of radicals, optionally over pi */
static int
_random_tower(gr_tower_t T, flint_rand_t state, int with_pi)
{
    slong steps = 1 + n_randint(state, 3), k;
    int status = GR_SUCCESS;

    if (with_pi)
        status |= gr_tower_adjoin_pi(T, NULL);

    for (k = 0; k < steps && status == GR_SUCCESS; k++)
    {
        gr_ctx_struct * top = gr_tower_field(T);
        gr_ptr x;
        ulong n = 2 + n_randint(state, 2);

        GR_TMP_INIT(x, top);
        if (k == 0 || n_randint(state, 2))
            status |= gr_set_ui(x, 2 + n_randint(state, 20), top);
        else
        {
            status |= gr_gen(x, top);
            status |= gr_add_ui(x, x, 1 + n_randint(state, 5), top);
        }
        status |= gr_tower_adjoin_root_ui(T, x, n, NULL);
        GR_TMP_CLEAR(x, top);
    }

    return status;
}

TEST_FUNCTION_START(gr_tower_no_roots, state)
{
    gr_ctx_t QQ;
    slong iter, excluded = 0, total = 0;

    gr_ctx_init_fmpq(QQ);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T;
        gr_ctx_struct * top;
        gr_poly_t f, h;
        gr_ptr r;
        int with_pi = n_randint(state, 3) == 0;
        int status;

        gr_tower_init(T, QQ);
        status = _random_tower(T, state, with_pi);
        if (status != GR_SUCCESS)
        {
            gr_tower_clear(T);
            continue;
        }

        top = gr_tower_field(T);
        gr_poly_init(f, top);
        gr_poly_init(h, top);
        GR_TMP_INIT(r, top);

        /* (x - r) h with h of degree 0..2 */
        status = gr_randtest(r, state, top);
        status |= gr_poly_set_coeff_si(f, 1, 1, top);
        status |= gr_neg(gr_poly_coeff_ptr(f, 0, top), r, top);
        status |= gr_poly_randtest(h, state, 1 + n_randint(state, 3), top);
        if (status == GR_SUCCESS && h->length > 0)
        {
            status |= gr_poly_mul(f, f, h, top);

            if (status == GR_SUCCESS && gr_tower_poly_no_roots_modular(f, T, 20))
            {
                flint_printf("FAIL: root excluded\n");
                gr_tower_print(T);
                flint_printf("r = "); gr_println(r, top);
                flint_printf("f = "); gr_poly_print(f, top); flint_printf("\n");
                flint_abort();
            }
        }

        /* x^2 - p for a prime p > 31 not dividing the radicands: no root
           (the tower has degree at most 27 and p does not ramify) */
        status = gr_poly_zero(f, top);
        status |= gr_poly_set_coeff_si(f, 2, 1, top);
        status |= gr_poly_set_coeff_si(f, 0, -(slong) n_nextprime(31 + n_randint(state, 1000), 1), top);
        if (status == GR_SUCCESS)
        {
            total++;
            excluded += gr_tower_poly_no_roots_modular(f, T, 20);
        }

        GR_TMP_CLEAR(r, top);
        gr_poly_clear(f, top);
        gr_poly_clear(h, top);
        gr_tower_clear(T);
    }

    if (excluded < total / 2)
    {
        flint_printf("FAIL: roots excluded in only %wd of %wd cases\n", excluded, total);
        flint_abort();
    }

    /* subset towers: sqrt(2), sqrt(3), sqrt(1 + sqrt(2)); the subset
       {sqrt(2), sqrt(1 + sqrt(2))} has degree 4 */
    {
        gr_tower_t T, S;
        gr_ctx_struct * top;
        gr_ptr x;
        int mark[3] = {1, 0, 1};
        acb_t a, b;

        gr_tower_init(T, QQ);
        gr_tower_init(S, QQ);
        top = gr_tower_field(T);
        GR_TMP_INIT(x, top);
        GR_MUST_SUCCEED(gr_set_ui(x, 2, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_root_ui(T, x, 2, "s2"));
        GR_TMP_CLEAR(x, top);
        top = gr_tower_field(T);
        GR_TMP_INIT(x, top);
        GR_MUST_SUCCEED(gr_set_ui(x, 3, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_root_ui(T, x, 2, "s3"));
        GR_TMP_CLEAR(x, top);
        top = gr_tower_field(T);
        GR_TMP_INIT(x, top);
        /* 1 + sqrt(2), through the generator of step 1 */
        {
            gr_ctx_struct * F1 = gr_tower_field_at(T, 1);
            gr_ptr y;
            GR_TMP_INIT(y, F1);
            GR_MUST_SUCCEED(gr_gen(y, F1));
            GR_MUST_SUCCEED(gr_add_ui(y, y, 1, F1));
            GR_MUST_SUCCEED(gr_tower_promote(x, y, 1, T->length, T));
            GR_TMP_CLEAR(y, F1);
        }
        GR_MUST_SUCCEED(gr_tower_adjoin_root_ui(T, x, 2, "r"));
        GR_TMP_CLEAR(x, top);

        gr_tower_set_subset(S, T, mark);

        acb_init(a);
        acb_init(b);
        if (S->num_gens != 2 || gr_tower_degree(S) != 4 ||
            strcmp(S->gens[0].name, "s2") != 0 || strcmp(S->gens[1].name, "r") != 0 ||
            S->gens[1].def_id != T->gens[2].def_id ||
            !acb_equal(&S->gens[1].enclosure, &T->gens[2].enclosure))
        {
            flint_printf("FAIL: subset tower\n");
            gr_tower_print(T);
            gr_tower_print(S);
            flint_abort();
        }

        /* the generator of S satisfies r^2 = 1 + sqrt(2) */
        {
            gr_ctx_struct * St = gr_tower_field(S);
            gr_ptr g, s2, t;
            GR_TMP_INIT3(g, s2, t, St);
            GR_MUST_SUCCEED(gr_gen(g, St));
            GR_MUST_SUCCEED(gr_sqr(t, g, St));
            GR_MUST_SUCCEED(gr_sub_ui(t, t, 1, St));
            GR_MUST_SUCCEED(gr_sqr(t, t, St));
            GR_MUST_SUCCEED(gr_sub_ui(t, t, 2, St));
            if (gr_is_zero(t, St) != T_TRUE)
            {
                flint_printf("FAIL: subset tower relation\n");
                gr_tower_print(S);
                flint_abort();
            }
            GR_TMP_CLEAR3(g, s2, t, St);
        }

        acb_clear(a);
        acb_clear(b);
        gr_tower_clear(S);
        gr_tower_clear(T);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
