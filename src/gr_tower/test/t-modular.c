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
#include "fmpz_poly.h"
#include "gr_poly.h"
#include "acb.h"
#include "qqbar.h"
#include "ulong_extras.h"
#include "gr_tower.h"

/*
    Modular irreducibility proofs: reducible moduli are never "proven"
    (X^n - b^n over towers containing b, X^2 - 2 over Q(sqrt 2), and
    minimal polynomials over Q of elements of the tower), while the
    genuinely irreducible moduli of random radical towers are proven
    most of the time.
*/

/* adjoins the principal n-th root of x without proof (as a dynamic step) */
static int
_adjoin_root_dynamic(gr_tower_t T, gr_srcptr x, ulong n)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_poly_t m;
    acb_t z;
    int status;

    gr_poly_init(m, top);
    acb_init(z);
    status = gr_poly_set_coeff_si(m, n, 1, top);
    status |= gr_neg(gr_poly_coeff_ptr(m, 0, top), x, top);
    _gr_poly_normalise(m, top);
    status |= gr_tower_get_acb(z, x, GR_TOWER_DEFAULT_PREC, T);
    acb_root_ui(z, z, n, GR_TOWER_DEFAULT_PREC);
    if (status == GR_SUCCESS)
        status = gr_tower_adjoin_algebraic(T, m, z, GR_TOWER_STATUS_DYNAMIC, NULL);
    gr_poly_clear(m, top);
    acb_clear(z);
    return status;
}

TEST_FUNCTION_START(gr_tower_modular, state)
{
    gr_ctx_t QQ;
    slong iter, proven = 0, total = 0;

    gr_ctx_init_fmpq(QQ);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T;
        slong steps = 1 + n_randint(state, 3), k;
        int status = GR_SUCCESS;

        gr_tower_init(T, QQ);

        /* a tower of radicals of small elements: sqrt(p), cbrt(1 + a) ...
           adjoined through gr_tower_adjoin_root_ui, which proves the
           steps by Capelli's theorem when it can */
        for (k = 0; k < steps && status == GR_SUCCESS; k++)
        {
            gr_ctx_struct * top = gr_tower_field(T);
            gr_ptr x;
            ulong n = 2 + n_randint(state, 3);

            GR_TMP_INIT(x, top);
            if (T->length == 0 || n_randint(state, 2))
                status |= gr_set_ui(x, 2 + n_randint(state, 20), top);
            else
            {
                /* 1 + (the last generator) + small integer */
                status |= gr_gen(x, top);
                status |= gr_add_ui(x, x, 1 + n_randint(state, 5), top);
            }
            status |= gr_tower_adjoin_root_ui(T, x, n, NULL);
            GR_TMP_CLEAR(x, top);
        }

        if (status != GR_SUCCESS)
        {
            gr_tower_clear(T);
            continue;
        }

        for (k = 1; k <= T->length; k++)
        {
            total++;
            proven += (GR_TOWER_STEP(T, k - 1)->status == GR_TOWER_STATUS_PROVEN);
        }

        /* reducible binomials: X^n - b^n for an element b of the tower */
        {
            gr_ctx_struct * top = gr_tower_field(T);
            gr_ptr b, bn;
            ulong n = 2 + n_randint(state, 4);
            slong len = T->length;

            GR_TMP_INIT2(b, bn, top);
            status = gr_randtest_not_zero(b, state, top);
            if (n_randint(state, 2))
                status |= gr_add_ui(b, b, 1, top);
            status |= gr_pow_ui(bn, b, n, top);

            if (status == GR_SUCCESS)
            {
                status = _adjoin_root_dynamic(T, bn, n);
                if (status == GR_SUCCESS)
                {
                    if (gr_tower_prove_step_modular(T, len + 1, 20))
                    {
                        flint_printf("FAIL: X^n - b^n proven irreducible\n");
                        gr_tower_print(T);
                        flint_printf("b = "); gr_println(b, top);
                        flint_abort();
                    }
                }
            }
            GR_TMP_CLEAR2(b, bn, top);
        }

        gr_tower_clear(T);
    }

    /* X^2 - 2 over Q(sqrt 2), and the minimal polynomial of sqrt 2 + 1
       over Q adjoined to Q(sqrt 2): reducible */
    {
        gr_tower_t T;
        gr_ctx_struct * top;
        gr_ptr x;
        qqbar_t q;

        gr_tower_init(T, QQ);
        top = gr_tower_field(T);
        GR_TMP_INIT(x, top);
        GR_MUST_SUCCEED(gr_set_ui(x, 2, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
        GR_TMP_CLEAR(x, top);

        top = gr_tower_field(T);
        GR_TMP_INIT(x, top);
        GR_MUST_SUCCEED(gr_set_ui(x, 2, top));
        GR_MUST_SUCCEED(_adjoin_root_dynamic(T, x, 2));
        GR_TMP_CLEAR(x, top);

        if (gr_tower_prove_step_modular(T, 2, 50))
        {
            flint_printf("FAIL: X^2 - 2 over Q(sqrt 2) proven irreducible\n");
            flint_abort();
        }

        qqbar_init(q);
        qqbar_set_ui(q, 2);
        qqbar_sqrt(q, q);
        qqbar_add_ui(q, q, 1);
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, q, NULL));
        if (GR_TOWER_STEP(T, 2)->status == GR_TOWER_STATUS_PROVEN || gr_tower_prove_step_modular(T, 3, 50))
        {
            flint_printf("FAIL: minimal polynomial of sqrt 2 + 1 proven irreducible over Q(sqrt 2)\n");
            flint_abort();
        }
        qqbar_clear(q);

        gr_tower_clear(T);
    }

    /* X^4 + 4 = (X^2 + 2X + 2)(X^2 - 2X + 2) over Q */
    {
        gr_tower_t T;
        gr_ctx_struct * top;
        gr_ptr x;

        gr_tower_init(T, QQ);
        top = gr_tower_field(T);
        GR_TMP_INIT(x, top);
        GR_MUST_SUCCEED(gr_set_si(x, -4, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_root_ui(T, x, 4, NULL));
        if (GR_TOWER_STEP(T, 0)->status == GR_TOWER_STATUS_PROVEN || gr_tower_prove_step_modular(T, 1, 50))
        {
            flint_printf("FAIL: X^4 + 4 proven irreducible\n");
            flint_abort();
        }
        GR_TMP_CLEAR(x, top);
        gr_tower_clear(T);
    }

    /* Trager's method: X^4 + 1 over Q(sqrt 2) is refined to a quadratic
       step (X^4 + 1 = (X^2 - sqrt2 X + 1)(X^2 + sqrt2 X + 1)), X^4 + 1
       over Q(sqrt 3) is proven irreducible, X^5 - X - 1 over
       Q(sqrt 2, sqrt 3) is proven irreducible, and X^2 - 6 over
       Q(sqrt 2, sqrt 3) is refined to a linear step */
    {
        gr_tower_t T;
        gr_ctx_struct * top;
        gr_ptr x;
        slong i;

        for (i = 0; i < 4; i++)
        {
            gr_tower_init(T, QQ);
            top = gr_tower_field(T);
            GR_TMP_INIT(x, top);
            GR_MUST_SUCCEED(gr_set_ui(x, (i == 1) ? 3 : 2, top));
            GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
            GR_TMP_CLEAR(x, top);

            if (i >= 2)
            {
                top = gr_tower_field(T);
                GR_TMP_INIT(x, top);
                GR_MUST_SUCCEED(gr_set_ui(x, 3, top));
                GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
                GR_TMP_CLEAR(x, top);
            }

            if (i == 3)
            {
                top = gr_tower_field(T);
                GR_TMP_INIT(x, top);
                GR_MUST_SUCCEED(gr_set_ui(x, 6, top));
                GR_MUST_SUCCEED(_adjoin_root_dynamic(T, x, 2));
                GR_TMP_CLEAR(x, top);
            }
            else
            {
                gr_poly_t m;
                acb_t z;

                top = gr_tower_field(T);
                gr_poly_init(m, top);
                acb_init(z);
                if (i == 2)
                {
                    /* X^5 - X - 1, root near 1.1673 */
                    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 5, 1, top));
                    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 1, -1, top));
                    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, -1, top));
                    acb_set_d(z, 1.1673039782614187);
                }
                else
                {
                    /* X^4 + 1, root exp(i pi/4) */
                    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 4, 1, top));
                    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, 1, top));
                    acb_set_d_d(z, 0.7071067811865476, 0.7071067811865476);
                }
                GR_MUST_SUCCEED(gr_tower_adjoin_algebraic(T, m, z, GR_TOWER_STATUS_DYNAMIC, NULL));
                gr_poly_clear(m, top);
                acb_clear(z);
            }

            if (GR_TOWER_STEP(T, T->length - 1)->status == GR_TOWER_STATUS_PROVEN)
            {
                flint_printf("FAIL: Trager case %wd: step proven before proving\n", i);
                flint_abort();
            }

            if (!gr_tower_prove_step_trager(T, T->length, GR_TOWER_OPTION(T, GR_TOWER_OPT_TRAGER_DEGREE_LIMIT)))
            {
                flint_printf("FAIL: Trager case %wd: not decided\n", i);
                gr_tower_print(T);
                flint_abort();
            }

            if (GR_TOWER_STEP(T, T->length - 1)->status != GR_TOWER_STATUS_PROVEN ||
                gr_tower_degree(T) != ((i == 0) ? 4 : (i == 1) ? 8 : (i == 2) ? 20 : 4))
            {
                flint_printf("FAIL: Trager case %wd: status %d, degree %wd\n", i,
                    GR_TOWER_STEP(T, T->length - 1)->status, gr_tower_degree(T));
                gr_tower_print(T);
                flint_abort();
            }

            gr_tower_clear(T);
        }
    }

    /* the irreducible steps of the random towers are proven most of the
       time (radicands which are p-th powers happen: 4^(1/2), 8^(1/3), or
       the square root of a square in the tower) */
    if (proven < total / 2)
    {
        flint_printf("FAIL: only %wd of %wd radical steps proven\n", proven, total);
        flint_abort();
    }

    /* Eisenstein in pi over an unproven step: a = root of x^2 - 1 near 1
       (dynamic), then X^2 - ((a - 1) + pi) pi, which looks Eisenstein in pi
       over Q(pi)[a] but is X^2 - pi^2 = (X - pi)(X + pi) in the field the
       tower refines to; the step must not be proven */
    {
        gr_tower_t T;
        gr_ctx_struct * K;
        gr_poly_t m;
        gr_ptr a;
        acb_t z;

        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
        K = gr_tower_field(T);
        gr_poly_init(m, K);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 2, 1, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, -1, K));
        acb_init(z);
        acb_one(z);
        mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -20);
        GR_MUST_SUCCEED(gr_tower_adjoin_algebraic(T, m, z, GR_TOWER_STATUS_DYNAMIC, "a"));
        gr_poly_clear(m, K);
        K = gr_tower_field(T);
        gr_poly_init(m, K);
        a = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_gen(a, K));
        GR_MUST_SUCCEED(gr_sub_ui(a, a, 1, K));
        {
            gr_ptr pi = gr_heap_init(K);
            GR_MUST_SUCCEED(gr_tower_gen_get(pi, T, 0));
            GR_MUST_SUCCEED(gr_add(a, a, pi, K));
            GR_MUST_SUCCEED(gr_mul(a, a, pi, K));
            gr_heap_clear(pi, K);
        }
        GR_MUST_SUCCEED(gr_neg(a, a, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 2, 1, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(m, 0, a, K));
        acb_const_pi(z, 64);
        mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -8);
        mag_set_ui_2exp_si(arb_radref(acb_imagref(z)), 1, -8);
        GR_MUST_SUCCEED(gr_tower_adjoin_algebraic(T, m, z, GR_TOWER_STATUS_DYNAMIC, "b"));
        if (gr_tower_prove_step_modular(T, 2, 6) || GR_TOWER_STEP(T, 1)->status == GR_TOWER_STATUS_PROVEN)
        {
            flint_printf("FAIL: X^2 - (a - 1) pi proven over a dynamic step\n");
            flint_abort();
        }
        gr_heap_clear(a, K);
        gr_poly_clear(m, K);
        acb_clear(z);
        gr_tower_clear(T);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
