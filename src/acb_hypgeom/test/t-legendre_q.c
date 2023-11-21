/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

/* these functions are not public for now */
void _acb_hypgeom_legendre_q_single(acb_t res, const acb_t n, const acb_t m,
    const acb_t z, slong prec);
void _acb_hypgeom_legendre_q_double(acb_t res, const acb_t n, const acb_t m,
    const acb_t z, slong prec);

TEST_FUNCTION_START(acb_hypgeom_legendre_q, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t n, m, z, res1, res2;
        slong prec1, prec2, ebits;

        acb_init(n);
        acb_init(m);
        acb_init(z);
        acb_init(res1);
        acb_init(res2);

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);
        ebits = 1 + n_randint(state, 10);

        if (n_randint(state, 2))
        {
            acb_set_si(m, n_randint(state, 20) - 10);
            acb_set_si(n, n_randint(state, 20) - 10);
        }
        else
        {
            acb_randtest_param(n, state, 1 + n_randint(state, 400), ebits);
            acb_randtest_param(m, state, 1 + n_randint(state, 400), ebits);
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 400), ebits);

        _acb_hypgeom_legendre_q_single(res1, n, m, z, prec1);
        _acb_hypgeom_legendre_q_double(res2, n, m, z, prec2);

        if (!acb_overlaps(res1, res2))
        {
            flint_printf("FAIL: consistency 1\n\n");
            flint_printf("iter = %wd, prec1 = %wd, prec2 = %wd\n\n", iter, prec1, prec2);
            flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
            flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(n);
        acb_clear(m);
        acb_clear(z);
        acb_clear(res1);
        acb_clear(res2);
    }

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t n, m, z, res1, res2, t, u;
        slong prec1, prec2, ebits;
        int type;

        acb_init(n);
        acb_init(m);
        acb_init(z);
        acb_init(res1);
        acb_init(res2);
        acb_init(t);
        acb_init(u);

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);
        ebits = 1 + n_randint(state, 10);

        if (n_randint(state, 2))
        {
            acb_set_si(m, n_randint(state, 20) - 10);
            acb_set_si(n, n_randint(state, 20) - 10);
        }
        else
        {
            acb_randtest_param(n, state, 1 + n_randint(state, 400), ebits);
            acb_randtest_param(m, state, 1 + n_randint(state, 400), ebits);
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 400), ebits);

        type = n_randint(state, 2);

        acb_hypgeom_legendre_q(res1, n, m, z, type, prec1);

        acb_neg(t, m);
        acb_hypgeom_legendre_p(res2, n, t, z, type, prec2);
        acb_add(u, m, n, prec2);
        acb_add_ui(u, u, 1, prec2);
        acb_gamma(u, u, prec2);
        acb_mul(res2, res2, u, prec2);
        acb_sub(u, n, m, prec2);
        acb_add_ui(u, u, 1, prec2);
        acb_rgamma(u, u, prec2);
        acb_mul(res2, res2, u, prec2);

        acb_hypgeom_legendre_p(t, n, m, z, type, prec2);

        if (type == 0)
        {
            acb_cos_pi(u, m, prec2);
            acb_mul(t, t, u, prec2);
        }
        acb_sub(res2, t, res2, prec2);

        if (type == 1)
        {
            acb_exp_pi_i(t, m, prec2);
            acb_mul(res2, res2, t, prec2);
        }

        acb_sin_pi(t, m, prec2);

        if (acb_contains_zero(t))
            acb_indeterminate(res2);
        else
            acb_div(res2, res2, t, prec2);

        acb_const_pi(t, prec2);
        acb_mul(res2, res2, t, prec2);
        acb_mul_2exp_si(res2, res2, -1);

        if (!acb_overlaps(res1, res2))
        {
            flint_printf("FAIL: consistency 2\n\n");
            flint_printf("iter = %wd, prec1 = %wd, prec2 = %wd\n\n", iter, prec1, prec2);
            flint_printf("type = %d\n\n", type);
            flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
            flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(n);
        acb_clear(m);
        acb_clear(z);
        acb_clear(res1);
        acb_clear(res2);
        acb_clear(t);
        acb_clear(u);
    }

    TEST_FUNCTION_END(state);
}
