/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "acb.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_ql_finite_diff, state)
{
    slong iter;

    /* Test: find correct coefficients for exp function */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong prec = 100 + n_randint(state, 1000);
        slong ord = n_randint(state, 4);
        slong g = 1 + n_randint(state, 4);
        slong b = ord + 1;
        slong nb_val = n_pow(b, g);
        slong nb_fd = acb_theta_jet_nb(ord, g);
        slong * tups;
        arb_t c, rho;
        arf_t eps, err;
        acb_ptr val, df, test;
        acb_t x, t;
        fmpz_t m;
        slong k, kk, j, i;

        tups = flint_malloc(g * nb_fd * sizeof(slong));
        arb_init(c);
        arb_init(rho);
        arf_init(eps);
        arf_init(err);
        val = _acb_vec_init(nb_val);
        df = _acb_vec_init(nb_fd);
        test = _acb_vec_init(nb_fd);
        acb_init(x);
        acb_init(t);
        fmpz_init(m);

        /* Get c, rho, eps, err */
        arb_one(rho);
        arb_set_si(c, g);
        arb_exp(c, c, prec);
        acb_theta_jet_ql_radius(eps, err, c, rho, ord, g, prec);

        /* Fill in values, apply jet_fd at 2*prec */
        for (k = 0; k < nb_val; k++)
        {
            acb_zero(x);
            kk = k;
            for (j = 0; j < g; j++)
            {
                acb_unit_root(t, b, 2 * prec);
                acb_pow_ui(t, t, (kk % b), 2 * prec);
                acb_add(x, x, t, 2 * prec);
                kk = kk / b;
            }
            acb_zero(t);
            arb_set_arf(acb_realref(t), eps);
            acb_mul(x, x, t, 2 * prec);
            acb_exp(&val[k], x, 2 * prec);
        }
        acb_theta_jet_ql_finite_diff(df, eps, err, rho, val, ord, g, 2 * prec);

        /* Fill in test */
        acb_theta_jet_tuples(tups, ord, g);
        for (j = 0; j < nb_fd; j++)
        {
            acb_one(x);
            for (i = 0; i < g; i++)
            {
                fmpz_fac_ui(m, tups[j * g + i]);
                acb_div_fmpz(x, x, m, prec);
            }
            acb_set(&test[j], x);
        }

        if (!_acb_vec_overlaps(df, test, nb_fd))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, ord = %wd, values:\n", g, ord);
            _acb_vec_printd(val, nb_val, 5);
            flint_printf("taylor coeffs:\n");
            _acb_vec_printd(df, nb_fd, 5);
            flint_printf("test:\n");
            _acb_vec_printd(test, nb_fd, 5);
            flint_printf("difference:\n");
            _acb_vec_sub(test, test, df, nb_fd, prec);
            _acb_vec_printd(test, nb_fd, 5);
            flint_abort();
        }

        flint_free(tups);
        arb_clear(c);
        arb_clear(rho);
        arf_clear(eps);
        arf_clear(err);
        _acb_vec_clear(val, nb_val);
        _acb_vec_clear(df, nb_fd);
        _acb_vec_clear(test, nb_fd);
        acb_clear(x);
        acb_clear(t);
        fmpz_clear(m);
    }

    TEST_FUNCTION_END(state);
}
