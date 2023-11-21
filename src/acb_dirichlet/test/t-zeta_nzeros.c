/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_zeta_nzeros, state)
{
    slong iter;

    for (iter = 0; iter < 130 + 20 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, t;
        arf_t f1, f2;
        fmpz_t n, m, k1, k2;
        acb_t z;
        slong prec1, prec2;

        arb_init(x);
        arb_init(y);
        arb_init(t);
        acb_init(z);
        arf_init(f1);
        arf_init(f2);
        fmpz_init(n);
        fmpz_init(m);
        fmpz_init(k1);
        fmpz_init(k2);

        if (iter < 130)
        {
            fmpz_set_si(n, iter + 1);
        }
        else
        {
            fmpz_randtest_unsigned(n, state, 20);
            fmpz_add_ui(n, n, 131);
        }

        prec1 = 2 + n_randtest(state) % 50;
        prec2 = 2 + n_randtest(state) % 200;

        acb_dirichlet_hardy_z_zero(t, n, prec1);

        acb_dirichlet_zeta_nzeros(x, t, prec2);

        fmpz_sub_ui(m, n, 1);
        if (!arb_contains_fmpz(x, n) || !arb_contains_fmpz(x, m))
        {
            flint_printf("FAIL: zero containment\n\n");
            flint_printf("n = "); fmpz_print(n);
            flint_printf("   prec1 = %wd  prec2 = %wd\n\n", prec1, prec2);
            flint_printf("t = "); arb_printn(t, 100, 0); flint_printf("\n\n");
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_get_lbound_arf(f1, t, prec1);
        arb_get_ubound_arf(f2, t, prec1);

        _acb_dirichlet_exact_zeta_nzeros(k1, f1);
        _acb_dirichlet_exact_zeta_nzeros(k2, f2);

        if (!(fmpz_cmp(k1, n) < 0 && fmpz_cmp(k2, n) >= 0))
        {
            flint_printf("FAIL: near a zero\n\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("t = "); arb_printn(t, 100, 0); flint_printf("\n\n");
            flint_printf("f1 = "); arf_print(f1); flint_printf("\n\n");
            flint_printf("f2 = "); arf_print(f2); flint_printf("\n\n");
            flint_printf("k1 = "); fmpz_print(k1); flint_printf("\n\n");
            flint_printf("k2 = "); fmpz_print(k2); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_sub_ui(m, n, 2);
        acb_dirichlet_gram_point(t, m, NULL, NULL, prec1);
        arb_get_lbound_arf(f1, t, prec1);
        arb_get_ubound_arf(f2, t, prec1);
        arb_set_interval_arf(t, f1, f2, prec1);
        acb_set_arb(z, t);
        acb_dirichlet_hardy_z(z, z, NULL, NULL, 1, prec1);
        if (!acb_contains_zero(z))
        {
            _acb_dirichlet_exact_zeta_nzeros(k1, f1);
            _acb_dirichlet_exact_zeta_nzeros(k2, f2);
            if (!fmpz_equal(k1, k2))
            {
                flint_printf("FAIL: near a gram point\n\n");
                flint_printf("m = "); fmpz_print(m); flint_printf("\n\n");
                flint_printf("t = "); arb_printn(t, 100, 0); flint_printf("\n\n");
                flint_printf("f1 = "); arf_print(f1); flint_printf("\n\n");
                flint_printf("f2 = "); arf_print(f2); flint_printf("\n\n");
                flint_printf("k1 = "); fmpz_print(k1); flint_printf("\n\n");
                flint_printf("k2 = "); fmpz_print(k2); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(t);
        acb_clear(z);
        arf_clear(f1);
        arf_clear(f2);
        fmpz_clear(n);
        fmpz_clear(m);
        fmpz_clear(k1);
        fmpz_clear(k2);
    }

    TEST_FUNCTION_END(state);
}
