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

TEST_FUNCTION_START(acb_dirichlet_zeta_zeros, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t v1, v2, z1, z2;
        fmpz_t n, m, k;
        acb_ptr p;
        slong prec1, prec2;
        slong i, len;
        const slong maxlen = 5;

        acb_init(v1);
        acb_init(v2);
        acb_init(z1);
        acb_init(z2);
        fmpz_init(n);
        fmpz_init(m);
        fmpz_init(k);
        p = _acb_vec_init(maxlen);

        fmpz_randtest_unsigned(n, state, 20);
        fmpz_add_ui(n, n, 1);
        prec1 = 2 + n_randtest(state) % 50;
        prec2 = 2 + n_randtest(state) % 200;

        len = 1 + n_randint(state, maxlen);
        i = n_randint(state, len);
        acb_dirichlet_zeta_zeros(p, n, len, prec1);
        acb_set(z1, p + i);

        fmpz_add_si(k, n, i);
        acb_dirichlet_zeta_zero(z2, k, prec2);

        acb_dirichlet_zeta(v1, z1, prec1 + 20);
        acb_dirichlet_zeta(v2, z2, prec2 + 20);

        if (!acb_overlaps(z1, z2) || !acb_contains_zero(v1) || !acb_contains_zero(v2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = "); fmpz_print(n);
            flint_printf("   prec1 = %wd  prec2 = %wd\n\n", prec1, prec2);
            flint_printf("z1 = "); acb_printn(z1, 100, 0); flint_printf("\n\n");
            flint_printf("z2 = "); acb_printn(z2, 100, 0); flint_printf("\n\n");
            flint_printf("v1 = "); acb_printn(v1, 100, 0); flint_printf("\n\n");
            flint_printf("v2 = "); acb_printn(v2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (acb_rel_accuracy_bits(z1) < prec1 - 3 || acb_rel_accuracy_bits(z2) < prec2 - 3)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("n = "); fmpz_print(n);
            flint_printf("   prec1 = %wd  prec2 = %wd\n\n", prec1, prec2);
            flint_printf("acc(z1) = %wd, acc(z2) = %wd\n\n", acb_rel_accuracy_bits(z1), acb_rel_accuracy_bits(z2));
            flint_printf("z1 = "); acb_printn(z1, 100, 0); flint_printf("\n\n");
            flint_printf("z2 = "); acb_printn(z2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z1);
        acb_clear(z2);
        acb_clear(v1);
        acb_clear(v2);
        fmpz_clear(n);
        fmpz_clear(m);
        fmpz_clear(k);
        _acb_vec_clear(p, maxlen);
    }

    TEST_FUNCTION_END(state);
}
