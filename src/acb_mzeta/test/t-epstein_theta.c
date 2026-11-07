/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"
#include "acb.h"
#include "acb_mzeta.h"

/* theta-Mellin against Crandall for random lattices, shifts and nu */
TEST_FUNCTION_START(acb_mzeta_epstein_theta, state)
{
    slong iter;

    for (iter = 0; iter < 3 * flint_test_multiplier(); iter++)
    {
        slong d, i, j, prec;
        arb_mat_t A;
        arb_ptr x, y;
        acb_t nu, z1, z2;
        int flags, ok1, ok2;

        d = 1 + n_randint(state, 3);
        prec = 32 + n_randint(state, 64);
        flags = n_randint(state, 3) == 0 ? ACB_MZETA_EPSTEIN_REGULARIZED : 0;

        arb_mat_init(A, d, d);
        x = _arb_vec_init(d);
        y = _arb_vec_init(d);
        acb_init(nu);
        acb_init(z1);
        acb_init(z2);

        for (i = 0; i < d; i++)
            for (j = 0; j < d; j++)
            {
                if (i == j)
                    arb_set_d(arb_mat_entry(A, i, j), 0.75 + 0.5 * (double) n_randint(state, 1000) / 1000);
                else if (j > i)
                    arb_set_d(arb_mat_entry(A, i, j), ((double) n_randint(state, 1000) - 500) / 1000);
            }
        for (i = 0; i < d; i++)
        {
            if (n_randint(state, 2))
                arb_set_d(x + i, ((double) n_randint(state, 64) - 32) / 64);
            if (n_randint(state, 2))
                arb_set_d(y + i, ((double) n_randint(state, 64) - 32) / 64);
        }
        acb_set_d(nu, ((double) n_randint(state, 64) - 16) / 8);
        if (n_randint(state, 4) == 0)
            arb_set_d(acb_imagref(nu), ((double) n_randint(state, 64) - 32) / 16);

        ok1 = acb_mzeta_epstein_theta(z1, nu, A, x, y, flags, NULL, prec);
        ok2 = acb_mzeta_epstein_crandall(z2, nu, A, x, y, flags, NULL, prec);

        if (ok1 && ok2 && !acb_overlaps(z1, z2))
        {
            flint_printf("FAIL (iter %wd)\n", iter);
            flint_printf("d = %wd, prec = %wd, flags = %d\n", d, prec, flags);
            flint_printf("A = "); arb_mat_printd(A, 10); flint_printf("\n");
            flint_printf("x = "); _arb_vec_printd(x, d, 10); flint_printf("\n");
            flint_printf("y = "); _arb_vec_printd(y, d, 10); flint_printf("\n");
            flint_printf("nu = "); acb_printd(nu, 10); flint_printf("\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n");
            flint_printf("z2 = "); acb_printd(z2, 30); flint_printf("\n");
            flint_abort();
        }

        if (ok1 && acb_is_finite(z2) && acb_rel_accuracy_bits(z2) > prec - 8
            && acb_rel_accuracy_bits(z1) < prec - 16)
        {
            flint_printf("FAIL (accuracy)\n");
            flint_printf("d = %wd, prec = %wd, flags = %d\n", d, prec, flags);
            flint_printf("nu = "); acb_printd(nu, 10); flint_printf("\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n");
            flint_printf("z2 = "); acb_printd(z2, 30); flint_printf("\n");
            flint_abort();
        }

        arb_mat_clear(A);
        _arb_vec_clear(x, d);
        _arb_vec_clear(y, d);
        acb_clear(nu);
        acb_clear(z1);
        acb_clear(z2);
    }

    TEST_FUNCTION_END(state);
}
