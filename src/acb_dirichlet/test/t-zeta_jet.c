/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_zeta_jet, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, one;
        acb_ptr u, v;
        slong prec, len, j;
        int deflate;

        prec = 2 + n_randint(state, n_randint(state, 4) ? 200 : 1500);
        len = 1 + n_randint(state, 8);
        deflate = n_randint(state, 2);

        acb_init(s);
        acb_init(one);
        u = _acb_vec_init(len);
        v = _acb_vec_init(len);
        acb_one(one);

        /* mostly around the cutoffs for the reflection formula */
        switch (n_randint(state, 4))
        {
            case 0:
                acb_set_si(s, -(slong) n_randint(state, 8));
                break;
            case 1:
                arb_set_si(acb_realref(s), (slong) n_randint(state, 64) - 48);
                arb_mul_2exp_si(acb_realref(s), acb_realref(s), -3);
                break;
            case 2:
                arb_set_si(acb_realref(s), (slong) n_randint(state, 64) - 48);
                arb_mul_2exp_si(acb_realref(s), acb_realref(s), -3);
                arb_set_si(acb_imagref(s), (slong) n_randint(state, 40) - 20);
                break;
            default:
                acb_randtest(s, state, 1 + n_randint(state, 200), 3);
        }

        if (!deflate && arb_contains_si(acb_realref(s), 1) && arb_contains_zero(acb_imagref(s)))
            deflate = 1;

        acb_dirichlet_zeta_jet(u, s, deflate, len, prec);

        /* Euler-Maclaurin summation directly at s, with plenty of
           guard bits */
        _acb_poly_zeta_cpx_series(v, s, one, deflate, len, prec + 200);

        for (j = 0; j < len; j++)
        {
            if (!acb_overlaps(u + j, v + j))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd  prec = %wd  len = %wd  deflate = %d  j = %wd\n\n", iter, prec, len, deflate, j);
                flint_printf("s = "); acb_printn(s, 30, 0); flint_printf("\n\n");
                flint_printf("u = "); acb_printn(u + j, 30, 0); flint_printf("\n\n");
                flint_printf("v = "); acb_printn(v + j, 30, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        /* the trivial zeros are exact */
        if (!deflate && acb_is_int(s) && arf_sgn(arb_midref(acb_realref(s))) < 0 &&
            arf_is_int_2exp_si(arb_midref(acb_realref(s)), 1) && !acb_is_zero(u))
        {
            flint_printf("FAIL: trivial zero\n\n");
            flint_printf("s = "); acb_printn(s, 30, 0); flint_printf("\n\n");
            flint_printf("u = "); acb_printn(u, 30, 0); flint_printf("\n\n");
            flint_abort();
        }

        /* accuracy for small exact s, where Euler-Maclaurin summation is
           used directly or after the reflection formula */
        if (prec >= 32 && acb_is_exact(s) && acb_is_real(s) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 3) < 0)
        {
            for (j = 0; j < len; j++)
            {
                if (!acb_is_zero(v + j) && !acb_contains_zero(v + j) &&
                    acb_rel_accuracy_bits(u + j) < prec - 3 * len - 12)
                {
                    flint_printf("FAIL: accuracy\n\n");
                    flint_printf("iter = %wd  prec = %wd  len = %wd  deflate = %d  j = %wd\n\n", iter, prec, len, deflate, j);
                    flint_printf("s = "); acb_printn(s, 30, 0); flint_printf("\n\n");
                    flint_printf("u = "); acb_printn(u + j, 30, 0); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        acb_clear(s);
        acb_clear(one);
        _acb_vec_clear(u, len);
        _acb_vec_clear(v, len);
    }

    TEST_FUNCTION_END(state);
}
