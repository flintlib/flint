/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"
#include "dfloat.h"

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

    /* the choice of the large height method and its working precision */
    {
        fmpz_t n;
        int ok = 1;
        fmpz_init(n);
        fmpz_ui_pow_ui(n, 10, 13);      /* L = log2(n) = 44 */
        ok = ok && _acb_dirichlet_hardy_z_zeros_use_platt(n, 200, 0) != 0;
        ok = ok && !_acb_dirichlet_hardy_z_zeros_use_platt(n, 5, 0);
        if (dfloat_is_supported())
        {
            ok = ok && _acb_dirichlet_hardy_z_zeros_use_platt(n, 14, 108) == 148;
            ok = ok && !_acb_dirichlet_hardy_z_zeros_use_platt(n, 13, 108);
            /* the arb sums (up to the accuracy of the method), else
               the dfloat sums and refinement */
            ok = ok && _acb_dirichlet_hardy_z_zeros_use_platt(n, 200, 180) == 228;
            ok = ok && _acb_dirichlet_hardy_z_zeros_use_platt(n, 200, 500) == 148;
            ok = ok && _acb_dirichlet_hardy_z_zeros_use_platt(n, 200, 30) == 100;
        }
        fmpz_sub_ui(n, n, 1);
        ok = ok && !_acb_dirichlet_hardy_z_zeros_use_platt(n, 15, 0);
        fmpz_ui_pow_ui(n, 10, 11);
        fmpz_sub_ui(n, n, 1);
        ok = ok && !_acb_dirichlet_hardy_z_zeros_use_platt(n, 1000, 0);
        fmpz_ui_pow_ui(n, 10, 20);
        ok = ok && _acb_dirichlet_hardy_z_zeros_use_platt(n, 1, 0) != 0;
        fmpz_ui_pow_ui(n, 10, 23);
        ok = ok && !_acb_dirichlet_hardy_z_zeros_use_platt(n, 1000, 0);
        if (!ok)
        {
            flint_printf("FAIL: choice of the method\n");
            flint_abort();
        }
        fmpz_clear(n);
    }

    /* the refinement of a zero from a ball (as for the zeros of the large
       height method beyond its accuracy) against Riemann-Siegel */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        fmpz_t n;
        arb_t z, w, v;
        slong prec0, prec;

        fmpz_init(n);
        arb_init(z);
        arb_init(w);
        arb_init(v);
        fmpz_randtest_unsigned(n, state, 30);
        fmpz_add_ui(n, n, 1);
        /* (a ball isolating the zero: well below the spacing) */
        prec0 = fmpz_bits(n) + 16 + n_randint(state, 100);
        prec = prec0 + n_randint(state, 300);

        acb_dirichlet_hardy_z_zero(z, n, prec0);
        _acb_dirichlet_refine_hardy_z_zero_ball(w, z, prec);
        acb_dirichlet_hardy_z_zero(v, n, prec);

        if (!arb_overlaps(w, v) || !arb_overlaps(w, z) ||
            arb_rel_accuracy_bits(w) < prec - 3)
        {
            flint_printf("FAIL: refinement from a ball\n\n");
            flint_printf("n = "); fmpz_print(n);
            flint_printf("  prec0 = %wd  prec = %wd\n\n", prec0, prec);
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_printf("w = "); arb_printn(w, 100, 0); flint_printf("\n\n");
            flint_printf("v = "); arb_printn(v, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(n);
        arb_clear(z);
        arb_clear(w);
        arb_clear(v);
    }

    /* the large height method (one run, a few seconds; only with dfloat,
       for which the choice is tuned) against Riemann-Siegel */
    if (dfloat_is_supported())
    {
        fmpz_t n, k;
        acb_ptr p;
        acb_t z;
        slong len = 15, prec, i, j;

        fmpz_init(n);
        fmpz_init(k);
        acb_init(z);
        p = _acb_vec_init(len);
        fmpz_ui_pow_ui(n, 10, 13);
        fmpz_add_ui(n, n, n_randint(state, 1000000));
        prec = 60 + n_randint(state, 48);

        if (!_acb_dirichlet_hardy_z_zeros_use_platt(n, len, prec))
        {
            flint_printf("FAIL: platt not chosen\n");
            flint_abort();
        }

        acb_dirichlet_zeta_zeros(p, n, len, prec);

        for (j = 0; j < 2; j++)
        {
            i = (j == 0) ? 0 : 1 + n_randint(state, len - 1);
            fmpz_add_si(k, n, i);
            acb_dirichlet_zeta_zero(z, k, prec);
            if (!acb_overlaps(z, p + i) || acb_rel_accuracy_bits(p + i) < prec - 3)
            {
                flint_printf("FAIL: large height method\n\n");
                flint_printf("n = "); fmpz_print(n); flint_printf("  i = %wd  prec = %wd\n\n", i, prec);
                flint_printf("z = "); acb_printn(z, 50, 0); flint_printf("\n\n");
                flint_printf("p = "); acb_printn(p + i, 50, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        fmpz_clear(n);
        fmpz_clear(k);
        acb_clear(z);
        _acb_vec_clear(p, len);
    }

    TEST_FUNCTION_END(state);
}
