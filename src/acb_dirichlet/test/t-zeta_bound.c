/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_zeta_bound, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, sb, z;
        arb_t x, y;
        slong prec, i;
        mag_t zbound;

        acb_init(s);
        acb_init(sb);
        acb_init(z);
        arb_init(x);
        arb_init(y);
        mag_init(zbound);

        acb_randtest(s, state, 2 + n_randint(state, 100), 4);
        if (n_randint(state, 2))
            arb_set_d(acb_realref(s), -0.5 + n_randint(state, 32) / 16.0);
        mag_zero(arb_radref(acb_realref(s)));
        mag_zero(arb_radref(acb_imagref(s)));

        /* create a larger interval containing s */
        acb_set(sb, s);
        for (i = 0; i < 2; i++)
        {
            arb_randtest(x, state, 2 + n_randint(state, 100), 2);
            arb_add(x, x, acb_realref(s), 2 * MAG_BITS);
            arb_intersection(acb_realref(sb), acb_realref(sb), x, 2 * MAG_BITS);

            arb_randtest(x, state, 2 + n_randint(state, 100), 2);
            arb_add(x, x, acb_imagref(s), 2 * MAG_BITS);
            arb_intersection(acb_imagref(sb), acb_imagref(sb), x, 2 * MAG_BITS);
        }

        acb_dirichlet_zeta_bound(zbound, sb);

        if (!mag_is_inf(zbound))
        {
            for (prec = 64; ; prec *= 2)
            {
                acb_zeta(z, s, prec);
                if (acb_rel_accuracy_bits(z) > MAG_BITS)
                    break;
            }

            acb_abs(x, z, 2 * MAG_BITS);

            arb_zero(y);
            arf_set_mag(arb_midref(y), zbound);

            if (!arb_le(x, y))
            {
                flint_printf("FAIL: bound\n\n");
                flint_printf("iter = %wd\n", iter);
                flint_printf("s = "); acb_printn(s, 50, 0); flint_printf("\n\n");
                flint_printf("sb = "); acb_printn(sb, 50, 0); flint_printf("\n\n");
                flint_printf("z = "); acb_printn(z, 50, 0); flint_printf("\n\n");
                flint_printf("zbound = "); mag_printd(zbound, 10); flint_printf("\n\n");
                flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
                flint_printf("y = "); arb_printn(y, 50, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(s);
        acb_clear(sb);
        acb_clear(z);
        arb_clear(x);
        arb_clear(y);
        mag_clear(zbound);
    }

    /* test deriv_bound */
    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s;
        acb_ptr v;
        mag_t b, b1, b2, m, m1, m2;

        acb_init(s);
        v = _acb_vec_init(3);
        mag_init(b);
        mag_init(b1);
        mag_init(b2);
        mag_init(m);
        mag_init(m1);
        mag_init(m2);

        acb_randtest(s, state, 2 + n_randint(state, 100), 2);
        arb_mul_ui(acb_realref(s), acb_realref(s), n_randtest(state) % 100, 100);
        arb_mul_ui(acb_imagref(s), acb_imagref(s), n_randtest(state) % 10000, 100);

        acb_dirichlet_zeta_bound(b, s);
        acb_dirichlet_zeta_deriv_bound(b1, b2, s);

        acb_get_mid(s, s);

        acb_dirichlet_zeta_jet(v, s, 0, 3, 53);

        acb_get_mag(m, v);
        acb_get_mag(m1, v + 1);
        acb_get_mag(m2, v + 2);

        if (mag_cmp(m, b) > 0 || mag_cmp(m1, b1) > 0 || mag_cmp(m2, b2) > 0)
        {
            flint_printf("FAIL\n\n");
            acb_printn(s, 30, 0); flint_printf("\n\n");
            mag_printd(m, 10); flint_printf("   "); mag_printd(b, 10); flint_printf("\n\n");
            mag_printd(m1, 10); flint_printf("   "); mag_printd(b1, 10); flint_printf("\n\n");
            mag_printd(m2, 10); flint_printf("   "); mag_printd(b2, 10); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(s);
        _acb_vec_clear(v, 3);

        mag_clear(b);
        mag_clear(b1);
        mag_clear(b2);
        mag_clear(m);
        mag_clear(m1);
        mag_clear(m2);
    }

    TEST_FUNCTION_END(state);
}
