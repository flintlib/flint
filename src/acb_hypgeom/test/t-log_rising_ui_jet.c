/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_log_rising_ui_jet, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z;
        acb_ptr s1, s2, s3, s4;
        slong prec, ebits, len, len1, len2, len3, len4;
        ulong n, m;

        prec = 2 + n_randint(state, 100);
        len = n_randint(state, 5);
        len1 = len + n_randint(state, 2);
        len2 = len + n_randint(state, 2);
        len3 = len + n_randint(state, 2);
        len4 = len + n_randint(state, 2);

        if (n_randint(state, 10) == 0)
            ebits = 100;
        else
            ebits = 10;
        ebits = 2;

        acb_init(z);
        s1 = _acb_vec_init(FLINT_MAX(1, len1));
        s2 = _acb_vec_init(FLINT_MAX(1, len2));
        s3 = _acb_vec_init(FLINT_MAX(1, len3));
        s4 = _acb_vec_init(FLINT_MAX(1, len4));

        acb_randtest(z, state, prec, ebits);
        acb_randtest(s1, state, prec, 10);
        acb_randtest(s2, state, prec, 10);
        n = n_randint(state, 8);
        m = n_randint(state, 8);

        acb_hypgeom_log_rising_ui_jet(s1, z, n + m, len1, prec);
        acb_hypgeom_log_rising_ui_jet(s2, z, n, len2, prec);
        acb_add_ui(s4, z, n, prec);
        acb_hypgeom_log_rising_ui_jet(s3, s4, m, len3, prec);
        _acb_vec_add(s4, s2, s3, len, prec);

        if (!_acb_poly_overlaps(s1, len, s4, len))
        {
            flint_printf("FAIL\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("n = %wu, m = %lu\n\n", n, m);
            flint_printf("z = "); acb_printn(z, 1000, 0); flint_printf("\n\n");
            flint_printf("s1 = "); _acb_vec_printn(s1, len, 30, 0); flint_printf("\n\n");
            flint_printf("s2 = "); _acb_vec_printn(s2, len, 30, 0); flint_printf("\n\n");
            flint_printf("s3 = "); _acb_vec_printn(s3, len, 30, 0); flint_printf("\n\n");
            flint_printf("s4 = "); _acb_vec_printn(s4, len, 30, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z);
        _acb_vec_clear(s1, FLINT_MAX(1, len1));
        _acb_vec_clear(s2, FLINT_MAX(1, len2));
        _acb_vec_clear(s3, FLINT_MAX(1, len3));
        _acb_vec_clear(s4, FLINT_MAX(1, len4));
    }

    TEST_FUNCTION_END(state);
}
