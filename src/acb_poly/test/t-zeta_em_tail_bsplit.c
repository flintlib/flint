/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_zeta_em_tail_bsplit, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t Na, s;
        acb_ptr z1, z2, Nasx;
        slong i, M, len, prec;

        prec = 2 + n_randint(state, 400);
        len = 1 + n_randint(state, 30);
        M = n_randint(state, 40);

        acb_init(Na);
        acb_init(s);

        Nasx = _acb_vec_init(len);
        z1 = _acb_vec_init(len);
        z2 = _acb_vec_init(len);

        acb_randtest(Na, state, prec, 4);
        acb_randtest(s, state, prec, 4);

        _acb_poly_acb_invpow_cpx(Nasx, Na, s, len, prec);

        _acb_poly_zeta_em_tail_naive(z1, s, Na, Nasx, M, len, prec);
        _acb_poly_zeta_em_tail_bsplit(z2, s, Na, Nasx, M, len, prec);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(z1 + i, z2 + i))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd\n", iter);
                flint_printf("prec = %wd, len = %wd, M = %wd\n", prec, len, M);
                flint_printf("s = "); acb_printd(s, prec / 3.33); flint_printf("\n\n");
                flint_printf("Na = "); acb_printd(Na, prec / 3.33); flint_printf("\n\n");
                flint_printf("z1 = "); acb_printd(z1 + i, prec / 3.33); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printd(z2 + i, prec / 3.33); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(Na);
        acb_clear(s);
        _acb_vec_clear(Nasx, len);
        _acb_vec_clear(z1, len);
        _acb_vec_clear(z2, len);
    }

    TEST_FUNCTION_END(state);
}
