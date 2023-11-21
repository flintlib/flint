/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_gamma_taylor_tab, state)
{

    {
        slong n, prec, maxprec;
        arb_t c;
        arf_t d;
        arb_ptr v, f;

        v = _arb_vec_init(ARB_HYPGEOM_GAMMA_TAB_NUM);
        f = _arb_vec_init(2);

        arb_one(f);
        arb_one(f + 1);
        arf_init(d);

        _arb_poly_rgamma_series(v, f, 2, ARB_HYPGEOM_GAMMA_TAB_NUM,
            ARB_HYPGEOM_GAMMA_TAB_PREC + 2 * ARB_HYPGEOM_GAMMA_TAB_NUM);

        for (n = 1; n < ARB_HYPGEOM_GAMMA_TAB_NUM; n++)
        {
            maxprec = arb_hypgeom_gamma_coeffs[n].nlimbs;
            maxprec *= 64;

            for (prec = 2; prec <= ARB_HYPGEOM_GAMMA_TAB_PREC; prec += (prec < maxprec - 64 ? n_randint(state, 64) : 1))
            {
                if (_arb_hypgeom_gamma_coeff_shallow(arb_midref(c), arb_radref(c), n, prec))
                {
                    if (!arb_contains(c, v + n))
                    {
                        flint_printf("FAIL\n\n");
                        flint_printf("prec = %wd, n = %wd\n\n", prec, n);
                        flint_printf("c = "); arb_printn(c, 1000, 0); flint_printf("\n\n");
                        flint_printf("v = "); arb_printn(v + n, 1000, 0); flint_printf("\n\n");
                        flint_abort();
                    }
                }
            }

            _arb_hypgeom_gamma_coeff_shallow(arb_midref(c), arb_radref(c), n, maxprec);
            arf_set_round(d, arb_midref(v + n), maxprec, ARF_RND_DOWN);

            if (!arf_equal(arb_midref(c), d))
            {
                flint_printf("FAIL\n\n");
                flint_printf("prec = %wd, n = %wd\n\n", prec, n);
                flint_printf("c = "); arb_printn(c, 1000, 0); flint_printf("\n\n");
                flint_printf("d = "); arf_printd(d, 1000); flint_printf("\n\n");
                flint_abort();
            }
        }

        _arb_vec_clear(v, ARB_HYPGEOM_GAMMA_TAB_NUM);
        _arb_vec_clear(f, 2);
        arf_clear(d);
    }

    TEST_FUNCTION_END(state);
}
