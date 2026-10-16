/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_scalar_mul_si, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t a, b, c;
        acb_t s;
        slong v;
        slong prec;

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_init(s);

        prec = 2 + n_randint(state, 200);
        acb_poly_randtest(a, state, 1 + n_randint(state, 10), prec, 10);
        v = (slong) n_randtest(state);

        acb_set_si(s, v);
        acb_poly_scalar_mul(b, a, s, prec);
        acb_poly_scalar_mul_si(c, a, v, prec);

        if (!acb_poly_overlaps(b, c))
        {
            flint_printf("FAIL\n\n");
            flint_abort();
        }

        acb_poly_scalar_mul_si(a, a, v, prec);

        if (!acb_poly_overlaps(a, c))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_clear(s);
    }

    TEST_FUNCTION_END(state);
}
