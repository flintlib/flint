/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Generated using Claude Opus 4.8 */

#include "test_helpers.h"
#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "acb_types.h"
#include "long_extras.h"

TEST_FUNCTION_START(acb_ode_poly_taylor_shift_aps_trunc, state)
{
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        acb_poly_t f, g, ref;
        acb_t a, b;
        slong prec, bits, n, len;

        acb_poly_init(f);
        acb_poly_init(g);
        acb_poly_init(ref);
        acb_init(a);
        acb_init(b);

        prec = 2 + n_randint(state, 200);
        bits = 1 + n_randint(state, 8);

        acb_poly_randtest(f, state, n_randint(state, 10), prec, bits);

        if (n_randint(state, 3) == 0)
            acb_zero(a);
        else
            acb_randtest(a, state, prec, bits);

        n = z_randint(state, 30);
        len = n_randint(state, 12);

        acb_add_si(b, a, n, prec);
        acb_poly_taylor_shift(ref, f, b, prec);
        acb_poly_truncate(ref, len);

        if (n_randint(state, 2))
            acb_ode_poly_taylor_shift_aps_trunc(g, f, a, n, len, prec);
        else
        {
            acb_poly_set(g, f);
            acb_ode_poly_taylor_shift_aps_trunc(g, g, a, n, len, prec);
        }

        if (!acb_poly_overlaps(g, ref))
        {
            flint_printf("FAIL (overlap)\n\n");
            flint_printf("f = %{acb_poly}\n", f);
            flint_printf("a = %{acb}\n", a);
            flint_printf("n = %wd, len = %wd\n", n, len);
            flint_printf("g   = %{acb_poly}\n", g);
            flint_printf("ref = %{acb_poly}\n", ref);
            flint_abort();
        }

        acb_clear(b);
        acb_clear(a);
        acb_poly_clear(ref);
        acb_poly_clear(g);
        acb_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
