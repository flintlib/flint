/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_pow_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d;
        fmpz_t e1, e2, e3;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        fmpz_init(e1);
        fmpz_init(e2);
        fmpz_init(e3);

        acb_randtest(a, state, 1 + n_randint(state, 1000), 200);
        acb_randtest(b, state, 1 + n_randint(state, 1000), 200);
        fmpz_randtest(e1, state, 200);
        fmpz_randtest(e2, state, 200);

        acb_pow_fmpz(b, a, e1, prec1);
        acb_pow_fmpz(c, a, e1, prec2);

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("e1 = "); fmpz_print(e1); flint_printf("\n\n");
            flint_abort();
        }

        /* check a^(e1+e2) = a^e1*a^e2 */
        acb_pow_fmpz(c, a, e2, prec1);
        acb_mul(d, b, c, prec1);
        fmpz_add(e3, e1, e2);
        acb_pow_fmpz(c, a, e3, prec1);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_printf("e1 = "); fmpz_print(e1); flint_printf("\n\n");
            flint_printf("e2 = "); fmpz_print(e2); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        fmpz_clear(e1);
        fmpz_clear(e2);
        fmpz_clear(e3);
    }

    TEST_FUNCTION_END(state);
}
