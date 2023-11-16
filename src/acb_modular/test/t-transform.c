/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_transform, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        psl2z_t g;
        acb_t z, w1, w2, t;
        slong prec;

        psl2z_init(g);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);

        psl2z_randtest(g, state, n_randint(state, 20));
        acb_randtest(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 20));
        acb_randtest(w1, state, 1 + n_randint(state, 200), 1 + n_randint(state, 20));
        acb_randtest(w2, state, 1 + n_randint(state, 200), 1 + n_randint(state, 20));

        acb_modular_transform(w1, g, z, 2 + n_randint(state, 200));

        prec = 2 + n_randint(state, 200);

        acb_mul_fmpz(t, z, &g->a, prec);
        acb_add_fmpz(t, t, &g->b, prec);

        acb_mul_fmpz(w2, z, &g->c, prec);
        acb_add_fmpz(w2, w2, &g->d, prec);

        acb_div(w2, t, w2, prec);

        if (!acb_overlaps(w1, w2))
        {
            flint_printf("FAIL\n");
            flint_printf("g = "); psl2z_print(g); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_abort();
        }

        psl2z_clear(g);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
