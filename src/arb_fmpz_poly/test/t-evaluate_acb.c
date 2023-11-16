/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "arb_fmpz_poly.h"

TEST_FUNCTION_START(arb_fmpz_poly_evaluate_acb, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t f, g, h;
        acb_t x, fx, gx, hx, fxgx;
        slong prec1, prec2, prec3;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        acb_init(x);
        acb_init(fx);
        acb_init(gx);
        acb_init(hx);
        acb_init(fxgx);

        fmpz_poly_randtest(f, state, 1 + n_randint(state, 50), 1 + n_randint(state, 1000));
        fmpz_poly_randtest(g, state, 1 + n_randint(state, 50), 1 + n_randint(state, 1000));
        fmpz_poly_add(h, f, g);

        acb_randtest(x, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));
        acb_randtest(fx, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));
        acb_randtest(gx, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));
        acb_randtest(hx, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));
        acb_randtest(fxgx, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));

        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);
        prec3 = 2 + n_randint(state, 2000);

        switch (n_randint(state, 6))
        {
            case 0:
                arb_fmpz_poly_evaluate_acb_horner(fx, f, x, prec1);
                break;
            case 1:
                arb_fmpz_poly_evaluate_acb_rectangular(fx, f, x, prec1);
                break;
            case 2:
                arb_fmpz_poly_evaluate_acb(fx, f, x, prec1);
                break;
            case 3:
                acb_set(fx, x);
                arb_fmpz_poly_evaluate_acb_horner(fx, f, fx, prec1);
                break;
            case 4:
                acb_set(fx, x);
                arb_fmpz_poly_evaluate_acb_rectangular(fx, f, fx, prec1);
                break;
            default:
                acb_set(fx, x);
                arb_fmpz_poly_evaluate_acb(fx, f, fx, prec1);
                break;
        }

        arb_fmpz_poly_evaluate_acb(gx, g, x, prec2);
        arb_fmpz_poly_evaluate_acb(hx, h, x, prec3);
        acb_add(fxgx, fx, gx, prec3);

        if (!acb_overlaps(fxgx, hx))
        {
            flint_printf("FAIL\n");
            fmpz_poly_print(f); flint_printf("\n\n");
            fmpz_poly_print(g); flint_printf("\n\n");
            fmpz_poly_print(h); flint_printf("\n\n");
            acb_printd(x, 30); flint_printf("\n\n");
            acb_printd(fx, 30); flint_printf("\n\n");
            acb_printd(gx, 30); flint_printf("\n\n");
            acb_printd(hx, 30); flint_printf("\n\n");
            acb_printd(fxgx, 30); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        acb_clear(x);
        acb_clear(fx);
        acb_clear(gx);
        acb_clear(hx);
        acb_clear(fxgx);
    }

    TEST_FUNCTION_END(state);
}
