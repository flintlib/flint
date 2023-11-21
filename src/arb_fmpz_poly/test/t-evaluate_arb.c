/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "arb_fmpz_poly.h"

TEST_FUNCTION_START(arb_fmpz_poly_evaluate_arb, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t f, g, h;
        arb_t x, fx, gx, hx, fxgx;
        slong prec1, prec2, prec3;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        arb_init(x);
        arb_init(fx);
        arb_init(gx);
        arb_init(hx);
        arb_init(fxgx);

        fmpz_poly_randtest(f, state, 1 + n_randint(state, 50), 1 + n_randint(state, 1000));
        fmpz_poly_randtest(g, state, 1 + n_randint(state, 50), 1 + n_randint(state, 1000));
        fmpz_poly_add(h, f, g);

        arb_randtest(x, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));
        arb_randtest(fx, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));
        arb_randtest(gx, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));
        arb_randtest(hx, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));
        arb_randtest(fxgx, state, 1 + n_randint(state, 2000), 1 + n_randint(state, 100));

        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);
        prec3 = 2 + n_randint(state, 2000);

        switch (n_randint(state, 6))
        {
            case 0:
                arb_fmpz_poly_evaluate_arb_horner(fx, f, x, prec1);
                break;
            case 1:
                arb_fmpz_poly_evaluate_arb_rectangular(fx, f, x, prec1);
                break;
            case 2:
                arb_fmpz_poly_evaluate_arb(fx, f, x, prec1);
                break;
            case 3:
                arb_set(fx, x);
                arb_fmpz_poly_evaluate_arb_horner(fx, f, fx, prec1);
                break;
            case 4:
                arb_set(fx, x);
                arb_fmpz_poly_evaluate_arb_rectangular(fx, f, fx, prec1);
                break;
            default:
                arb_set(fx, x);
                arb_fmpz_poly_evaluate_arb(fx, f, fx, prec1);
                break;
        }

        arb_fmpz_poly_evaluate_arb(gx, g, x, prec2);
        arb_fmpz_poly_evaluate_arb(hx, h, x, prec3);
        arb_add(fxgx, fx, gx, prec3);

        if (!arb_overlaps(fxgx, hx))
        {
            flint_printf("FAIL\n");
            fmpz_poly_print(f); flint_printf("\n\n");
            fmpz_poly_print(g); flint_printf("\n\n");
            fmpz_poly_print(h); flint_printf("\n\n");
            arb_printd(x, 30); flint_printf("\n\n");
            arb_printd(fx, 30); flint_printf("\n\n");
            arb_printd(gx, 30); flint_printf("\n\n");
            arb_printd(hx, 30); flint_printf("\n\n");
            arb_printd(fxgx, 30); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        arb_clear(x);
        arb_clear(fx);
        arb_clear(gx);
        arb_clear(hx);
        arb_clear(fxgx);
    }

    TEST_FUNCTION_END(state);
}
