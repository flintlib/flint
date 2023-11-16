/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arith.h"
#include "arb.h"

TEST_FUNCTION_START(arb_euler_number_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t b1, b2;
        fmpz_t n, bv;
        slong prec1, prec2, acc1, acc2;

        fmpz_init(n);
        fmpz_init(bv);
        arb_init(b1);
        arb_init(b2);

        if (iter % 100 == 0)
        {
            fmpz_randtest(n, state, 1 + n_randint(state, 1000));
            fmpz_abs(n, n);
        }
        else
        {
            fmpz_randtest(n, state, 1 + n_randint(state, 20));
            fmpz_abs(n, n);
        }

        prec1 = 2 + n_randint(state, 1000);

        if (iter % 100 == 0)
        {
            fmpz_one(n);
            fmpz_mul_2exp(n, n, 64);
            arb_euler_number_fmpz(b1, n, prec1);
            arb_set_str(b2, "3.69772629631570244101368532955e+343764398125858671130", 64);
            prec2 = 64;
        }
        else if (iter % 100 == 1)
        {
            fmpz_one(n);
            fmpz_mul_2exp(n, n, 64);
            fmpz_add_ui(n, n, 2);
            arb_euler_number_fmpz(b1, n, prec1);
            arb_set_str(b2, "-5.09958051083474382772175848886e+343764398125858671168", 64);
            prec2 = 64;
        }
        else
        {
            prec2 = 2 + n_randint(state, 1000);
            arb_euler_number_fmpz(b1, n, prec1);
            arb_euler_number_fmpz(b2, n, prec2);
        }

        if (fmpz_cmp_ui(n, 300) < 0)
        {
            arith_euler_number(bv, fmpz_get_ui(n));
        }

        if (!arb_overlaps(b1, b2) || (fmpz_cmp_ui(n, 100) < 0
            && !arb_contains_fmpz(b1, bv)))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("b1 = "); arb_printn(b1, 50, 0); flint_printf("\n\n");
            flint_printf("b2 = "); arb_printn(b2, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        acc1 = arb_rel_accuracy_bits(b1);
        acc2 = arb_rel_accuracy_bits(b2);

        if (acc1 < prec1 - 4 || acc2 < prec2 - 4)
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd\n", prec1, acc1);
            flint_printf("prec2 = %wd, acc2 = %wd\n", prec2, acc2);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("b1 = "); arb_printn(b1, 50, 0); flint_printf("\n\n");
            flint_printf("b2 = "); arb_printn(b2, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(b1);
        arb_clear(b2);
        fmpz_clear(n);
        fmpz_clear(bv);
    }

    TEST_FUNCTION_END(state);
}
