/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "qqbar.h"

int
qqbar_equal_fmpq_poly_val2(const qqbar_t x, const fmpq_poly_t f, const qqbar_t y)
{
    int found;
    qqbar_t v;
    qqbar_init(v);
    qqbar_evaluate_fmpq_poly(v, f, y);
    found = qqbar_equal(v, x);
    qqbar_clear(v);
    return found;
}

TEST_FUNCTION_START(qqbar_equal_fmpq_poly_val, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        qqbar_ptr v;
        fmpq_poly_t f;
        int equal1, equal2;
        slong d;

        qqbar_init(x);
        qqbar_init(y);
        fmpq_poly_init(f);

        if (n_randint(state, 2))
            qqbar_randtest(y, state, 4, 10);
        else
            qqbar_randtest(y, state, 2, 100);

        do {
            fmpq_poly_randtest(f, state, 10, 10);
        } while (f->length == 0);

        switch (n_randint(state, 3))
        {
            case 0:
                qqbar_randtest(x, state, 4, 20);
                break;
            case 1:
                qqbar_evaluate_fmpq_poly(x, f, y);
                break;
            default:
                qqbar_evaluate_fmpq_poly(x, f, y);
                d = qqbar_degree(x);
                v = _qqbar_vec_init(d);
                qqbar_conjugates(v, x);
                qqbar_set(x, v + n_randint(state, d));
                _qqbar_vec_clear(v, d);
        }

        equal1 = qqbar_equal_fmpq_poly_val(x, f, y);
        equal2 = qqbar_equal_fmpq_poly_val2(x, f, y);

        if (equal1 != equal2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("f = "); fmpq_poly_print(f); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        fmpq_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
