/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_express_in_field, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_poly_t f, g;
        qqbar_t x, alpha;
        slong prec, max_bits;

        fmpq_poly_init(f);
        fmpq_poly_init(g);

        qqbar_init(x);
        qqbar_init(alpha);

        if (n_randint(state, 2))
            qqbar_randtest(alpha, state, 4, 100);
        else
            qqbar_randtest(alpha, state, 6, 10);

        if (n_randint(state, 4) == 0)
            fmpq_poly_randtest(f, state, qqbar_degree(alpha), 100);
        else
            fmpq_poly_randtest(f, state, qqbar_degree(alpha), 8);

        max_bits = _fmpz_vec_max_bits(f->coeffs, f->length);
        max_bits = FLINT_ABS(max_bits);
        max_bits = FLINT_MAX(max_bits, fmpz_bits(f->den));

        qqbar_evaluate_fmpq_poly(x, f, alpha);

        for (prec = 64; ; prec *= 2)
        {
            if (qqbar_express_in_field(g, alpha, x, max_bits, 0, prec))
                break;

            if (prec > 10000)
            {
                flint_printf("FAIL!\n");
                flint_printf("alpha = "); qqbar_print(alpha); flint_printf("\n\n");
                flint_printf("f = "); fmpq_poly_print(f); flint_printf("\n\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("g = "); fmpq_poly_print(g); flint_printf("\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_abort();
            }
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);

        qqbar_clear(alpha);
        qqbar_clear(x);
    }

    TEST_FUNCTION_END(state);
}
