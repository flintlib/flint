/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "calcium.h"
#include "fexpr.h"

TEST_FUNCTION_START(fexpr_call_vec, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fexpr_t f, f2, f3, g, v, a;
        fexpr_ptr args;
        fmpz_t c;
        slong i, len;

        len = n_randint(state, 20);
        fmpz_init(c);
        fexpr_init(f);
        fexpr_init(f2);
        fexpr_init(g);
        fexpr_init(a);
        args = _fexpr_vec_init(len);

        if (n_randint(state, 2))
            fexpr_set_symbol_str(f, "Foo");
        else
            fexpr_set_symbol_str(f, "FooVeryLongSymbolName");

        for (i = 0; i < len; i++)
        {
            fmpz_randtest(c, state, 200);
            fexpr_set_fmpz(args + i, c);
        }

        fexpr_call_vec(g, f, args, len);

        for (i = 0; i < len; i++)
        {
            fexpr_func(f2, g);
            fexpr_view_func(f3, g);

            fexpr_arg(a, g, i);
            fexpr_view_arg(v, g, i);

            if (fexpr_nargs(g) != len || !fexpr_equal(f, f2) || !fexpr_equal(f, f3) || !fexpr_equal(a, args + i) || !fexpr_equal(v, args + i))
            {
                flint_printf("FAIL\n\n");
                flint_printf("len = %wd\n", len);
                flint_printf("g = "); fexpr_print(g); printf("\n");
                flint_printf("arg = "); fexpr_print(args + i); printf("\n");
                flint_printf("a = "); fexpr_print(a); printf("\n");
                flint_printf("v = "); fexpr_print(v); printf("\n");
                flint_printf("f2 = "); fexpr_print(f2); printf("\n");
                flint_abort();
            }
        }

        fexpr_clear(f);
        fexpr_clear(f2);
        fexpr_clear(g);
        fexpr_clear(a);
        _fexpr_vec_clear(args, len);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
