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
#include "fexpr_builtin.h"

#ifdef __GNUC__
# define strcmp __builtin_strcmp
#else
# include <string.h>
#endif

TEST_FUNCTION_START(fexpr_builtins, state)
{
    {
        slong i, j;

        for (i = 0; i < FEXPR_BUILTIN_LENGTH; i++)
        {
            if (fexpr_builtin_table[i].symbol != i)
            {
                flint_printf("FAIL (missing)\n");
                flint_printf("%s\n", fexpr_builtin_name(i));
                flint_abort();
            }

            if (i < FEXPR_BUILTIN_LENGTH - 1 &&
                strcmp(fexpr_builtin_name(i), fexpr_builtin_name(i + 1)) >= 0)
            {
                flint_printf("FAIL (order)\n");
                flint_printf("%s\n", fexpr_builtin_name(i));
                flint_printf("%s\n", fexpr_builtin_name(i + 1));
                flint_abort();
            }

            j = fexpr_builtin_lookup(fexpr_builtin_name(i));

            if (i != j)
            {
                flint_printf("FAIL (lookup)\n");
                flint_printf("%s\n", fexpr_builtin_name(i));
                flint_abort();
            }
        }
    }

    if (fexpr_builtin_lookup("") != -1 || fexpr_builtin_lookup("FooBarBaz") != -1 ||
        fexpr_builtin_lookup("ZZZZZZZ") != -1)
    {
        flint_printf("FAIL (lookup)\n");
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
