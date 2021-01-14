/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "calcium.h"
#include "fexpr.h"

int main()
{
    flint_rand_t state;

    flint_printf("builtins...");
    fflush(stdout);

    flint_randinit(state);

    {
        slong i, j;

        for (i = 0; i < FEXPR_BUILTIN_LENGTH; i++)
        {
            if (fexpr_builtins[i].symbol != i)
            {
                flint_printf("FAIL (missing)\n");
                flint_printf("%s\n", fexpr_builtins[i].string);
                flint_abort();
            }

            if (i < FEXPR_BUILTIN_LENGTH - 1 &&
                strcmp(fexpr_builtins[i].string, fexpr_builtins[i + 1].string) >= 0)
            {
                flint_printf("FAIL (order)\n");
                flint_printf("%s\n", fexpr_builtins[i].string);
                flint_printf("%s\n", fexpr_builtins[i + 1].string);
                flint_abort();
            }

            j = fexpr_get_builtin_str(fexpr_builtins[i].string);

            if (i != j)
            {
                flint_printf("FAIL (lookup)\n");
                flint_printf("%s\n", fexpr_builtins[i].string);
                flint_abort();
            }
        }
    }

    if (fexpr_get_builtin_str("") != -1 || fexpr_get_builtin_str("FooBarBaz") != -1 ||
        fexpr_get_builtin_str("ZZZZZZZ") != -1)
    {
        flint_printf("FAIL (lookup)\n");
        flint_abort();
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
