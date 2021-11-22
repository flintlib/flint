/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int result;
    char *str;
    fmpz_poly_t a;

    FLINT_TEST_INIT(state);
    
    flint_printf("get_str_pretty....");
    fflush(stdout);

    fmpz_poly_init(a);

    str = fmpz_poly_get_str_pretty(a, "t");
    result = strcmp(str, "0") == 0;
    if (!result)
    {
        flint_printf("FAIL:\n");
        flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n");
        flint_printf("str(a) = {%s}\n", str);
        fflush(stdout);
        flint_abort();
    }
    flint_free(str);

    fmpz_poly_set_si(a, -2);
    str = fmpz_poly_get_str_pretty(a, "t");
    result = strcmp(str, "-2") == 0;
    if (!result)
    {
        flint_printf("FAIL:\n");
        flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n");
        flint_printf("str(a) = {%s}\n", str);
        fflush(stdout);
        flint_abort();
    }
    flint_free(str);

    fmpz_poly_set_coeff_si(a, 3, 1);
    str = fmpz_poly_get_str_pretty(a, "t");
    result = strcmp(str, "t^3-2") == 0;
    if (!result)
    {
        flint_printf("FAIL:\n");
        flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n");
        flint_printf("str(a) = {%s}\n", str);
        fflush(stdout);
        flint_abort();
    }
    flint_free(str);
    fmpz_poly_clear(a);
    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
