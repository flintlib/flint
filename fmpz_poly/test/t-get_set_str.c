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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("get_set_str....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        char * str;
        int ans;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 100), n_randint(state, 200));

        str = fmpz_poly_get_str(a);
        ans = fmpz_poly_set_str(b, str);

        result = (ans == 0 && fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        flint_free(str);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
