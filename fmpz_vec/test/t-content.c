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
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("content....");
    fflush(stdout);

    

    /* Check that content(a f) = abs(a) content(f) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c, d;
        fmpz * f, *g;
        slong len = n_randint(state, 100);

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);
        f = _fmpz_vec_init(len);
        g = _fmpz_vec_init(len);
        _fmpz_vec_randtest(f, state, len, 200);
        fmpz_randtest(a, state, 100);

        _fmpz_vec_content(c, f, len);
        _fmpz_vec_scalar_mul_fmpz(g, f, len, a);
        fmpz_abs(a, a);
        fmpz_mul(c, a, c);
        _fmpz_vec_content(d, g, len);

        result = (fmpz_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("d = "), fmpz_print(d), flint_printf("\n\n");
            flint_printf("vec = "), _fmpz_vec_print(f, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
        _fmpz_vec_clear(f, len);
        _fmpz_vec_clear(g, len);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
