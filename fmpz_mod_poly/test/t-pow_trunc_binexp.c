/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Lina Kulakova

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
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);
    

    flint_printf("pow_trunc_binexp....");
    fflush(stdout);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c;
        fmpz_t p;
        slong e, trunc;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 30));
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        fmpz_mod_poly_set(c, a);

        fmpz_mod_poly_pow_trunc_binexp(b, a, e, trunc);
        fmpz_mod_poly_pow_trunc_binexp(c, c, e, trunc);

        result = (fmpz_mod_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, p = %wu, exp = %wd, trunc = %wd\n",
                a->length, a->p, e, trunc);
            flint_printf("a:\n"); fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
    }

    /* Check powering against naive method */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c;
        fmpz_t p;
        slong e, trunc;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 30));
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        fmpz_mod_poly_pow_trunc_binexp(b, a, e, trunc);
        fmpz_mod_poly_pow(c, a, e);
        fmpz_mod_poly_truncate(c, trunc);
        
        result = (fmpz_mod_poly_equal(b, c)
            || (a->length == 0 && e == 0 && c->length == 1 && c->coeffs[0] == 1));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, p = %wu, exp = %wd, trunc = %wd\n",
                a->length, a->p, e, trunc);
            flint_printf("a:\n"); fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b:\n"); fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("c:\n"); fmpz_mod_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
