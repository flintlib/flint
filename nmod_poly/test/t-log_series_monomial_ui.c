/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

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
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result = 1;
    FLINT_TEST_INIT(state);
    

    flint_printf("log_series_monomial_ui....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, logA, res;
        slong n;
        mp_limb_t mod;
        ulong power;
        mp_limb_t coeff;

        mod = n_randtest_prime(state, 0);
        n = n_randtest(state) % 100;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(logA, mod);
        nmod_poly_init(res, mod);

        coeff = n_randlimb(state) % mod;
        power = 1 + n_randint(state, 2*n + 1);

        nmod_poly_set_coeff_ui(A, 0, UWORD(1));
        nmod_poly_set_coeff_ui(A, power, coeff);

        nmod_poly_log_series(logA, A, n);
        nmod_poly_log_series_monomial_ui(res, coeff, power, n);

        result = nmod_poly_equal(logA, res);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd, mod = %wu\n", n, mod);
            flint_printf("power = %wu, coeff = %wu\n", power, coeff);
            flint_printf("A: "); nmod_poly_print(A), flint_printf("\n\n");
            flint_printf("log(A): "); nmod_poly_print(logA), flint_printf("\n\n");
            flint_printf("res: "); nmod_poly_print(res), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(A);
        nmod_poly_clear(logA);
        nmod_poly_clear(res);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
