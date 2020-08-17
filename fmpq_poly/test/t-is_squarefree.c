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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("is_squarefree....");
    fflush(stdout);

    

    /* Check that polynomials of degree <= 1 are square-free */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f;

        fmpq_poly_init(f);
        fmpq_poly_randtest(f, state, n_randint(state, 2), 100);

        result = (fmpq_poly_is_squarefree(f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n");
            abort();
        }

        fmpq_poly_clear(f);
    }

    /* Check that a^2 f is not square-free */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, f;

        fmpq_poly_init(a);
        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1, 40);
        if (a->length < 2)
        {
            fmpq_poly_clear(a);
            continue;
        }
        fmpq_poly_init(f);
        fmpq_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1, 40);

        fmpq_poly_mul(a, a, a);
        fmpq_poly_mul(f, a, f);

        result = (!fmpq_poly_is_squarefree(f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(f);
    }

    /* Check that f + N*(x^M + 1) is square-free, for N >> f, M > deg(f) */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, f;
        fmpz_t N;

        fmpq_poly_init(a);
        fmpq_poly_set_coeff_si(a, 0, WORD(1));
        fmpq_poly_set_coeff_si(a, n_randint(state, 20), WORD(1));
        if (a->length < 2)
        {
            fmpq_poly_clear(a);
            continue;
        }
        fmpq_poly_init(f);
        fmpq_poly_randtest(f, state, a->length - 2, 40);

        fmpz_init_set_ui(N, UWORD(1));
        fmpz_mul_2exp(N, N, 45 + a->length);

        fmpq_poly_scalar_mul_fmpz(a, a, N);
        fmpq_poly_add(f, a, f);

        result = fmpq_poly_is_squarefree(f);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n");
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(f);
        fmpz_clear(N);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
