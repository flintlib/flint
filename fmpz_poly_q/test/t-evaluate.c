/*
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "fmpz_poly_q.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("evaluate...");
    fflush(stdout);

    

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        int ans1, ans2;
        mpq_t a, b;
        fmpz_t num, den;
        fmpz_poly_q_t f;

        mpq_init(a);
        mpq_init(b);
        fmpz_init(num);
        fmpz_init(den);
        fmpz_poly_q_init(f);
        fmpz_poly_q_randtest(f, state, n_randint(state, 10), 10, n_randint(state, 10), 10);

        fmpz_randtest(num, state, 50);
        fmpz_randtest_not_zero(den, state, 50);
        fmpz_get_mpz(mpq_numref(a), num);
        fmpz_get_mpz(mpq_denref(a), den);
        mpq_canonicalize(a);

        ans1 = fmpz_poly_q_evaluate(b, f, a);
        ans2 = fmpz_poly_q_evaluate(a, f, a);

        result = (ans1 == ans2) && mpq_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_poly_q_print(f), flint_printf("\n");
            flint_printf("num = "), fmpz_print(num), flint_printf("\n");
            flint_printf("den = "), fmpz_print(den), flint_printf("\n");
            gmp_printf("a = %Qd\n", a);
            gmp_printf("b = %Qd\n", b);
            flint_printf("ans1 = %d\n", ans1);
            flint_printf("ans2 = %d\n", ans2);
            fflush(stdout);
            flint_abort();
        }

        mpq_clear(a);
        mpq_clear(b);
        fmpz_clear(num);
        fmpz_clear(den);
        fmpz_poly_q_clear(f);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
