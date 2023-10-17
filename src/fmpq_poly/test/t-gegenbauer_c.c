/*
    Copyright (C) 2016  Ralf Stephan

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpq_poly.h"

#define NRATS 20

TEST_FUNCTION_START(fmpq_poly_gegenbauer_c, state)
{
    fmpq_poly_t T0, T1, T2, t, tt;
    fmpq_t a, rat;
    fmpq *rats;
    slong n, d;
    flint_rand_t rand_state;

    fmpq_poly_init(T0);
    fmpq_poly_init(T1);
    fmpq_poly_init(T2);
    fmpq_poly_init(t);
    fmpq_poly_init(tt);
    fmpq_init(a);
    fmpq_init(rat);

    rats = _fmpq_vec_init(NRATS);
    fmpq_set_si(rats, 0, 1);
    for (d = 1; d < 16; d++)
        fmpq_set_si(rats + d, 1, d);
    flint_randinit(rand_state);
    for (d = 16; d < NRATS; d++)
        fmpq_randtest_not_zero(rats + d, rand_state, 32);
    flint_randclear(rand_state);

    for (d = 0; d < NRATS; d++)
    {
        fmpq_set(a, rats + d);
        fmpq_poly_gegenbauer_c(T0, 0, a);
        fmpq_poly_gegenbauer_c(T1, 1, a);

        for (n = 1; n <= 20; n++)
        {
            fmpq_poly_gegenbauer_c(T2, n+1, a);
            fmpq_poly_set(t, T1);

            /* Verify (n+1)C^a_{n+1} = 2x(n+a) C^a_n - (n+2a-1)C^a_{n-1} */
            fmpq_poly_shift_left(t, t, 1);
            fmpq_set(rat, a);
            fmpq_add_si(rat, rat, n);
            fmpq_mul_2exp(rat, rat, 1);
            fmpq_poly_scalar_mul_fmpq(t, t, rat);

            fmpq_set(rat, a);
            fmpq_mul_2exp(rat, rat, 1);
            fmpq_add_si(rat, rat, n-1);
            fmpq_poly_scalar_mul_fmpq(tt, T0, rat);
            fmpq_poly_sub(t, t, tt);
            fmpq_poly_scalar_mul_si(tt, T2, n+1);

            if (!fmpq_poly_equal(t, tt))
            {
                flint_printf("\nFAIL: n = %wd, a = ", n);
                fmpq_print(a); flint_printf("\n");
                flint_printf("t: "); fmpq_poly_print_pretty(t, "x"); flint_printf("\n");
                flint_printf("tt: "); fmpq_poly_print_pretty(tt, "x"); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_poly_swap(T0, T1);
            fmpq_poly_swap(T1, T2);
        }
    }

    fmpq_poly_clear(T0);
    fmpq_poly_clear(T1);
    fmpq_poly_clear(T2);
    fmpq_poly_clear(t);
    fmpq_poly_clear(tt);
    _fmpq_vec_clear(rats, NRATS);
    fmpq_clear(a);
    fmpq_clear(rat);

    TEST_FUNCTION_END(state);
}
