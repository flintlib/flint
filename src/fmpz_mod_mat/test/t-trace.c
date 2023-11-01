/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_trace, state)
{
    slong i;

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mat_t A, B, AB, BA;
        fmpz_t trab, trba, mod;
        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_init(mod);
        fmpz_randtest_not_zero(mod, state, 100);
        fmpz_abs(mod, mod);

        fmpz_mod_mat_init(A, m, n, mod);
        fmpz_mod_mat_init(B, n, m, mod);
        fmpz_mod_mat_init(AB, m, m, mod);
        fmpz_mod_mat_init(BA, n, n, mod);

        fmpz_init(trab);
        fmpz_init(trba);

        fmpz_mod_mat_randtest(A, state);
        fmpz_mod_mat_randtest(B, state);

        fmpz_mod_mat_mul(AB, A, B);
        fmpz_mod_mat_mul(BA, B, A);

        fmpz_mod_mat_trace(trab, AB);
        fmpz_mod_mat_trace(trba, BA);

        if (!fmpz_equal(trab, trba))
        {
            flint_printf("FAIL:\n");
            fmpz_mod_mat_print_pretty(A), flint_printf("\n");
            fmpz_mod_mat_print_pretty(B), flint_printf("\n");
            fmpz_mod_mat_print_pretty(AB), flint_printf("\n");
            fmpz_mod_mat_print_pretty(BA), flint_printf("\n");
            flint_printf("tr(AB): "),  fmpz_print(trab),    flint_printf("\n");
            flint_printf("tr(BA): "),  fmpz_print(trba),    flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(AB);
        fmpz_mod_mat_clear(BA);
        fmpz_clear(trab);
        fmpz_clear(trba);
        fmpz_clear(mod);
    }

    TEST_FUNCTION_END(state);
}
