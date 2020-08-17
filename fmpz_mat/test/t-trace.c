/*
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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("trace....");
    fflush(stdout);

    

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, AB, BA;
        fmpz_t trab, trba;
        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, m);
        fmpz_mat_init(AB, m, m);
        fmpz_mat_init(BA, n, n);

        fmpz_init(trab);
        fmpz_init(trba);

        fmpz_mat_randtest(A, state, 1 + n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1 + n_randint(state, 100));

        fmpz_mat_mul(AB, A, B);
        fmpz_mat_mul(BA, B, A);

        fmpz_mat_trace(trab, AB);
        fmpz_mat_trace(trba, BA);

        if (!fmpz_equal(trab, trba))
        {
            flint_printf("FAIL:\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            fmpz_mat_print_pretty(B), flint_printf("\n");
            fmpz_mat_print_pretty(AB), flint_printf("\n");
            fmpz_mat_print_pretty(BA), flint_printf("\n");
            flint_printf("tr(AB): "),  fmpz_print(trab),    flint_printf("\n");
            flint_printf("tr(BA): "),  fmpz_print(trba),    flint_printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(AB);
        fmpz_mat_clear(BA);
        fmpz_clear(trab);
        fmpz_clear(trba);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
