/*
    Copyright (C) 2020, 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("atan_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B;
        slong n1;

        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);

        ca_poly_randtest_rational(A, state, 10, 10, ctx);

        n1 = n_randint(state, 10);

        if (n_randint(state, 2))
        {
            ca_poly_atan_series(B, A, n1, ctx);
        }
        else
        {
            ca_poly_set(B, A, ctx);
            ca_poly_atan_series(B, B, n1, ctx);
        }

        if (!(A->length == 0 || ca_check_is_zero(A->coeffs, ctx) != T_FALSE))
        {
            /* ... todo: test something here ... */
        }

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);

        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
