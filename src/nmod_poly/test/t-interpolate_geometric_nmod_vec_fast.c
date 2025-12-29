/*
    Copyright (C) 2025, Vincent Neiger, Éric Schost
    Copyright (C) 2025, Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_interpolate_geometric_nmod_vec_fast, state)
{
    int i, result = 1;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P, Q;
        nn_ptr y;
        ulong mod, r;
        slong n, npoints;

        npoints = 1 + (i == 0 ? 1 : n_randint(state, 100));
        n = 1 + n_randint(state, npoints);
        do 
        { 
            mod = n_randtest_prime(state, 1); 
        }
        while (mod <= 2*FLINT_MAX(npoints, n) + 1); // minimum limit for maximum order r

        nmod_poly_init(P, mod);
        nmod_poly_init(Q, mod);

        r = n_primitive_root_prime(mod);
        y = _nmod_vec_init(npoints);

        nmod_poly_randtest(P, state, n);

        nmod_poly_evaluate_geometric_nmod_vec_fast(y, P, r, npoints);
        nmod_poly_interpolate_geometric_nmod_vec_fast(Q, r, y, npoints);

        result = nmod_poly_equal(P, Q);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("mod=%wu, n=%wd, npoints=%wd\n\n", mod, n, npoints);
            nmod_poly_print_pretty(P, "x"), flint_printf("\n\n");
            nmod_poly_print_pretty(Q, "x"), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(P);
        nmod_poly_clear(Q);
        _nmod_vec_clear(y);
    }

    TEST_FUNCTION_END(state);
}
