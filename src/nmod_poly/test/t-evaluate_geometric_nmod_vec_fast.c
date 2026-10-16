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
#include "nmod.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(nmod_poly_evaluate_geometric_nmod_vec_fast, state)
{
    int i, result = 1;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P;
        nn_ptr y, z;
        ulong mod, r;
        ulong n, npoints;

        /* will do a few tests with repeated points */
        int repeated_points = (i % 10 == 5);

        /* number of points */
        npoints = (i < 20) ? i : n_randint(state, 500);

        /* poly length */
        n = n_randint(state, 1000);

        /* modulus and geometric progression */
        if (repeated_points) 
        {
            mod = 2 + n_randint(state, 5);
            do
            {
                r = n_randint(state, mod);
            }
            while (n_gcd(r, mod) != 1);
        }
        else
        {
            do 
            { 
                mod = n_randtest_prime(state, 1); 
            }
            while (mod <= 2*FLINT_MAX(npoints, n) + 1);
            r = n_primitive_root_prime(mod);
        }

        nmod_poly_init(P, mod);
        y = _nmod_vec_init(npoints);
        z = _nmod_vec_init(npoints);

        nmod_poly_randtest(P, state, n);

        nmod_poly_evaluate_geometric_nmod_vec_iter(y, P, r, npoints);
        nmod_poly_evaluate_geometric_nmod_vec_fast(z, P, r, npoints);

        result = _nmod_vec_equal(y, z, npoints);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("mod=%wu, n=%wd, npoints=%wd\n\n", mod, n, npoints);
            flint_printf("P: "); nmod_poly_print_pretty(P, "x"); flint_printf("\n\n");
            flint_printf("y: "); _nmod_vec_print_pretty(y, npoints, P->mod); flint_printf("\n\n");
            flint_printf("z: "); _nmod_vec_print_pretty(z, npoints, P->mod); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(P);
        _nmod_vec_clear(y);
        _nmod_vec_clear(z);
    }

    TEST_FUNCTION_END(state);
}
