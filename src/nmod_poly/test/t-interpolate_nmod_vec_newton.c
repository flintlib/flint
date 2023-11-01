/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_interpolate_nmod_vec_newton, state)
{
    int i, result = 1;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P, Q;
        mp_ptr x, y;
        mp_limb_t mod;
        slong j, n, npoints;

        mod = n_randtest_prime(state, 0);
        npoints = n_randint(state, FLINT_MIN(100, mod));
        n = n_randint(state, npoints + 1);

        nmod_poly_init(P, mod);
        nmod_poly_init(Q, mod);
        x = _nmod_vec_init(npoints);
        y = _nmod_vec_init(npoints);

        nmod_poly_randtest(P, state, n);

        for (j = 0; j < npoints; j++)
            x[j] = j;

        nmod_poly_evaluate_nmod_vec(y, P, x, npoints);
        nmod_poly_interpolate_nmod_vec_newton(Q, x, y, npoints);

        result = nmod_poly_equal(P, Q);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("mod=%wu, n=%wd, npoints=%wd\n\n", mod, n, npoints);
            nmod_poly_print(P), flint_printf("\n\n");
            nmod_poly_print(Q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(P);
        nmod_poly_clear(Q);
        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
    }

    TEST_FUNCTION_END(state);
}
