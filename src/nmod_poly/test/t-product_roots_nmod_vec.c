/*
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

TEST_FUNCTION_START(nmod_poly_product_roots_nmod_vec, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P, Q, tmp;
        mp_ptr x;
        mp_limb_t mod;
        slong j, n;

        n = n_randint(state, 100);
        mod = n_randtest_prime(state, 0);

        nmod_poly_init(P, mod);
        nmod_poly_init(Q, mod);
        nmod_poly_init(tmp, mod);
        x = _nmod_vec_init(n);
        _nmod_vec_randtest(x, state, n, P->mod);

        nmod_poly_product_roots_nmod_vec(P, x, n);

        nmod_poly_set_coeff_ui(Q, 0, UWORD(1));

        for (j = 0; j < n; j++)
        {
            nmod_poly_zero(tmp);
            nmod_poly_set_coeff_ui(tmp, 1, UWORD(1));
            nmod_poly_set_coeff_ui(tmp, 0, n_negmod(x[j], mod));
            nmod_poly_mul(Q, Q, tmp);
        }

        result = (nmod_poly_equal(P, Q));
        if (!result)
        {
            flint_printf("FAIL (P != Q):\n");
            nmod_poly_print(P), flint_printf("\n\n");
            nmod_poly_print(Q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(P);
        nmod_poly_clear(Q);
        nmod_poly_clear(tmp);
        _nmod_vec_clear(x);
    }

    TEST_FUNCTION_END(state);
}
