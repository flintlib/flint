/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_charpoly_berkowitz, state)
{
    slong m, n, rep, i;
    ulong mod;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C, D;
        nmod_poly_t f, g;

        m = n_randint(state, 10);
        n = m;

        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(C, m, m, mod);
        nmod_mat_init(D, n, n, mod);
        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_mul(C, A, B);
        nmod_mat_mul(D, B, A);

        nmod_mat_charpoly_berkowitz(f, C);
        nmod_mat_charpoly_berkowitz(g, D);

        if (!nmod_poly_equal(f, g))
            TEST_FUNCTION_FAIL(
                    "charpoly(AB) != charpoly(BA)\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "cp(AB) = %{nmod_poly}\n"
                    "cp(BA) = %{nmod_poly}\n",
                    A, B, f, g);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A;
        nmod_poly_t f, g;

        m = n_randint(state, 10);
        n = m;

        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);

        nmod_mat_randtest(A, state);

        nmod_mat_charpoly_berkowitz(f, A);

        for (i = 0; i < 10; i++)
           nmod_mat_similarity(A, n_randint(state, m), n_randint(state, mod));

        nmod_mat_charpoly(g, A);

        if (!nmod_poly_equal(f, g))
            TEST_FUNCTION_FAIL(
                    "charpoly(P^{-1}AP) != charpoly(A)\n"
                    "A = %{nmod_mat}\n"
                    "cp(A) = %{nmod_poly}\n"
                    "cp(P^{-1}AP) = %{nmod_poly}\n",
                    A, f, g);

        nmod_mat_clear(A);
        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
