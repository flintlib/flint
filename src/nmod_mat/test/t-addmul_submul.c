/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_addmul_submul, state)
{
    slong i;

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D, T, E;
        ulong mod = n_randtest_not_zero(state);
        int type, operation;

        slong m, k, n;

        m = n_randint(state, 100);
        k = n_randint(state, 100);
        n = n_randint(state, 100);

        /* Force Strassen test */
        if (n_randint(state, 5) == 0)
        {
            m += 300;
            k += 300;
            n += 300;
        }

        type = n_randint(state, 2);
        operation = n_randint(state, 2);

        nmod_mat_init(A, m, k, mod);
        nmod_mat_init(B, k, n, mod);
        nmod_mat_init(C, m, n, mod);
        nmod_mat_init(D, m, n, mod);
        nmod_mat_init(T, m, n, mod);
        nmod_mat_init(E, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);
        nmod_mat_randtest(C, state);

        nmod_mat_mul(T, A, B);

        if (operation)
        {
            /* addmul */
            nmod_mat_add(E, C, T);

            if (type)
            {
                /* without aliasing */
                nmod_mat_addmul(D, C, A, B);
            }
            else
            {
                /* with aliasing */
                nmod_mat_set(D, C);
                nmod_mat_addmul(D, D, A, B);
            }
        }
        else
        {
            /* submul */
            nmod_mat_sub(E, C, T);
            
            if (type)
            {
                /* without aliasing */
                nmod_mat_submul(D, C, A, B);
            }
            else
            {
                /* with aliasing */
                nmod_mat_set(D, C);
                nmod_mat_submul(D, D, A, B);
            }
        }

        if (!nmod_mat_equal(D, E))
            TEST_FUNCTION_FAIL(
                    "Results not equal\n"
                    "type: %d\n"
                    "operation: %d\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "C = %{nmod_mat}\n"
                    "D = %{nmod_mat}\n"
                    "E = %{nmod_mat}\n",
                    type, operation, A, B, C, D, E);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
        nmod_mat_clear(E);
        nmod_mat_clear(T);
    }

    TEST_FUNCTION_END(state);
}
