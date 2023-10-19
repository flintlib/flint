/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_fmpz_vec_mul, state)
{
    slong i;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C;
        fmpz * a;
        fmpq * c;
        fmpz ** aa;
        fmpq ** cc;
        slong j, m, n, alen;

        m = n_randint(state, 40);
        n = n_randint(state, 40);
        alen = n_randint(state, 40);

        fmpq_mat_init(C, 1, n);
        fmpq_mat_init(A, 1, m);
        fmpq_mat_init(B, m, n);
        c = _fmpq_vec_init(n);
        a = _fmpz_vec_init(alen);

        fmpq_mat_randtest(B, state, n_randint(state, 200) + 1);
        _fmpq_vec_randtest(c, state, n, n_randint(state, 200) + 1);
        _fmpz_vec_randtest(a, state, alen, n_randint(state, 200) + 1);

        cc = FLINT_ARRAY_ALLOC(n, fmpq*);
        for (j = 0; j < n; j++)
        {
            cc[j] = FLINT_ARRAY_ALLOC(1, fmpq);
            fmpq_init(cc[j]);
            fmpq_set(cc[j], c + j);
        }

        aa = FLINT_ARRAY_ALLOC(alen, fmpz*);
        for (j = 0; j < alen; j++)
        {
            aa[j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init_set(aa[j], a + j);
        }

        fmpq_mat_fmpz_vec_mul(c, a, alen, B);
        fmpq_mat_fmpz_vec_mul_ptr(cc, (const fmpz * const *)aa, alen, B);

        /* supposed to match mul of the chopped or zero-extended a */
        for (j = 0; j < m && j < alen; j++)
            fmpq_set_fmpz(fmpq_mat_entry(A, 0, j), a + j);

        fmpq_mat_mul(C, A, B);

        for (j = 0; j < n; j++)
        {
            if (!fmpq_equal(fmpq_mat_entry(C, 0, j), c + j) ||
                !fmpq_equal(fmpq_mat_entry(C, 0, j), cc[j]))
            {
                flint_printf("FAIL: wrong answer\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        _fmpq_vec_clear(c, n);
        _fmpz_vec_clear(a, alen);

        for (j = 0; j < n; j++)
        {
            fmpq_clear(cc[j]);
            flint_free(cc[j]);
        }
        flint_free(cc);

        for (j = 0; j < alen; j++)
        {
            fmpz_clear(aa[j]);
            flint_free(aa[j]);
        }
        flint_free(aa);
    }

    TEST_FUNCTION_END(state);
}
