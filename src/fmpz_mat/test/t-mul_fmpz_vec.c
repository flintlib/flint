/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_mul_fmpz_vec, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        fmpz * b, * c;
        fmpz ** bb, ** cc;
        slong j, m, n, blen;

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        blen = n_randint(state, 50);

        fmpz_mat_init(C, m, 1);
        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, 1);
        c = _fmpz_vec_init(m);
        b = _fmpz_vec_init(blen);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        _fmpz_vec_randtest(c, state, m, n_randint(state, 200) + 1);
        _fmpz_vec_randtest(b, state, blen, n_randint(state, 200) + 1);

        cc = FLINT_ARRAY_ALLOC(m, fmpz*);
        for (j = 0; j < m; j++)
        {
            cc[j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init_set(cc[j], c + j);
        }

        bb = FLINT_ARRAY_ALLOC(blen, fmpz*);
        for (j = 0; j < blen; j++)
        {
            bb[j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init_set(bb[j], b + j);
        }

        fmpz_mat_mul_fmpz_vec(c, A, b, blen);
        fmpz_mat_mul_fmpz_vec_ptr(cc, A, (const fmpz * const *)bb, blen);

        /* supposed to match mul of the chopped or zero-extended b */
        for (j = 0; j < n && j < blen; j++)
            fmpz_set(fmpz_mat_entry(B, j, 0), b + j);

        fmpz_mat_mul(C, A, B);

        for (j = 0; j < m; j++)
        {
            if (!fmpz_equal(fmpz_mat_entry(C, j, 0), c + j) ||
                !fmpz_equal(fmpz_mat_entry(C, j, 0), cc[j]))
            {
                flint_printf("FAIL: wrong answer\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        _fmpz_vec_clear(c, m);
        _fmpz_vec_clear(b, blen);

        for (j = 0; j < m; j++)
        {
            fmpz_clear(cc[j]);
            flint_free(cc[j]);
        }
        flint_free(cc);

        for (j = 0; j < blen; j++)
        {
            fmpz_clear(bb[j]);
            flint_free(bb[j]);
        }
        flint_free(bb);
    }

    TEST_FUNCTION_END(state);
}
