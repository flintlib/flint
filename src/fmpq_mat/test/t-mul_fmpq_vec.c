/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_mul_fmpq_vec, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C;
        fmpq * b;
        fmpq * c;
        fmpq ** bb;
        fmpq ** cc;
        slong j, m, n, blen;

        m = n_randint(state, 40);
        n = n_randint(state, 40);
        blen = n_randint(state, 40);

        fmpq_mat_init(C, m, 1);
        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, 1);
        c = _fmpq_vec_init(m);
        b = _fmpq_vec_init(blen);

        fmpq_mat_randtest(A, state, n_randint(state, 200) + 1);
        _fmpq_vec_randtest(c, state, m, n_randint(state, 200) + 1);
        _fmpq_vec_randtest(b, state, blen, n_randint(state, 200) + 1);

        cc = FLINT_ARRAY_ALLOC(m, fmpq*);
        for (j = 0; j < m; j++)
        {
            cc[j] = FLINT_ARRAY_ALLOC(1, fmpq);
            fmpq_init(cc[j]);
            fmpq_set(cc[j], c + j);
        }

        bb = FLINT_ARRAY_ALLOC(blen, fmpq*);
        for (j = 0; j < blen; j++)
        {
            bb[j] = FLINT_ARRAY_ALLOC(1, fmpq);
            fmpq_init(bb[j]);
            fmpq_set(bb[j], b + j);
        }

        fmpq_mat_mul_fmpq_vec(c, A, b, blen);
        fmpq_mat_mul_fmpq_vec_ptr(cc, A, (const fmpq * const *)bb, blen);

        /* supposed to match mul of the chopped or zero-extended b */
        for (j = 0; j < n && j < blen; j++)
            fmpq_set(fmpq_mat_entry(B, j, 0), b + j);

        fmpq_mat_mul(C, A, B);

        for (j = 0; j < m; j++)
        {
            if (!fmpq_equal(fmpq_mat_entry(C, j, 0), c + j) ||
                !fmpq_equal(fmpq_mat_entry(C, j, 0), cc[j]))
            {
                flint_printf("FAIL: wrong answer\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        _fmpq_vec_clear(c, m);
        _fmpq_vec_clear(b, blen);

        for (j = 0; j < m; j++)
        {
            fmpq_clear(cc[j]);
            flint_free(cc[j]);
        }
        flint_free(cc);

        for (j = 0; j < blen; j++)
        {
            fmpq_clear(bb[j]);
            flint_free(bb[j]);
        }
        flint_free(bb);
    }

    TEST_FUNCTION_END(state);
}
