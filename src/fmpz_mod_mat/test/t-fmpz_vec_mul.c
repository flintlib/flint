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
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_fmpz_vec_mul, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mat_t A, B, C;
        fmpz * a, * c;
        fmpz ** aa, ** cc;
        slong j, m, n, alen;
        fmpz_t mod;

        fmpz_init(mod);
        fmpz_randtest_not_zero(mod, state, 200);
        fmpz_abs(mod, mod);

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        alen = n_randint(state, 50);

        fmpz_mod_mat_init(C, 1, n, mod);
        fmpz_mod_mat_init(A, 1, m, mod);
        fmpz_mod_mat_init(B, m, n, mod);
        c = _fmpz_vec_init(n);
        a = _fmpz_vec_init(alen);

        fmpz_mod_mat_randtest(B, state);
        _fmpz_vec_randtest(c, state, n, n_randint(state, 200) + 1);
        _fmpz_vec_randtest(a, state, alen, n_randint(state, 200) + 1);

        cc = FLINT_ARRAY_ALLOC(n, fmpz*);
        for (j = 0; j < n; j++)
        {
            cc[j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init_set(cc[j], c + j);
        }

        aa = FLINT_ARRAY_ALLOC(alen, fmpz*);
        for (j = 0; j < alen; j++)
        {
            aa[j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init_set(aa[j], a + j);
        }

        fmpz_mod_mat_fmpz_vec_mul(c, a, alen, B);
        fmpz_mod_mat_fmpz_vec_mul_ptr(cc, (const fmpz * const *)aa, alen, B);

        /* supposed to match mul of the chopped or zero-extended a */
        for (j = 0; j < m && j < alen; j++)
            fmpz_mod(fmpz_mod_mat_entry(A, 0, j), a + j, mod);

        fmpz_mod_mat_mul(C, A, B);

        for (j = 0; j < n; j++)
        {
            if (!fmpz_equal(fmpz_mod_mat_entry(C, 0, j), c + j) ||
                !fmpz_equal(fmpz_mod_mat_entry(C, 0, j), cc[j]))
            {
                flint_printf("FAIL: wrong answer\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(mod);
        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(C);
        _fmpz_vec_clear(c, n);
        _fmpz_vec_clear(a, alen);

        for (j = 0; j < n; j++)
        {
            fmpz_clear(cc[j]);
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
