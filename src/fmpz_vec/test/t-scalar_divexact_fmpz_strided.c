/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_scalar_divexact_fmpz_strided, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz *Q, *A, *B, *C;
        slong r, c, i, j, Astride, Bstride, Alen, Blen, qbits, xbits;
        int alias;
        fmpz_t x;

        r = n_randint(state, 6);
        c = n_randint(state, 6);
        alias = n_randint(state, 2);
        Astride = c + n_randint(state, 3);
        Bstride = alias ? Astride : c + n_randint(state, 3);
        Alen = r * Astride;
        Blen = r * Bstride;

        /* divisors from one limb to many; quotients much smaller,
           comparable or larger */
        xbits = n_randint(state, 4) == 0 ? 1 + n_randint(state, 3000) : 1 + n_randint(state, 300);
        qbits = n_randint(state, 3) == 0 ? 1 + n_randint(state, 3000) : 1 + n_randint(state, 200);

        fmpz_init(x);
        fmpz_randtest_not_zero(x, state, xbits);

        /* Q: quotients; A = x Q; B: output with random entries in the
           gaps between rows, which must be left unchanged; C: expected
           output (per-row division) */
        Q = _fmpz_vec_init(Alen);
        A = _fmpz_vec_init(Alen);
        B = _fmpz_vec_init(Blen);
        C = _fmpz_vec_init(Blen);

        _fmpz_vec_randtest(Q, state, Alen, qbits);
        for (i = 0; i < r; i++)
            if (n_randint(state, 4) == 0)
                for (j = 0; j < c; j++)
                    fmpz_set_si(Q + i * Astride + j, (slong) n_randint(state, 5) - 2);

        _fmpz_vec_scalar_mul_fmpz(A, Q, Alen, x);
        _fmpz_vec_randtest(B, state, Blen, 100);
        _fmpz_vec_set(C, B, Blen);

        for (i = 0; i < r; i++)
            _fmpz_vec_scalar_divexact_fmpz(C + i * Bstride, A + i * Astride, c, x);

        if (alias)
        {
            /* expected: rows divided in place, gaps unchanged */
            fmpz * D = _fmpz_vec_init(Alen);
            _fmpz_vec_set(D, A, Alen);
            for (i = 0; i < r; i++)
                _fmpz_vec_scalar_divexact_fmpz(D + i * Astride, D + i * Astride, c, x);

            _fmpz_vec_scalar_divexact_fmpz_strided(A, Astride, A, Astride, r, c, x);

            if (!_fmpz_vec_equal(A, D, Alen))
                TEST_FUNCTION_FAIL("aliasing: r = %wd, c = %wd, Astride = %wd\n"
                    "x = %{fmpz}\nA = %{fmpz*}\nD = %{fmpz*}\n",
                    r, c, Astride, x, A, Alen, D, Alen);

            _fmpz_vec_clear(D, Alen);
        }
        else
        {
            _fmpz_vec_scalar_divexact_fmpz_strided(B, Bstride, A, Astride, r, c, x);

            if (!_fmpz_vec_equal(B, C, Blen))
                TEST_FUNCTION_FAIL("r = %wd, c = %wd, Astride = %wd, Bstride = %wd\n"
                    "x = %{fmpz}\nB = %{fmpz*}\nC = %{fmpz*}\n",
                    r, c, Astride, Bstride, x, B, Blen, C, Blen);
        }

        _fmpz_vec_clear(Q, Alen);
        _fmpz_vec_clear(A, Alen);
        _fmpz_vec_clear(B, Blen);
        _fmpz_vec_clear(C, Blen);
        fmpz_clear(x);
    }

    TEST_FUNCTION_END(state);
}
