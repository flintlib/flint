/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_dot_general, state)
{
    int iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        fmpz * a, * b;
        fmpz_t s, t, c;
        slong n, bits1, bits2;
        int alias, negate, initial, reverse;

        initial = n_randint(state, 2);
        alias = n_randint(state, 2);
        negate = n_randint(state, 2);
        reverse = n_randint(state, 2);
        n = n_randint(state, 8);

        if (n_randint(state, 30) == 0)
            bits1 = 2 + n_randint(state, 20000);
        else
            bits1 = 2 + n_randint(state, 1000);

        bits2 = 2 + n_randint(state, 1000);

        a = _fmpz_vec_init(n);
        b = _fmpz_vec_init(n);
        fmpz_init(s);
        fmpz_init(t);
        fmpz_init(c);

        fmpz_randtest(c, state, bits1);
        _fmpz_vec_randtest(a, state, n, bits1);
        _fmpz_vec_randtest(b, state, n, bits2);

        if (initial && alias)
        {
            fmpz_set(s, c);
            fmpz_set(t, c);
            _fmpz_vec_dot_general(s, s, negate, a, b, reverse, n);
            _fmpz_vec_dot_general_naive(t, t, negate, a, b, reverse, n);
        }
        else
        {
            _fmpz_vec_dot_general(s, initial ? c : NULL, negate, a, b, reverse, n);
            _fmpz_vec_dot_general_naive(t, initial ? c : NULL, negate, a, b, reverse, n);
        }

        if (!fmpz_equal(s, t) || !_fmpz_is_canonical(s))
        {
            flint_printf("negate = %d, initial = %d, reverse = %d, alias = %d\n", negate, initial, reverse, alias);
            flint_printf("c = %{fmpz}\n\n", c);
            flint_printf("a = %{fmpz*}\n\n", a, n);
            flint_printf("b = %{fmpz*}\n\n", b, n);
            flint_printf("s = %{fmpz}\n\n", s);
            flint_printf("t = %{fmpz}\n\n", t);
            flint_abort();
        }

        fmpz_clear(s);
        fmpz_clear(t);
        fmpz_clear(c);
        _fmpz_vec_clear(a, n);
        _fmpz_vec_clear(b, n);
    }

    TEST_FUNCTION_END(state);
}
