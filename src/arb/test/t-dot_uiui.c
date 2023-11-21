/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_dot_uiui, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_ptr x, y;
        ulong * w;
        arb_t s1, s2, z;
        fmpz_t c;
        slong i, len, prec;
        int initial, subtract, revx, revy;

        len = n_randint(state, 5);
        prec = 2 + n_randint(state, 200);

        initial = n_randint(state, 2);
        subtract = n_randint(state, 2);
        revx = n_randint(state, 2);
        revy = n_randint(state, 2);

        x = _arb_vec_init(len);
        y = _arb_vec_init(len);
        w = flint_malloc(2 * sizeof(ulong) * len);
        arb_init(s1);
        arb_init(s2);
        arb_init(z);
        fmpz_init(c);

        for (i = 0; i < len; i++)
        {
            arb_randtest(x + i, state, 2 + n_randint(state, 200), 10);
            w[2 * i] = n_randtest(state);
            w[2 * i + 1] = n_randtest(state);
            fmpz_set_uiui(c, w[2 * i + 1], w[2 * i]);
            arb_set_fmpz(y + i, c);
        }

        arb_randtest(s1, state, 200, 10);
        arb_randtest(s2, state, 200, 10);
        arb_randtest(z, state, 200, 10);

        arb_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);

        arb_dot_uiui(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (w + 2 * len - 2) : w, revy ? -1 : 1,
            len, prec);

        if (!arb_overlaps(s1, s2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("iter = %wd, len = %wd, prec = %wd\n\n", iter, len, prec);

            if (initial)
            {
                flint_printf("z = ", i); arb_printn(z, 100, ARB_STR_MORE); flint_printf(" (%wd)\n\n", arb_bits(z));
            }

            for (i = 0; i < len; i++)
            {
                flint_printf("x[%wd] = ", i); arb_printn(x + i, 100, ARB_STR_MORE); flint_printf(" (%wd)\n", arb_bits(x + i));
                flint_printf("y[%wd] = ", i); arb_printn(y + i, 100, ARB_STR_MORE); flint_printf(" (%wd)\n", arb_bits(y + i));
            }
            flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 100, ARB_STR_MORE); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 100, ARB_STR_MORE); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(s1);
        arb_clear(s2);
        arb_clear(z);
        _arb_vec_clear(x, len);
        _arb_vec_clear(y, len);
        flint_free(w);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
