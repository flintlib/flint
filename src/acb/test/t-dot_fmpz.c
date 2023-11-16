/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "acb.h"

TEST_FUNCTION_START(acb_dot_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_ptr x, y;
        fmpz * w;
        acb_t s1, s2, z;
        slong i, len, prec;
        int initial, subtract, revx, revy;

        len = n_randint(state, 5);
        prec = 2 + n_randint(state, 200);

        initial = n_randint(state, 2);
        subtract = n_randint(state, 2);
        revx = n_randint(state, 2);
        revy = n_randint(state, 2);

        x = _acb_vec_init(len);
        y = _acb_vec_init(len);
        w = _fmpz_vec_init(len);
        acb_init(s1);
        acb_init(s2);
        acb_init(z);

        for (i = 0; i < len; i++)
        {
            acb_randtest(x + i, state, 2 + n_randint(state, 200), 10);
            fmpz_randtest(w + i, state, 1 + n_randint(state, 200));
            acb_set_fmpz(y + i, w + i);
        }

        acb_randtest(s1, state, 200, 10);
        acb_randtest(s2, state, 200, 10);
        acb_randtest(z, state, 200, 10);

        acb_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);

        acb_dot_fmpz(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (w + len - 1) : w, revy ? -1 : 1,
            len, prec);

        if (!acb_overlaps(s1, s2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("iter = %wd, len = %wd, prec = %wd\n\n", iter, len, prec);

            if (initial)
            {
                flint_printf("z = ", i); acb_printn(z, 100, ARB_STR_MORE); flint_printf(" (%wd)\n\n", acb_bits(z));
            }

            for (i = 0; i < len; i++)
            {
                flint_printf("x[%wd] = ", i); acb_printn(x + i, 100, ARB_STR_MORE); flint_printf(" (%wd)\n", acb_bits(x + i));
                flint_printf("y[%wd] = ", i); acb_printn(y + i, 100, ARB_STR_MORE); flint_printf(" (%wd)\n", acb_bits(y + i));
            }
            flint_printf("\n\n");
            flint_printf("s1 = "); acb_printn(s1, 100, ARB_STR_MORE); flint_printf("\n\n");
            flint_printf("s2 = "); acb_printn(s2, 100, ARB_STR_MORE); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(s1);
        acb_clear(s2);
        acb_clear(z);
        _acb_vec_clear(x, len);
        _acb_vec_clear(y, len);
        _fmpz_vec_clear(w, len);
    }

    TEST_FUNCTION_END(state);
}
