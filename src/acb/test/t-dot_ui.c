/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("dot_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        acb_ptr x, y;
        ulong * w;
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
        w = flint_malloc(sizeof(ulong) * len);
        acb_init(s1);
        acb_init(s2);
        acb_init(z);

        for (i = 0; i < len; i++)
        {
            acb_randtest(x + i, state, 2 + n_randint(state, 200), 10);
            w[i] = n_randtest(state);
            acb_set_ui(y + i, w[i]);
        }

        acb_randtest(s1, state, 200, 10);
        acb_randtest(s2, state, 200, 10);
        acb_randtest(z, state, 200, 10);

        acb_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);

        acb_dot_ui(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (w + len - 1) : w, revy ? -1 : 1,
            len, prec);

        if (!acb_equal(s1, s2))
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
        flint_free(w);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

