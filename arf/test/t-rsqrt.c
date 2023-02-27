/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_rsqrt_naive(arf_t z, const arf_t x, slong prec, arf_rnd_t rnd)
{
    fmpr_t a;
    slong r;

    fmpr_init(a);

    arf_get_fmpr(a, x);

    r = fmpr_rsqrt(a, a, prec, rnd);

    arf_set_fmpr(z, a);

    fmpr_clear(a);

    return (r == FMPR_RESULT_EXACT) ? 0 : 1;
}

int main()
{
    slong iter, iter2;
    flint_rand_t state;

    flint_printf("rsqrt....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arf_t x, z, v;
        slong prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(z);
        arf_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 100);
            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 4))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                default: rnd = ARF_RND_CEIL; break;
            }

            switch (n_randint(state, 2))
            {
            case 0:
                r1 = arf_rsqrt(z, x, prec, rnd);
                r2 = arf_rsqrt_naive(v, x, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                r2 = arf_rsqrt_naive(v, x, prec, rnd);
                r1 = arf_rsqrt(x, x, prec, rnd);
                if (!arf_equal(v, x) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;
            }
        }

        arf_clear(x);
        arf_clear(z);
        arf_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
