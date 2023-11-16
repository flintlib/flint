/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

FLINT_DLL extern slong acb_dot_gauss_dot_cutoff;

TEST_FUNCTION_START(acb_dot, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_ptr x, y;
        acb_t s1, s2, z;
        slong i, len, prec, xbits, ybits, ebits;
        int ok, initial, subtract, revx, revy;

        if (n_randint(state, 100) == 0)
            len = n_randint(state, 50);
        else if (n_randint(state, 10) == 0)
            len = n_randint(state, 5);
        else
            len = n_randint(state, 3);

        acb_dot_gauss_dot_cutoff = 3 + n_randint(state, 30);

        if (n_randint(state, 10) != 0 || len > 10)
        {
            prec = 2 + n_randint(state, 500);
            xbits = 2 + n_randint(state, 500);
            ybits = 2 + n_randint(state, 500);
        }
        else
        {
            prec = 2 + n_randint(state, 5000);
            xbits = 2 + n_randint(state, 5000);
            ybits = 2 + n_randint(state, 5000);
        }

        if (n_randint(state, 100) == 0)
            ebits = 2 + n_randint(state, 100);
        else
            ebits = 2 + n_randint(state, 10);

        initial = n_randint(state, 2);
        subtract = n_randint(state, 2);
        revx = n_randint(state, 2);
        revy = n_randint(state, 2);

        x = _acb_vec_init(len);
        y = _acb_vec_init(len);
        acb_init(s1);
        acb_init(s2);
        acb_init(z);

        switch (n_randint(state, 3))
        {
            case 0:
                for (i = 0; i < len; i++)
                {
                    acb_randtest(x + i, state, xbits, ebits);
                    acb_randtest(y + i, state, ybits, ebits);
                }
                break;

            /* Test with cancellation */
            case 1:
                for (i = 0; i < len; i++)
                {
                    if (i <= len / 2)
                    {
                        acb_randtest(x + i, state, xbits, ebits);
                        acb_randtest(y + i, state, ybits, ebits);
                    }
                    else
                    {
                        acb_neg(x + i, x + len - i - 1);
                        acb_set(y + i, y + len - i - 1);
                    }
                }
                break;

            default:
                for (i = 0; i < len; i++)
                {
                    if (i <= len / 2)
                    {
                        acb_randtest(x + i, state, xbits, ebits);
                        acb_randtest(y + i, state, ybits, ebits);
                    }
                    else
                    {
                        acb_neg_round(x + i, x + len - i - 1, 2 + n_randint(state, 500));
                        acb_set_round(y + i, y + len - i - 1, 2 + n_randint(state, 500));
                    }
                }
                break;
        }

        acb_randtest(s1, state, 200, 100);
        acb_randtest(s2, state, 200, 100);
        acb_randtest(z, state, xbits, ebits);

        acb_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);

        acb_dot_precise(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, ebits <= 12 ? ARF_PREC_EXACT : 2 * prec + 100);

        if (ebits <= 12)
            ok = acb_contains(s1, s2);
        else
            ok = acb_overlaps(s1, s2);

        if (!ok)
        {
            flint_printf("FAIL\n\n");
            flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd, subtract = %d\n\n", iter, len, prec, ebits, subtract);

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

        /* With the fast algorithm, we expect identical results when
           reversing the vectors. */
        if (ebits <= 12)
        {
            revx ^= 1;
            revy ^= 1;

            acb_dot(s2, initial ? z : NULL, subtract,
                revx ? (x + len - 1) : x, revx ? -1 : 1,
                revy ? (y + len - 1) : y, revy ? -1 : 1,
                len, prec);

            if (!acb_equal(s1, s2))
            {
                flint_printf("FAIL (reversal)\n\n");
                flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd\n\n", iter, len, prec, ebits);

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
        }

        acb_clear(s1);
        acb_clear(s2);
        acb_clear(z);
        _acb_vec_clear(x, len);
        _acb_vec_clear(y, len);
    }

    TEST_FUNCTION_END(state);
}
