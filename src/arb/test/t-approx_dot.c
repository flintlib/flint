/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_approx_dot, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_ptr x, y;
        arb_t s1, s2, z;
        slong i, len, prec, xbits, ybits, ebits;
        int initial, subtract, revx, revy;

        if (n_randint(state, 100) == 0)
            len = n_randint(state, 100);
        else if (n_randint(state, 10) == 0)
            len = n_randint(state, 10);
        else
            len = n_randint(state, 3);

        if (n_randint(state, 10) != 0 || len > 10)
        {
            prec = 2 + n_randint(state, 500);
            xbits = 2 + n_randint(state, 500);
            ybits = 2 + n_randint(state, 500);
        }
        else
        {
            prec = 2 + n_randint(state, 4000);
            xbits = 2 + n_randint(state, 4000);
            ybits = 2 + n_randint(state, 4000);
        }

        if (n_randint(state, 100) == 0)
            ebits = 2 + n_randint(state, 100);
        else
            ebits = 2 + n_randint(state, 10);

        initial = n_randint(state, 2);
        subtract = n_randint(state, 2);
        revx = n_randint(state, 2);
        revy = n_randint(state, 2);

        x = _arb_vec_init(len);
        y = _arb_vec_init(len);
        arb_init(s1);
        arb_init(s2);
        arb_init(z);

        switch (n_randint(state, 3))
        {
            case 0:
                for (i = 0; i < len; i++)
                {
                    arb_randtest(x + i, state, xbits, ebits);
                    arb_randtest(y + i, state, ybits, ebits);
                }
                break;

            /* Test with cancellation */
            case 1:
                for (i = 0; i < len; i++)
                {
                    if (i <= len / 2)
                    {
                        arb_randtest(x + i, state, xbits, ebits);
                        arb_randtest(y + i, state, ybits, ebits);
                    }
                    else
                    {
                        arb_neg(x + i, x + len - i - 1);
                        arb_set(y + i, y + len - i - 1);
                    }
                }
                break;

            default:
                for (i = 0; i < len; i++)
                {
                    if (i <= len / 2)
                    {
                        arb_randtest(x + i, state, xbits, ebits);
                        arb_randtest(y + i, state, ybits, ebits);
                    }
                    else
                    {
                        arb_neg_round(x + i, x + len - i - 1, 2 + n_randint(state, 500));
                        arb_set_round(y + i, y + len - i - 1, 2 + n_randint(state, 500));
                    }
                }
                break;
        }

        arb_randtest(s1, state, 200, 100);
        arb_randtest(s2, state, 200, 100);
        arb_randtest(z, state, xbits, ebits);

        arb_approx_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);
        mag_zero(arb_radref(s1));

        /* With the fast algorithm, we expect identical results when
           reversing the vectors. */
        if (ebits <= 12)
        {
            arb_approx_dot(s2, initial ? z : NULL, subtract,
                !revx ? (x + len - 1) : x, !revx ? -1 : 1,
                !revy ? (y + len - 1) : y, !revy ? -1 : 1,
                len, prec);
            mag_zero(arb_radref(s2));

            if (!arb_equal(s1, s2))
            {
                flint_printf("FAIL (reversal)\n\n");
                flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd\n\n", iter, len, prec, ebits);

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
        }

        /* Verify that radii are ignored */
        for (i = 0; i < len; i++)
        {
            arb_get_mid_arb(x + i, x + i);
            arb_get_mid_arb(y + i, y + i);
        }
        arb_get_mid_arb(z, z);

        arb_approx_dot(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);
        mag_zero(arb_radref(s2));

        if (!arb_equal(s1, s2))
        {
            flint_printf("FAIL (radii)\n\n");
            flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd\n\n", iter, len, prec, ebits);

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

        /* Compare with arb_dot */
        arb_approx_dot(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);

        {
            mag_t err, xx, yy;

            mag_init(err);
            mag_init(xx);
            mag_init(yy);

            if (initial)
                arb_get_mag(err, z);

            for (i = 0; i < len; i++)
            {
                arb_get_mag(xx, revx ? x + len - 1 - i : x + i);
                arb_get_mag(yy, revx ? y + len - 1 - i : y + i);
                mag_addmul(err, xx, yy);
            }

            mag_mul_2exp_si(err, err, -prec + 2);
            arb_add_error_mag(s2, err);

            if (!arb_contains(s2, s1))
            {
                flint_printf("FAIL (inclusion)\n\n");
                flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd\n\n", iter, len, prec, ebits);

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

            mag_clear(err);
            mag_clear(xx);
            mag_clear(yy);
        }

        arb_clear(s1);
        arb_clear(s2);
        arb_clear(z);
        _arb_vec_clear(x, len);
        _arb_vec_clear(y, len);
    }

    TEST_FUNCTION_END(state);
}
