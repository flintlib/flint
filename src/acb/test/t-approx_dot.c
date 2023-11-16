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

TEST_FUNCTION_START(acb_approx_dot, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_ptr x, y;
        acb_t s1, s2, z;
        slong i, len, prec, xbits, ybits, ebits;
        int initial, subtract, revx, revy;

        if (n_randint(state, 100) == 0)
            len = n_randint(state, 100);
        else if (n_randint(state, 10) == 0)
            len = n_randint(state, 10);
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

        acb_approx_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);
        mag_zero(arb_radref(acb_realref(s1)));
        mag_zero(arb_radref(acb_imagref(s1)));

        /* With the fast algorithm, we expect identical results when
           reversing the vectors. */
        if (ebits <= 12)
        {
            acb_approx_dot(s2, initial ? z : NULL, subtract,
                !revx ? (x + len - 1) : x, !revx ? -1 : 1,
                !revy ? (y + len - 1) : y, !revy ? -1 : 1,
                len, prec);
            mag_zero(arb_radref(acb_realref(s2)));
            mag_zero(arb_radref(acb_imagref(s2)));

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

        /* Verify that radii are ignored */
        for (i = 0; i < len; i++)
        {
            arb_get_mid_arb(acb_realref(x + i), acb_realref(x + i));
            arb_get_mid_arb(acb_imagref(x + i), acb_imagref(x + i));
            arb_get_mid_arb(acb_realref(y + i), acb_realref(y + i));
            arb_get_mid_arb(acb_imagref(y + i), acb_imagref(y + i));
        }
        arb_get_mid_arb(acb_realref(z), acb_realref(z));
        arb_get_mid_arb(acb_imagref(z), acb_imagref(z));

        acb_approx_dot(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);
        mag_zero(arb_radref(acb_realref(s2)));
        mag_zero(arb_radref(acb_imagref(s2)));

        if (!acb_equal(s1, s2))
        {
            flint_printf("FAIL (radii)\n\n");
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

        /* Compare with acb_dot */
        acb_approx_dot(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);

        {
            mag_t err, xx, yy;

            mag_init(err);
            mag_init(xx);
            mag_init(yy);

            if (initial)
                acb_get_mag(err, z);

            for (i = 0; i < len; i++)
            {
                acb_get_mag(xx, revx ? x + len - 1 - i : x + i);
                acb_get_mag(yy, revx ? y + len - 1 - i : y + i);
                mag_addmul(err, xx, yy);
            }

            mag_mul_2exp_si(err, err, -prec + 2);
            acb_add_error_mag(s2, err);

            if (!acb_contains(s2, s1))
            {
                flint_printf("FAIL (inclusion)\n\n");
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

            mag_clear(err);
            mag_clear(xx);
            mag_clear(yy);
        }

        acb_clear(s1);
        acb_clear(s2);
        acb_clear(z);
        _acb_vec_clear(x, len);
        _acb_vec_clear(y, len);
    }

    TEST_FUNCTION_END(state);
}
