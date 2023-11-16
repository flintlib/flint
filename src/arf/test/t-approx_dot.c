/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "arb.h"

TEST_FUNCTION_START(arf_approx_dot, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_ptr x, y;
        arf_t s1, s2, z;
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

        x = _arf_vec_init(len);
        y = _arf_vec_init(len);
        arf_init(s1);
        arf_init(s2);
        arf_init(z);

        switch (n_randint(state, 3))
        {
            case 0:
                for (i = 0; i < len; i++)
                {
                    arf_randtest(x + i, state, xbits, ebits);
                    arf_randtest(y + i, state, ybits, ebits);
                }
                break;

            /* Test with cancellation */
            case 1:
                for (i = 0; i < len; i++)
                {
                    if (i <= len / 2)
                    {
                        arf_randtest(x + i, state, xbits, ebits);
                        arf_randtest(y + i, state, ybits, ebits);
                    }
                    else
                    {
                        arf_neg(x + i, x + len - i - 1);
                        arf_set(y + i, y + len - i - 1);
                    }
                }
                break;

            default:
                for (i = 0; i < len; i++)
                {
                    if (i <= len / 2)
                    {
                        arf_randtest(x + i, state, xbits, ebits);
                        arf_randtest(y + i, state, ybits, ebits);
                    }
                    else
                    {
                        arf_neg_round(x + i, x + len - i - 1, 2 + n_randint(state, 500), ARF_RND_NEAR);
                        arf_set_round(y + i, y + len - i - 1, 2 + n_randint(state, 500), ARF_RND_NEAR);
                    }
                }
                break;
        }

        arf_randtest(s1, state, 200, 100);
        arf_randtest(s2, state, 200, 100);
        arf_randtest(z, state, xbits, ebits);

        arf_approx_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec, ARB_RND);

        /* With the fast algorithm, we expect identical results when
           reversing the vectors. */
        if (ebits <= 12)
        {
            arf_approx_dot(s2, initial ? z : NULL, subtract,
                !revx ? (x + len - 1) : x, !revx ? -1 : 1,
                !revy ? (y + len - 1) : y, !revy ? -1 : 1,
                len, prec, ARB_RND);

            if (!arf_equal(s1, s2))
            {
                flint_printf("FAIL (reversal)\n\n");
                flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd\n\n", iter, len, prec, ebits);

                if (initial)
                {
                    flint_printf("z = ", i); arf_printd(z, 100); flint_printf(" (%wd)\n\n", arf_bits(z));
                }

                for (i = 0; i < len; i++)
                {
                    flint_printf("x[%wd] = ", i); arf_printd(x + i, 100); flint_printf(" (%wd)\n", arf_bits(x + i));
                    flint_printf("y[%wd] = ", i); arf_printd(y + i, 100); flint_printf(" (%wd)\n", arf_bits(y + i));
                }
                flint_printf("\n\n");
                flint_printf("s1 = "); arf_printd(s1, 100); flint_printf("\n\n");
                flint_printf("s2 = "); arf_printd(s2, 100); flint_printf("\n\n");
                flint_abort();
            }
        }

        /* Compare with arb_dot */
        {
            arb_ptr ax, ay;
            arb_t as2, az;

            ax = _arb_vec_init(len);
            ay = _arb_vec_init(len);
            arb_init(as2);
            arb_init(az);

            for (i = 0; i < len; i++)
            {
                arb_set_arf(ax + i, x + i);
                arb_set_arf(ay + i, y + i);
            }

            if (initial)
                arb_set_arf(az, z);

            arb_dot(as2, initial ? az : NULL, subtract,
                revx ? (ax + len - 1) : ax, revx ? -1 : 1,
                revy ? (ay + len - 1) : ay, revy ? -1 : 1,
                len, prec);

            {
                mag_t err, xx, yy;

                mag_init(err);
                mag_init(xx);
                mag_init(yy);

                if (initial)
                    arf_get_mag(err, z);

                for (i = 0; i < len; i++)
                {
                    arf_get_mag(xx, revx ? x + len - 1 - i : x + i);
                    arf_get_mag(yy, revx ? y + len - 1 - i : y + i);
                    mag_addmul(err, xx, yy);
                }

                mag_mul_2exp_si(err, err, -prec + 2);
                arb_add_error_mag(as2, err);

                if (!arb_contains_arf(as2, s1))
                {
                    flint_printf("FAIL (inclusion)\n\n");
                    flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd\n\n", iter, len, prec, ebits);

                    if (initial)
                    {
                        flint_printf("z = ", i); arf_printd(z, 100); flint_printf(" (%wd)\n\n", arf_bits(z));
                    }

                    for (i = 0; i < len; i++)
                    {
                        flint_printf("x[%wd] = ", i); arf_printd(x + i, 100); flint_printf(" (%wd)\n", arf_bits(x + i));
                        flint_printf("y[%wd] = ", i); arf_printd(y + i, 100); flint_printf(" (%wd)\n", arf_bits(y + i));
                    }
                    flint_printf("\n\n");
                    flint_printf("s1 = "); arf_printd(s1, 100); flint_printf("\n\n");
                    flint_printf("s2 = "); arf_printd(s2, 100); flint_printf("\n\n");
                    flint_abort();
                }

                mag_clear(err);
                mag_clear(xx);
                mag_clear(yy);
            }

            _arb_vec_clear(ax, len);
            _arb_vec_clear(ay, len);
            arb_clear(as2);
            arb_clear(az);
        }

        arf_clear(s1);
        arf_clear(s2);
        arf_clear(z);
        _arf_vec_clear(x, len);
        _arf_vec_clear(y, len);
    }

    TEST_FUNCTION_END(state);
}
