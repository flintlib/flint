/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acf.h"
#include "acb.h"

TEST_FUNCTION_START(acf_approx_dot, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acf_ptr x, y;
        acf_t s1, s2, z;
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

        x = _acf_vec_init(len);
        y = _acf_vec_init(len);
        acf_init(s1);
        acf_init(s2);
        acf_init(z);

        switch (n_randint(state, 3))
        {
            case 0:
                for (i = 0; i < len; i++)
                {
                    acf_randtest(x + i, state, xbits, ebits);
                    acf_randtest(y + i, state, ybits, ebits);
                }
                break;

            /* Test with cancellation */
            case 1:
                for (i = 0; i < len; i++)
                {
                    if (i <= len / 2)
                    {
                        acf_randtest(x + i, state, xbits, ebits);
                        acf_randtest(y + i, state, ybits, ebits);
                    }
                    else
                    {
                        acf_neg(x + i, x + len - i - 1);
                        acf_set(y + i, y + len - i - 1);
                    }
                }
                break;

            default:
                for (i = 0; i < len; i++)
                {
                    if (i <= len / 2)
                    {
                        acf_randtest(x + i, state, xbits, ebits);
                        acf_randtest(y + i, state, ybits, ebits);
                    }
                    else
                    {
                        acf_neg_round(x + i, x + len - i - 1, 2 + n_randint(state, 500), ARF_RND_NEAR);
                        acf_set_round(y + i, y + len - i - 1, 2 + n_randint(state, 500), ARF_RND_NEAR);
                    }
                }
                break;
        }

        acf_randtest(s1, state, 200, 100);
        acf_randtest(s2, state, 200, 100);
        acf_randtest(z, state, xbits, ebits);

        acf_approx_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec, ARB_RND);

        /* With the fast algorithm, we expect identical results when
           reversing the vectors. */
        if (ebits <= 12 && 0)
        {
            acf_approx_dot(s2, initial ? z : NULL, subtract,
                !revx ? (x + len - 1) : x, !revx ? -1 : 1,
                !revy ? (y + len - 1) : y, !revy ? -1 : 1,
                len, prec, ARB_RND);

            if (!acf_equal(s1, s2))
            {
                flint_printf("FAIL (reversal)\n\n");
                flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd\n\n", iter, len, prec, ebits);

                if (initial)
                {
                    flint_printf("z = ", i); acf_printd(z, 100); flint_printf(" (%wd)\n\n", acf_bits(z));
                }

                for (i = 0; i < len; i++)
                {
                    flint_printf("x[%wd] = ", i); acf_printd(x + i, 100); flint_printf(" (%wd)\n", acf_bits(x + i));
                    flint_printf("y[%wd] = ", i); acf_printd(y + i, 100); flint_printf(" (%wd)\n", acf_bits(y + i));
                }
                flint_printf("\n\n");
                flint_printf("s1 = "); acf_printd(s1, 100); flint_printf("\n\n");
                flint_printf("s2 = "); acf_printd(s2, 100); flint_printf("\n\n");
                flint_abort();
            }
        }

        /* Compare with acb_dot */
        {
            acb_ptr ax, ay;
            acb_t as2, az;

            ax = _acb_vec_init(len);
            ay = _acb_vec_init(len);
            acb_init(as2);
            acb_init(az);

            for (i = 0; i < len; i++)
            {
                arb_set_arf(acb_realref(ax + i), acf_realref(x + i));
                arb_set_arf(acb_imagref(ax + i), acf_imagref(x + i));
                arb_set_arf(acb_realref(ay + i), acf_realref(y + i));
                arb_set_arf(acb_imagref(ay + i), acf_imagref(y + i));
            }

            if (initial)
            {
                arb_set_arf(acb_realref(az), acf_realref(z));
                arb_set_arf(acb_imagref(az), acf_imagref(z));
            }

            acb_dot(as2, initial ? az : NULL, subtract,
                revx ? (ax + len - 1) : ax, revx ? -1 : 1,
                revy ? (ay + len - 1) : ay, revy ? -1 : 1,
                len, prec);

            {
                mag_t err, xx, yy;

                mag_init(err);
                mag_init(xx);
                mag_init(yy);

                if (initial)
                    acf_get_mag(err, z);

                for (i = 0; i < len; i++)
                {
                    acf_get_mag(xx, revx ? x + len - 1 - i : x + i);
                    acf_get_mag(yy, revx ? y + len - 1 - i : y + i);
                    mag_addmul(err, xx, yy);
                }

                mag_mul_2exp_si(err, err, -prec + 2);
                acb_add_error_mag(as2, err);

                if (!arb_contains_arf(acb_realref(as2), acf_realref(s1)) || !arb_contains_arf(acb_imagref(as2), acf_imagref(s1)))
                {
                    flint_printf("FAIL (inclusion)\n\n");
                    flint_printf("iter = %wd, len = %wd, prec = %wd, ebits = %wd\n\n", iter, len, prec, ebits);

                    if (initial)
                    {
                        flint_printf("z = ", i); acf_printd(z, 100); flint_printf(" (%wd)\n\n", acf_bits(z));
                    }

                    for (i = 0; i < len; i++)
                    {
                        flint_printf("x[%wd] = ", i); acf_printd(x + i, 100); flint_printf(" (%wd)\n", acf_bits(x + i));
                        flint_printf("y[%wd] = ", i); acf_printd(y + i, 100); flint_printf(" (%wd)\n", acf_bits(y + i));
                    }
                    flint_printf("\n\n");
                    flint_printf("s1 = "); acf_printd(s1, 100); flint_printf("\n\n");
                    flint_printf("s2 = "); acb_printd(as2, 100); flint_printf("\n\n");
                    flint_abort();
                }

                mag_clear(err);
                mag_clear(xx);
                mag_clear(yy);
            }

            _acb_vec_clear(ax, len);
            _acb_vec_clear(ay, len);
            acb_clear(as2);
            acb_clear(az);
        }

        acf_clear(s1);
        acf_clear(s2);
        acf_clear(z);
        _acf_vec_clear(x, len);
        _acf_vec_clear(y, len);
    }

    TEST_FUNCTION_END(state);
}
