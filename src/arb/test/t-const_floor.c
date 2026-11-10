/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "mpn_extras.h"
#include "mp_real/impl.h"
#include "arb/elem.h"

/* _arb_const_set_floor gives the radius half an ulp, which is valid only
   when the midpoint is the correctly rounded value: the tabulated floor Y
   (of w limbs, as the function reads it) and Y + 1 must round alike at
   every precision it serves.  Also compares with the constants at a
   higher precision. */
TEST_FUNCTION_START(arb_const_floor, state)
{
    static const int ids[2] = { MP_REAL_CONST_ID_PI4, MP_REAL_CONST_ID_LOG2 };
    slong prec, maxprec = FLINT_BITS * (MP_REAL_CONST_STATIC_N - 1) - 16;
    int t;

    for (t = 0; t < 2; t++)
    {
        for (prec = 2; prec <= maxprec; prec++)
        {
            slong w = (prec + 16 + FLINT_BITS - 1) / FLINT_BITS, f1, f2;
            ulong y[MP_REAL_CONST_STATIC_N + 1];
            arf_t a, b;
            int carry;

            arf_init(a);
            arf_init(b);

            flint_mpn_copyi(y, _mp_real_const_ptr(ids[t], w), w);
            _arf_set_round_mpn(a, &f1, y, w, 0, prec, ARF_RND_NEAR);
            y[w] = mpn_add_1(y, y, w, 1);
            carry = (y[w] != 0);
            _arf_set_round_mpn(b, &f2, y, w + carry, 0, prec, ARF_RND_NEAR);
            fmpz_set_si(ARF_EXPREF(a), f1);
            fmpz_set_si(ARF_EXPREF(b), f2);

            if (!arf_equal(a, b))
            {
                flint_printf("FAIL: the table %d does not round correctly at prec = %wd\n", t, prec);
                flint_abort();
            }

            arf_clear(a);
            arf_clear(b);
        }
    }

    /* the balls against higher-precision values */
    for (prec = 2; prec <= maxprec + 200; prec += 1 + n_randint(state, 30))
    {
        arb_t x, y;
        arb_init(x);
        arb_init(y);

        for (t = 0; t < 2; t++)
        {
            if (t == 0)
            {
                arb_const_pi(x, prec);
                arb_const_pi(y, prec + 100);
            }
            else
            {
                arb_const_log2(x, prec);
                arb_const_log2(y, prec + 100);
            }

            if (!arb_contains(x, y))
            {
                flint_printf("FAIL: containment, t = %d, prec = %wd\n", t, prec);
                flint_abort();
            }

            /* the radius at most half an ulp plus 2^-30 ulp */
            {
                mag_t u;
                mag_init(u);
                arf_get_mag(u, arb_midref(x));
                mag_mul_2exp_si(u, u, -prec);
                mag_mul_ui(u, u, (UWORD(1) << 29) + 1);
                mag_mul_2exp_si(u, u, -29);
                if (mag_cmp(arb_radref(x), u) > 0)
                {
                    flint_printf("FAIL: radius, t = %d, prec = %wd\n", t, prec);
                    arb_printd(x, 30); flint_printf("\n");
                    flint_abort();
                }
                mag_clear(u);
            }
        }

        arb_clear(x);
        arb_clear(y);
    }

    TEST_FUNCTION_END(state);
}
