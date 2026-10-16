/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"
#include "ball_helpers.h"

/* the series constants pi and log 2 at up to 1000 limbs against arb,
   at a radius within a couple of guard limbs */

static void
test_log2(void)
{
    slong ns[] = { 2, 3, 5, 16, 100, 500 };
    slong i;

    for (i = 0; i < 6; i++)
    {
        slong n = ns[i];
        mp_real_t v;
        arb_t p1, p2;

        mp_real_init(v);
        arb_init(p1);
        arb_init(p2);

        mp_real_const_log2(v, n, 0);
        mp_real_get_arb(p1, v);
        arb_const_log2(p2, FLINT_BITS * n + 64);

        if (!arb_overlaps(p1, p2))
        {
            flint_printf("FAIL: log2 n=%wd\n", n);
            mp_real_print(v);
            flint_abort();
        }

        {
            mag_t rad;
            mag_init(rad);
            mag_set_ui(rad, FLINT_MAX(v->err, 1));
            mag_mul_2exp_si(rad, rad,
                FLINT_BITS * (v->exp - v->size));
            if (mag_cmp_2exp_si(rad, -FLINT_BITS * (n - 4)) > 0)
            {
                flint_printf("FAIL: log2 radius too large n=%wd "
                    "err=%wu size=%wd\n", n, v->err, v->size);
                flint_abort();
            }
            mag_clear(rad);
        }

        mp_real_clear(v);
        arb_clear(p1);
        arb_clear(p2);
    }
}

static void
test_pi(void)
{
    slong ns[] = {2, 3, 4, 5, 8, 16, 33, 100, 331, 1000};
    slong i;

    for (i = 0; i < 10; i++)
    {
        slong n = ns[i];
        mp_real_t pi;
        arb_t p1, p2;

        mp_real_init(pi);
        arb_init(p1);
        arb_init(p2);

        mp_real_const_pi4(pi, n, 0);
        mp_real_mul_2exp_si(pi, pi, 2);
        mp_real_get_arb(p1, pi);
        arb_const_pi(p2, FLINT_BITS * n + 64);

        if (!arb_overlaps(p1, p2))
        {
            flint_printf("FAIL: pi n=%wd\n", n);
            mp_real_print(pi);
            flint_abort();
        }

        /* the radius should be within a couple of guard limbs */
        {
            mag_t rad;
            mag_init(rad);
            mag_set_ui(rad, FLINT_MAX(pi->err, 1));
            mag_mul_2exp_si(rad, rad,
                FLINT_BITS * (pi->exp - pi->size));
            if (mag_cmp_2exp_si(rad, -FLINT_BITS * (n - 4)) > 0)
            {
                flint_printf("FAIL: pi radius too large n=%wd err=%wu "
                    "size=%wd\n", n, pi->err, pi->size);
                flint_abort();
            }
            mag_clear(rad);
        }

        mp_real_clear(pi);
        arb_clear(p1);
        arb_clear(p2);
    }
}

TEST_FUNCTION_START(mp_real_const_series, state)
{
    test_pi();
    test_log2();

    TEST_FUNCTION_END(state);
}
