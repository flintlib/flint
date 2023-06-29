/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("eld_interval....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        int a = n_randint(state, 2);
        slong prec = ACB_THETA_ELD_DEFAULT_PREC;
        slong mag_bits = n_randint(state, 5);
        int guaranteed_pt = iter % 2;

        slong min, max, mid;
        arb_t ctr, tmax, tmin;
        arf_t rad;
        arf_t pos;
        int fail;

        arb_init(ctr);
        arf_init(rad);
        arb_init(tmax);
        arb_init(tmin);
        arf_init(pos);

        arb_randtest(ctr, state, prec, mag_bits);
        arf_randtest(rad, state, prec, mag_bits);

        if (guaranteed_pt)
        {
            arf_set_si(arb_midref(ctr), a);
            arf_abs(rad, rad);
        }

        acb_theta_eld_interval(&min, &mid, &max, ctr, rad, a, prec);
        arb_set_si(tmax, max + 3);
        arb_sub_arf(tmax, tmax, rad, prec);
        arb_set_si(tmin, min - 3);
        arb_add_arf(tmin, tmin, rad, prec);

        fail = ((min > max) && guaranteed_pt)
            || ((min <= max) &&
                (FLINT_ABS(min) % 2 != a
                 || FLINT_ABS(mid) % 2 != a
                 || FLINT_ABS(max) % 2 != a
                 || mid < min
                 || mid > max || !arb_gt(tmax, ctr) || !arb_lt(tmin, ctr)));

        if (fail)
        {
            flint_printf("FAIL\n");
            flint_printf("Center, radius:\n");
            arb_printd(ctr, 10);
            flint_printf("\n");
            arf_printd(rad, 10);
            flint_printf("\n");
            flint_printf("a = %i, min = %wd, mid = %wd, max = %wd\n", a, min,
                         mid, max);
            fflush(stdout);
            flint_abort();
        }

        arb_clear(ctr);
        arf_clear(rad);
        arb_clear(tmax);
        arb_clear(tmin);
        arf_clear(pos);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
