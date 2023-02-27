/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/double_extras.h"
#include "mag.h"

/* XXX: d_randtest is not good enough */

#define EXP_MINUS_32 2.3283064365386962891e-10
#define EXP_MINUS_64 5.42101086242752217e-20

double
d_randtest2(flint_rand_t state)
{
    mp_limb_t m1, m2;
    double t;

    if (FLINT_BITS == 64)
    {
        m1 = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));

        t = ((double) m1) * EXP_MINUS_64;
    }
    else
    {
        m1 = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));
        m2 = n_randtest(state);

        t = ((double) m1) * EXP_MINUS_32 +
            ((double) m2) * EXP_MINUS_64;
    }

    return t;
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("d_log_upper_bound....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        mpfr_t t;
        double x, y, z;

        mpfr_init2(t, 53);

        x = d_randtest2(state);
        x = ldexp(x, 100 - n_randint(state, 200));

        switch (n_randint(state, 10))
        {
            case 0:
                x = 1.0 + x;
                break;
            case 1:
                x = 1.0 - x;
                break;
            case 2:
                x = D_INF;
                break;
            case 3:
                x = 0.0;
                break;
            case 4:
                x = D_NAN;
                break;
            default:
                break;
        }

        y = mag_d_log_upper_bound(x);

        mpfr_set_d(t, x, MPFR_RNDD);
        mpfr_log(t, t, MPFR_RNDU);
        z = mpfr_get_d(t, MPFR_RNDD);

        if (y < z || fabs(y-z) > 0.000001 * fabs(z))
        {
            flint_printf("FAIL\n");
            flint_printf("x = %.20g\n", x);
            flint_printf("y = %.20g\n", y);
            flint_printf("z = %.20g\n", z);
            flint_abort();
        }

        mpfr_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

