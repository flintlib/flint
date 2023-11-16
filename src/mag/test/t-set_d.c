/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "double_extras.h"
#include "mag.h"

/* Defined in t-d_log_lower_bound.c, t-d_log_upper_bound.c, t-set_d_2exp_fmpz.c
 * and t-set_d.c */
#ifndef d_randtest2
#define d_randtest2 d_randtest2
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
#endif

TEST_FUNCTION_START(mag_set_d, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t a, b, c;
        mag_t m;
        double x;

        arf_init(a);
        arf_init(b);
        arf_init(c);
        mag_init(m);

        x = d_randtest2(state);
        x = ldexp(x, 1100 - n_randint(state, 2200));

        if (n_randint(state, 100) == 0) x = 0.0;
        if (n_randint(state, 100) == 0) x = D_INF;
        if (n_randint(state, 100) == 0) x = D_NAN;
        if (n_randint(state, 2)) x = -x;

        arf_set_d(a, x);
        mag_set_d(m, x);

        arf_set_mag(b, m);

        arf_set(c, a);
        arf_mul_ui(c, c, 1025, MAG_BITS, ARF_RND_UP);
        arf_mul_2exp_si(c, c, -10);

        MAG_CHECK_BITS(m)

        if (!(arf_cmpabs(a, b) <= 0 && arf_cmpabs(b, c) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("c = "); arf_print(c); flint_printf("\n\n");
            flint_abort();
        }

        mag_set_d_lower(m, x);
        arf_set_mag(b, m);

        arf_set(c, a);
        arf_mul_ui(c, c, 1023, MAG_BITS, ARF_RND_DOWN);
        arf_mul_2exp_si(c, c, -10);

        MAG_CHECK_BITS(m)

        if (!(arf_cmpabs(c, b) <= 0 && arf_cmpabs(b, a) <= 0))
        {
            flint_printf("FAIL (lower)\n\n");
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("c = "); arf_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(a);
        arf_clear(b);
        arf_clear(c);
        mag_clear(m);
    }

    TEST_FUNCTION_END(state);
}
