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

    flint_printf("set_d_2exp_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpr_t a, b, c;
        fmpz_t e;
        mag_t m;
        double x;

        fmpr_init(a);
        fmpr_init(b);
        fmpr_init(c);
        fmpz_init(e);
        mag_init(m);

        x = d_randtest2(state);
        x = ldexp(x, 1100 - n_randint(state, 2200));

        if (n_randint(state, 100) == 0) x = 0.0;
        if (n_randint(state, 100) == 0) x = D_INF;
        if (n_randint(state, 100) == 0) x = D_NAN;
        if (n_randint(state, 2)) x = -x;

        fmpz_randtest(e, state, 100);

        fmpr_set_d(a, x);
        fmpr_abs(a, a);
        fmpr_mul_2exp_fmpz(a, a, e);

        mag_set_d_2exp_fmpz(m, x, e);
        mag_get_fmpr(b, m);

        fmpr_set(c, a);
        fmpr_mul_ui(c, c, 1025, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_2exp_si(c, c, -10);

        MAG_CHECK_BITS(m)

        if (!(fmpr_cmpabs(a, b) <= 0 && fmpr_cmpabs(b, c) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = %g\n\n", x);
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            flint_printf("a = "); fmpr_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpr_print(b); flint_printf("\n\n");
            flint_printf("c = "); fmpr_print(c); flint_printf("\n\n");
            flint_abort();
        }

        mag_set_d_2exp_fmpz_lower(m, x, e);
        mag_get_fmpr(b, m);

        fmpr_set(c, a);
        fmpr_mul_ui(c, c, 1023, MAG_BITS, FMPR_RND_DOWN);
        fmpr_mul_2exp_si(c, c, -10);

        MAG_CHECK_BITS(m)

        if (!(fmpr_cmpabs(c, b) <= 0 && fmpr_cmpabs(b, a) <= 0))
        {
            flint_printf("FAIL (lower)\n\n");
            flint_printf("a = "); fmpr_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpr_print(b); flint_printf("\n\n");
            flint_printf("c = "); fmpr_print(c); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(a);
        fmpr_clear(b);
        fmpr_clear(c);
        fmpz_clear(e);
        mag_clear(m);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

