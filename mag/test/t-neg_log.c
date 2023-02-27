/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("neg_log....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y, z, z2;
        fmpz_t n;
        mag_t xb, yb;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(z2);

        mag_init(xb);
        mag_init(yb);

        fmpz_init(n);

        mag_randtest_special(xb, state, 25);
        mag_randtest_special(yb, state, 25);

        if (n_randint(state, 2))
        {
            mag_neg_log(yb, xb);

            mag_get_fmpr(x, xb);
            if (mag_cmp_2exp_si(xb, 0) >= 0)
                fmpr_zero(z);
            else
            {
                fmpr_log(z, x, MAG_BITS, FMPR_RND_UP);
                fmpr_neg(z, z);
            }
        }
        else
        {
            mag_get_fmpr(x, xb);

            fmpz_randtest(n, state, 100);
            mag_mul_2exp_fmpz(xb, xb, n);
            mag_neg_log(yb, xb);

            if (mag_cmp_2exp_si(xb, 0) >= 0)
            {
                fmpr_zero(z);
            }
            else
            {
                fmpr_log(z, x, 2 * MAG_BITS, FMPR_RND_UP);
                fmpr_neg(z, z);
                fmpr_set_ui(z2, 2);
                fmpr_log(z2, z2, 2 * MAG_BITS, FMPR_RND_UP);
                fmpr_submul_fmpz(z, z2, n, 2 * MAG_BITS, FMPR_RND_UP);
            }
        }

        mag_get_fmpr(y, yb);

        fmpr_mul_ui(z2, z, 1025, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_2exp_si(z2, z2, -10);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(z, y) <= 0 && fmpr_cmpabs(y, z2) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); fmpr_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        mag_neg_log(xb, xb);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(z2);

        mag_clear(xb);
        mag_clear(yb);

        fmpz_clear(n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

