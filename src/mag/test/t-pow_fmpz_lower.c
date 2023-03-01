/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"
#include "mag.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("pow_fmpz_lower....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y, w;
        mag_t xb, yb;
        fmpz_t e;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(w);
        fmpz_init(e);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 80);
        mag_randtest_special(yb, state, 80);
        fmpz_randtest(e, state, 200);

        mag_get_fmpr(x, xb);

        fmpr_pow_sloppy_fmpz(y, x, e, 2 * MAG_BITS, FMPR_RND_DOWN);
        if (fmpr_is_nan(y))
            fmpr_pos_inf(y);

        mag_pow_fmpz_lower(yb, xb, e);
        mag_get_fmpr(w, yb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(w, y) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("w = "); fmpr_print(w); flint_printf("\n\n");
            flint_abort();
        }

        mag_pow_fmpz_lower(xb, xb, e);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            mag_print(xb); flint_printf("\n\n");
            mag_print(yb); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(w);
        fmpz_clear(e);

        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

