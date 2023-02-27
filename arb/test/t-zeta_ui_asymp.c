/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("zeta_ui_asymp....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t r;
        ulong n;
        mpfr_t s;
        slong prec;

        prec = 2 + n_randint(state, 1 << n_randint(state, 10));

        arb_init(r);
        mpfr_init2(s, prec + 100);

        n = 2 + n_randint(state, 1 << n_randint(state, 10));

        arb_zeta_ui_asymp(r, n, prec);
        mpfr_zeta_ui(s, n, MPFR_RNDN);

        if (!arb_contains_mpfr(r, s))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
            flint_printf("s = "); mpfr_printf("%.275Rf\n", s); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(r);
        mpfr_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
