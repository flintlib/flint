/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("sdiv_qrnnd....");
    fflush(stdout);

    for (i = 0; i < 1000000*flint_test_multiplier(); i++)
    {
        mp_limb_signed_t d, nh, nl, q, r, ph, pl;

        do
        {
            d = n_randtest_not_zero(state);
            nh = n_randtest(state);
        } while ((FLINT_ABS(nh) >= FLINT_ABS(d)/2) || (nh == WORD_MIN));

        nl = n_randtest(state);

        sdiv_qrnnd(q, r, nh, nl, d);
        smul_ppmm(ph, pl, d, q);
        add_ssaaaa(ph, pl, ph, pl, FLINT_SIGN_EXT(r), r);

        if (ph != nh ||
            pl != nl ||
            (d > 0 && r >= d) ||
            (d < 0 && r <= d))
        {
            flint_printf("FAIL:\n");
            flint_printf("nh = %wd, nl = %wd\n", nh, nl);
            flint_printf("d = %wd\n", d);
            flint_printf("q = %wd\n", q);
            flint_printf("r = %wd\n", r);
            flint_printf("ph = %wu, pl = %wu\n", ph, pl);
            flint_abort();
        }
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
