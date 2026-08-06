/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "profiler.h"
#include "radix.h"

#include "fmpz.h"
#include "arb.h"

/* s = floor(sqrt(a)) and perfect-square detection via mpn, mirroring
   what int radix_sqrt computes */
int
flint_mpn_sqrt(nn_ptr s, nn_srcptr a, slong an, nn_ptr scratch)
{
    return mpn_sqrtrem(s, scratch, a, an) == 0;
}

int main()
{
    radix_t radix;
    slong digits, nd, n1, n2;
    nn_ptr a, c, d;
    int FLINT_SET_BUT_UNUSED(exact);
    double tmpn, tradix, FLINT_SET_BUT_UNUSED(tt);

    flint_rand_t state;
    flint_rand_init(state);

    flint_printf("   decimal      mpn  decimal        time        time    relative\n");
    flint_printf("    digits    limbs    limbs    mpn_sqrtrem radix_sqrt  time\n\n");

    for (nd = 1; nd <= 100000000; nd = FLINT_MAX(nd + 1, nd * 1.5))
    {
        radix_init(radix, 10, 0);

        digits = nd * 19;

        /* Number of full-word limbs */
        n1 = (slong) (digits * (log(10) / (FLINT_BITS * log(2))) + 1.0);
        /* Number of radix 10^e limbs */
        n2 = (digits + radix->exp - 1) / radix->exp;

        a = flint_malloc(2 * n2 * sizeof(ulong));
        c = flint_malloc(n2 * sizeof(ulong));
        d = flint_malloc(2 * n2 * sizeof(ulong));   /* mpn_sqrtrem needs n limbs of remainder space */

        flint_mpn_urandomb(a, state, 2 * n1 * FLINT_BITS);
        a[2 * n1 - 1] |= (UWORD(1) << (FLINT_BITS - 1));

        exact = flint_mpn_sqrt(c, a, 2 * n1, d);
        TIMEIT_START;
        exact = flint_mpn_sqrt(c, a, 2 * n1, d);
        TIMEIT_STOP_VALUES(tt, tmpn);

        radix_rand_limbs(a, state, 2 * n2, radix);
        if (a[2 * n2 - 1] == 0)
            a[2 * n2 - 1] = 1;

        exact = radix_sqrt(c, a, 2 * n2, radix);
        TIMEIT_START;
        exact = radix_sqrt(c, a, 2 * n2, radix);
        TIMEIT_STOP_VALUES(tt, tradix);

        flint_printf("%10wd %8wd %8wd    %8g    %8g    %.3fx\n",
            digits, n1, n2, tmpn, tradix, tradix / tmpn);

        flint_free(a);
        flint_free(c);
        flint_free(d);

        radix_clear(radix);
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
