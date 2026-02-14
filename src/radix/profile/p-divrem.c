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

/* there is no flint_mpn_divrem yet, so mock one up */
void
flint_mpn_divrem(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn)
{
    if (bn < 1000)
    {
        mpn_tdiv_qr(q, r, 0, a, an, b, bn);
    }
    else
    {
        fmpz_t fa, fb, fq, fr;

        fmpz_init(fa);
        fmpz_init(fb);
        fmpz_init(fq);
        fmpz_init(fr);

        fmpz_set_ui_array(fa, a, an);
        fmpz_set_ui_array(fb, b, bn);

        fmpz_tdiv_qr(fq, fr, fa, fb);

        //fmpz_get_ui_array(q, an - bn + 1, fq);
        //fmpz_get_ui_array(r, bn, fr);

        fmpz_clear(fa);
        fmpz_clear(fb);
        fmpz_clear(fq);
        fmpz_clear(fr);
    }
}

int main()
{
    radix_t radix;
    slong digits, nd, n1, n2;
    nn_ptr a, b, c, d;
    double tmpn, tradix, FLINT_SET_BUT_UNUSED(tt);

    flint_rand_t state;
    flint_rand_init(state);

    flint_printf("   decimal      mpn  decimal        time        time    relative\n");
    flint_printf("    digits    limbs    limbs flint_mpn_divrem radix_divrem  time\n\n");

    for (nd = 1; nd <= 100000000; nd = FLINT_MAX(nd + 1, nd * 1.5))
    {
        radix_init(radix, 10, 0);

        digits = nd * 19;

        /* Number of full-word limbs */
        n1 = (slong) (digits * (log(10) / (FLINT_BITS * log(2))) + 1.0);
        /* Number of radix 10^e limbs */
        n2 = (digits + radix->exp - 1) / radix->exp;

        a = flint_malloc(2 * n2 * sizeof(ulong));
        b = flint_malloc(n2 * sizeof(ulong));
        c = flint_malloc((n2 + 1) * sizeof(ulong));
        d = flint_malloc(n2 * sizeof(ulong));

        flint_mpn_urandomb(a, state, 2 * n1 * FLINT_BITS);
        flint_mpn_urandomb(b, state, n1 * FLINT_BITS);

        flint_mpn_divrem(c, d, a, 2 * n1, b, n1);
        TIMEIT_START;
        flint_mpn_divrem(c, d, a, 2 * n1, b, n1);
        TIMEIT_STOP_VALUES(tt, tmpn);

        radix_rand_limbs(a, state, 2 * n2, radix);
        radix_rand_limbs(b, state, n2, radix);

        radix_divrem(c, d, a, 2 * n2, b, n2, radix);
        TIMEIT_START;
        radix_divrem(c, d, a, 2 * n2, b, n2, radix);
        TIMEIT_STOP_VALUES(tt, tradix);

        flint_printf("%10wd %8wd %8wd    %8g    %8g    %.3fx\n",
            digits, n1, n2, tmpn, tradix, tradix / tmpn);

        flint_free(a);
        flint_free(b);
        flint_free(c);
        flint_free(d);

        radix_clear(radix);
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}

