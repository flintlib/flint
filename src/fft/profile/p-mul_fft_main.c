/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fft.h"
#include "profiler.h"

int main(void)
{
    mp_ptr r, x, y;
    mp_size_t n;
    slong i, nt;
    double t, FLINT_SET_BUT_UNUSED(_);

    flint_printf("         n     mpn_mul    1 thread    2 threads    4 threads    8 threads\n");

    for (n = 1000; n <= 10000000; n *= 1.025)
    {
        r = flint_malloc(2 * n * sizeof(mp_limb_t));
        x = flint_malloc(n * sizeof(mp_limb_t));
        y = flint_malloc(n * sizeof(mp_limb_t));

        for (i = 0; i < n; i++)
        {
            x[i] = UWORD_MAX - i - 1;
            y[i] = UWORD_MAX - i - 2;
        }

        flint_printf("%10wd", n); fflush(stdout);

        TIMEIT_START
        mpn_mul(r, x, n, y, n);
        TIMEIT_STOP_VALUES(_, t)
        flint_printf("%12g", t); fflush(stdout);

        for (nt = 1; nt <= 8; nt *= 2)
        {
            flint_set_num_threads(nt);
            TIMEIT_START
            flint_mpn_mul_fft_main(r, x, n, y, n);
            TIMEIT_STOP_VALUES(_, t)
            flint_printf("%12g", t); fflush(stdout);
        }

        flint_printf("\n");

        flint_free(r);
        flint_free(x);
        flint_free(y);
    }

   return 0;
}
