/*
    Copyright (C) 2023 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_HAVE_FFT_SMALL

#include "fft_small.h"
#include "profiler.h"

#define N_MIN_MUL 1500
#define N_MAX_MUL 1600

#define N_MIN_SQR 2900
#define N_MAX_SQR 3200

int main(void)
{
    mp_ptr x, y, r, s;
    slong n;

    x = flint_malloc(sizeof(mp_limb_t) * FLINT_MAX(N_MAX_MUL, N_MAX_SQR));
    y = flint_malloc(sizeof(mp_limb_t) * FLINT_MAX(N_MAX_MUL, N_MAX_SQR));
    r = flint_malloc(2 * sizeof(mp_limb_t) * FLINT_MAX(N_MAX_MUL, N_MAX_SQR));
    s = flint_malloc(2 * sizeof(mp_limb_t) * FLINT_MAX(N_MAX_MUL, N_MAX_SQR));

    flint_mpn_rrandom(x, state, FLINT_MAX(N_MAX_MUL, N_MAX_SQR));
    flint_mpn_rrandom(y, state, FLINT_MAX(N_MAX_MUL, N_MAX_SQR));

    flint_printf("mpn_mul_n vs fft_small\n\n");

    for (n = N_MIN_MUL; n <= N_MAX_MUL; n += 5)
    {
        double t1, t2, FLINT_SET_BUT_UNUSED(__);

        flint_printf("n = %4wd: ", n);

        TIMEIT_START
        mpn_mul_n(r, x, y, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        mpn_mul_default_mpn_ctx(s, x, n, y, n);
        TIMEIT_STOP_VALUES(__, t2)

        flint_printf("%.2f\n", t1 / t2);
    }

    flint_printf("mpn_sqr vs fft_small\n\n");

    for (n = N_MIN_SQR; n <= N_MAX_SQR; n += 10)
    {
        double t1, t2, FLINT_SET_BUT_UNUSED(__);

        flint_printf("n = %4wd: ", n);

        TIMEIT_START
        mpn_sqr(r, x, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        mpn_mul_default_mpn_ctx(s, x, n, x, n);
        TIMEIT_STOP_VALUES(__, t2)

        flint_printf("%.2f\n", t1 / t2);
    }

    flint_free(x);
    flint_free(y);
    flint_free(r);
    flint_free(s);

    flint_cleanup_master();

    return 0;
}
#else
int main(void) { return 0; }
#endif
