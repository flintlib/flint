/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "profiler.h"

#if !FLINT_HAVE_ADX
# error
#endif

#define mpn_mul_basecase __gmpn_mul_basecase
void mpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
void mpfr_mulhigh_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

int main(void)
{
#define N_MIN 1
#define N_MAX 64

    mp_limb_t rf[N_MAX];
    mp_limb_t rg[2 * N_MAX];
    mp_limb_t rm[2 * N_MAX];
    mp_limb_t xp[N_MAX];
    mp_limb_t yp[N_MAX];
    mp_size_t n;

    flint_printf("Best of mpn_mul_n, mpn_mul_basecase and mpfr_mulhigh_n\nagainst flint_mpn_mulhigh_n_basecase:\n");
    for (n = N_MIN; n <= N_MAX; n++)
    {
        double t1, t2, t3, t4, tmin, __attribute__((unused)) __;
        int type;
        flint_printf("n = %2wd:", n);

        mpn_random2(xp, n);
        mpn_random2(yp, n);

        TIMEIT_START
        flint_mpn_mulhigh_n(rf, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        mpn_mul_basecase(rg, xp, n, yp, n);
        TIMEIT_STOP_VALUES(__, t2)

        TIMEIT_START
        mpn_mul_n(rg, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t3)

        TIMEIT_START
        mpfr_mulhigh_n(rm, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t4)

        if (t2 < t3)
        {
            tmin = t2;
            type = 2;
        }
        else
        {
            tmin = t3;
            type = 3;
        }

        if (tmin < t4)
            ;
        else
        {
            tmin = t4;
            type = 4;
        }

        flint_printf("%7.2fx  (%s)\n",
                tmin / t1,
                (type == 2) ? "mpn_mul_basecase" : (type == 3) ? "mpn_mul_n" : "mpfr_mulhigh_n");
    }

    flint_cleanup_master();

    return 0;
}

#undef N_MIN
#undef N_MAX
