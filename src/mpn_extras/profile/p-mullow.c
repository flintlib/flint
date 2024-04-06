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

#define mpn_mullo_basecase __gmpn_mullo_basecase
void mpn_mullo_basecase(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

#define N_MIN 1
#define N_MAX 30

#if N_MAX < 9
# define P 1
#elif N_MAX < 99
# define P 2
#elif N_MAX < 999
# define P 3
#elif N_MAX < 9999
# define P 4
#else
# define P 5
#endif

#define _STR(x) #x
#define STR(x) _STR(x)

int
main(void)
{
    mp_limb_t rp[2 * N_MAX];
    mp_limb_t xp[N_MAX];
    mp_limb_t yp[N_MAX];
    mp_size_t n;

    flint_printf("%.*s      mullo_basecase / FLINT || mul / mullow\n", P, "                              ");
    for (n = N_MIN; n <= N_MAX; n++)
    {
        double t1, t2, t3, FLINT_SET_BUT_UNUSED(__);
        flint_printf("n = %" STR(P) "wd:", n);

        mpn_random2(xp, n);
        mpn_random2(yp, n);

        TIMEIT_START
        mpn_mullo_basecase(rp, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        flint_mpn_mul_n(rp, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t2)

        TIMEIT_START
        flint_mpn_mullow_n(rp, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t3)

        flint_printf("         %7.2fx       ||  %7.2f\n", t1 / t3, t2 / t3);
    }

    flint_cleanup_master();

    return 0;
}

#undef N_MIN
#undef N_MAX
