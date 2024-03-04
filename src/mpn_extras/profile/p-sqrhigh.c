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
#include "templates.h"

/* TODO: Remove me when fully implemented */
#if FLINT_HAVE_NATIVE_MPN_MULHIGH_BASECASE
void mpfr_sqrhigh_n(mp_ptr, mp_srcptr, mp_size_t);

#define N_MIN 1
#define N_MAX 2100

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

static mp_size_t inc(mp_size_t n)
{
    return FLINT_MAX(n + 1, n * 1.01);
}

int main(void)
{
    mp_limb_t rf[N_MAX];
    mp_limb_t rg[2 * N_MAX];
    mp_limb_t rm[2 * N_MAX];
    mp_limb_t xp[N_MAX];
    mp_size_t n;

    flint_printf("%.*s        sqr / sqrhigh || mpfr / flint\n", P, "                              ");
    for (n = N_MIN; n <= N_MAX; n = inc(n))
    {
        double t1, t2, t3, FLINT_SET_BUT_UNUSED(__);
        flint_printf("n = %" STR(P) "wd:", n);

        mpn_random2(xp, n);

        TIMEIT_START
        flint_mpn_sqrhigh(rf, xp, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        flint_mpn_sqr(rg, xp, n);
        TIMEIT_STOP_VALUES(__, t2)

        TIMEIT_START
        mpfr_sqrhigh_n(rm, xp, n);
        TIMEIT_STOP_VALUES(__, t3)

        flint_printf("   %7.3fx      ||  %7.3f\n", t2 / t1, t3 / t1);
    }

    flint_cleanup_master();

    return 0;
}

#undef N_MIN
#undef N_MAX
#else
int main(void)
{
    return 0;
}
#endif
