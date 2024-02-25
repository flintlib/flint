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

#if FLINT_HAVE_NATIVE_2ADD_N_INPLACE
void mpfr_mulhigh_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

#define N_MIN 4
#define N_MAX 950

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
    if (n < 16)
        return n + 1;
    else if (n < 32)
        return n + 3;
    else if (n < 64)
        return n + 9;
    else if (n < 128)
        return n + 17;
    else if (n < 256)
        return n + 27;
    else
        return n + 57;
}

int main(void)
{
    mp_limb_t rf[N_MAX];
    mp_limb_t rg[N_MAX];
    mp_limb_t xp[N_MAX];
    mp_limb_t yp[N_MAX];
    mp_size_t n;

    flint_printf("2 mpn_add_n / flint_mpn_2add_n_inplace\n");
    for (n = N_MIN; n <= N_MAX; n = inc(n))
    {
        double t1, t2, FLINT_SET_BUT_UNUSED(__);

        flint_printf("n = %" STR(P) "wd:", n);

        mpn_random2(rf, n);
        mpn_random2(xp, n);
        mpn_random2(yp, n);
        flint_mpn_copyi(rg, rf, n);

        TIMEIT_START
        flint_mpn_2add_n_inplace(rf, xp, yp, n);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        mpn_add_n(rg, rg, xp, n);
        mpn_add_n(rg, rg, yp, n);
        TIMEIT_STOP_VALUES(__, t2)

        flint_printf("%7.2fx\n", t2 / t1);
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
