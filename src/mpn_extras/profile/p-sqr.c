/*
    Copyright (C) 2023 Fredrik Johansson
    Copyright (C) 2023, 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "mpn_extras.h"
#include "profiler.h"

#define N_MAX FLINT_MPN_SQR_FUNC_TAB_WIDTH

#if N_MAX > 0

#define mpn_sqr_basecase __gmpn_sqr_basecase
void mpn_sqr_basecase(mp_ptr, mp_srcptr, mp_size_t);

static void measure(mp_ptr rp1, mp_ptr rp2, mp_ptr ap, slong mx)
{
    double t1, t2, FLINT_SET_BUT_UNUSED(__);

    if (FLINT_HAVE_SQR_FUNC(mx))
    {
        flint_printf("m = %2wd:", mx);

        mpn_random2(ap, mx);

        TIMEIT_START
        mpn_sqr_basecase(rp1, ap, mx);
        TIMEIT_STOP_VALUES(__, t1)

        TIMEIT_START
        flint_mpn_sqr(rp2, ap, mx);
        TIMEIT_STOP_VALUES(__, t2)

        flint_printf("%7.2fx\n", t1 / t2);
    }
}

int main(int argc, char ** argv)
{
    mp_limb_t res1[2 * N_MAX];
    mp_limb_t res2[2 * N_MAX];
    mp_limb_t ap[N_MAX];
    slong mx;

    if (argc == 2)
    {
        measure(res1, res2, ap, strtol(argv[1], NULL, 10));
        goto end;
    }

    flint_printf("mpn_sqr_basecase / flint_mpn_sqr\n");
    for (mx = 1; mx <= N_MAX; mx++)
        measure(res1, res2, ap, mx);

end:
    flint_cleanup_master();

    return 0;
}
#else
int main(void) { return 0; }
#endif
