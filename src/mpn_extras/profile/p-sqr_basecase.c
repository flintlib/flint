/*
    Copyright (C) 2023 Fredrik Johansson
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "profiler.h"

#define mpn_sqr_basecase __gmpn_sqr_basecase
void mpn_sqr_basecase(mp_ptr, mp_srcptr, mp_size_t);

int main(void)
{
    mp_limb_t res1[28];
    mp_limb_t res2[28];
    mp_limb_t ap[14];
    slong mx;

    for (mx = 1; mx <= 14; mx++)
    {
        double t1, t2, t3, __attribute__((unused)) __;

        if (FLINT_HAVE_SQR_FUNC(mx))
        {
            flint_printf("m = %2wd:", mx);

            mpn_random2(ap, mx);

            TIMEIT_START
            mpn_sqr_basecase(res1, ap, mx);
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            flint_mpn_sqr_basecase(res2, ap, mx);
            TIMEIT_STOP_VALUES(__, t2)

            TIMEIT_START
            flint_mpn_mul_basecase(res2, ap, ap, mx, mx);
            TIMEIT_STOP_VALUES(__, t3)

            flint_printf("%7.2fx", t1 / t2);
            flint_printf("    %7.2fx", t1 / t3);

            /* if (mpn_cmp(res1, res2, 2 * mx) != 0) */
            /*     flint_abort(); */

            flint_printf("\n\n");
        }
    }

    flint_cleanup_master();

    return 0;
}
