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

#define mpn_mul_basecase __gmpn_mul_basecase
void mpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

int main(void)
{
    mp_limb_t res1[24];
    mp_limb_t res2[24];
    mp_limb_t ap[16];
    mp_limb_t bp[8];
    slong mx, nx;

    flint_printf("       ");
    for (nx = 1; nx <= 8; nx++)
        flint_printf("%6wd  ", nx);
    flint_printf("\n\n");

    for (mx = 1; mx <= 16; mx++)
    {
        if (FLINT_HAVE_MUL_FUNC(mx, 1))
        {
            flint_printf("m = %2wd:", mx);

            for (nx = 1; nx <= FLINT_MIN(mx, WORD(8)); nx++)
            {
                double t1, t2, __attribute__((unused)) __;

                if (FLINT_HAVE_MUL_FUNC(mx, nx))
                {
                    mpn_random2(ap, mx);
                    mpn_random2(bp, nx);

                    TIMEIT_START
                    mpn_mul_basecase(res1, ap, mx, bp, nx);
                    TIMEIT_STOP_VALUES(__, t1)

                    TIMEIT_START
                    flint_mpn_mul_basecase(res2, ap, bp, mx, nx);
                    TIMEIT_STOP_VALUES(__, t2)

                    flint_printf("%7.2fx", t1 / t2);

                    if (mpn_cmp(res1, res2, mx + nx) != 0)
                        flint_abort();
                }
            }

            flint_printf("\n\n");
        }
    }

    flint_cleanup_master();

    return 0;
}
