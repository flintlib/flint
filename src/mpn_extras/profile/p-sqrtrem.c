/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "mpn_extras.h"
#include "profiler.h"
#include "ulong_extras.h"

#define TIME(expr, res) \
    do { timeit_t __t; slong __r; expr; \
         TIMEIT_REPEAT(__t, __r) { expr; } TIMEIT_END_REPEAT(__t, __r); \
         res = (double) __t->wall * 0.001 / __r; } while (0)

int main(int argc, char * argv[])
{
    flint_rand_t state;
    mp_size_t an, sn;
    double tg, tn, tg2, tn2, tb;
    flint_rand_init(state);

    flint_printf("Small square roots: an limbs, times in seconds (full-size random input)\n");
    flint_printf("%8s %10s %10s %10s %10s\n", "an", "gmp_qr", "flint_qr", "gmp_q", "flint_q");

    for (an = 1; an <= 6; an++)
    {
        mp_limb_t a[6 * 64], s[4 * 64], r[7 * 64];
        slong i;

        /* 64 inputs per timing, to average over the data dependent paths */
        flint_mpn_rrandom(a, state, an * 64);
        for (i = 0; i < 64; i++)
            a[i * an + an - 1] |= UWORD(1) << (FLINT_BITS - 1 - n_randint(state, 4));

        sn = (an + 1) / 2;
        TIME(for (i = 0; i < 64; i++) mpn_sqrtrem(s + i * sn, r + i * (an + 1), a + i * an, an), tg);
        TIME(for (i = 0; i < 64; i++) flint_mpn_sqrtrem(s + i * sn, r + i * (an + 1), a + i * an, an), tn);
        TIME(for (i = 0; i < 64; i++) mpn_sqrtrem(s + i * sn, NULL, a + i * an, an), tg2);
        TIME(for (i = 0; i < 64; i++) flint_mpn_sqrtrem(s + i * sn, NULL, a + i * an, an), tn2);

        flint_printf("%8wd %10.3e %10.3e %10.3e %10.3e\n", an, tg / 64, tn / 64, tg2 / 64, tn2 / 64);
    }

    flint_printf("\nSquare root: an limbs, times in seconds\n");
    flint_printf("%8s %10s %10s %10s %10s %10s\n", "an", "gmp_qr", "newton_qr", "gmp_q", "newton_q", "bsqrt");

    for (an = 32; an <= 200000; an = an * 3 / 2)
    {
        mp_ptr a, s, r, t;
        sn = (an + 1) / 2;
        a = flint_malloc(an * sizeof(mp_limb_t));
        s = flint_malloc((sn + 1) * sizeof(mp_limb_t));
        r = flint_malloc(an * sizeof(mp_limb_t));
        t = flint_malloc(an * sizeof(mp_limb_t));
        flint_mpn_rrandom(a, state, an);
        a[an - 1] |= 1;
        a[0] = (a[0] & ~UWORD(7)) | 1;

        TIME(mpn_sqrtrem(s, r, a, an), tg);
        TIME(_flint_mpn_sqrtrem_newton(s, r, a, an), tn);
        TIME(mpn_sqrtrem(s, NULL, a, an), tg2);
        TIME(_flint_mpn_sqrtrem_newton(s, NULL, a, an), tn2);
        TIME(flint_mpn_bsqrt(t, a, an, an), tb);

        flint_printf("%8wd %10.3e %10.3e %10.3e %10.3e %10.3e\n", an, tg, tn, tg2, tn2, tb);
        flint_free(a); flint_free(s); flint_free(r); flint_free(t);
    }

    flint_rand_clear(state);
    return 0;
}
