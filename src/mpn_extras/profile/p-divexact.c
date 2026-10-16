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

#define TIME(expr, res) \
    do { timeit_t __t; slong __r; expr; \
         TIMEIT_REPEAT(__t, __r) { expr; } TIMEIT_END_REPEAT(__t, __r); \
         res = (double) __t->wall * 0.001 / __r; } while (0)

int main(int argc, char * argv[])
{
    flint_rand_t state;
    mp_size_t n, bn, ratio;
    double tg, th, tt, tc, tk, ti;
    flint_rand_init(state);

    flint_printf("Exact division a = q b, q: n limbs, b: bn limbs; Hensel qr with n quotient limbs\n");
    flint_printf("%8s %8s %10s %10s %10s %10s %10s %10s\n", "bn", "n", "gmp_dex", "hensel", "tdiv_q", "bdiv_cl", "bdiv_km", "binv");

    for (ratio = 1; ratio <= 16; ratio *= 4)
    {
        for (bn = 32; bn <= 32768; bn *= 2)
        {
            mp_ptr a, b, q, q2, r, y;
            mp_size_t an;
            n = ratio * bn;
            if (n + bn > 400000) continue;

            b = flint_malloc(bn * sizeof(mp_limb_t));
            q = flint_malloc((n + 1) * sizeof(mp_limb_t));
            q2 = flint_malloc((n + 1) * sizeof(mp_limb_t));
            a = flint_malloc((n + bn) * sizeof(mp_limb_t));
            r = flint_malloc(bn * sizeof(mp_limb_t));
            y = flint_malloc(n * sizeof(mp_limb_t));
            flint_mpn_rrandom(b, state, bn);
            flint_mpn_rrandom(q2, state, n);
            b[0] |= 1;
            b[bn - 1] |= 1;
            q2[n - 1] |= 1;
            if (n >= bn) flint_mpn_mul(a, q2, n, b, bn); else flint_mpn_mul(a, b, bn, q2, n);
            an = n + bn - (a[n + bn - 1] == 0);

#if FLINT_HAVE_NATIVE_mpn_divexact
            TIME(mpn_divexact(q, a, an, b, bn), tg);
#else
            tg = 0;
#endif
            TIME(_flint_mpn_divexact_hensel(q, a, an, b, bn), th);
            TIME(flint_mpn_tdiv_q(q, a, an, b, bn), tt);
            TIME(flint_mpn_bdiv_qr_classical(q, r, a, an, b, bn, n), tc);
            TIME(flint_mpn_bdiv_qr_karp_markstein(q, r, a, an, b, bn, n), tk);
            TIME(flint_mpn_binv(y, b, bn, n), ti);

            flint_printf("%8wd %8wd %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", bn, n, tg, th, tt, tc, tk, ti);
            flint_free(a); flint_free(b); flint_free(q); flint_free(q2); flint_free(r); flint_free(y);
        }
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
