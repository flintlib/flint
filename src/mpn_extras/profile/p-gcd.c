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

/* Usage: p-gcd [max_limbs]. Compares _flint_mpn_gcd_inplace and
   _flint_mpn_gcdext_inplace with GMP's mpn_gcd and mpn_gcdext on random
   operands of n limbs. */

#define TIME(expr, res) \
    do { timeit_t __t; slong __r; expr; \
         TIMEIT_REPEAT(__t, __r) { expr; } TIMEIT_END_REPEAT(__t, __r); \
         res = (double) __t->wall * 0.001 / __r; } while (0)

int main(int argc, char * argv[])
{
    flint_rand_t state;
    slong n, maxn = (argc > 1) ? atol(argv[1]) : 100000;
    double t1, t2, t3, t4;

    flint_rand_init(state);

    flint_printf("gcd and gcdext of n-limb operands, times in seconds\n");
    flint_printf("%8s %10s %10s %7s %10s %10s %7s\n", "n", "mpn_gcd", "flint", "ratio",
        "mpn_gcdext", "flint", "ratio");

    for (n = 2; n <= maxn; n = (n < 16) ? n * 2 : (slong) (n * 1.6))
    {
        mp_ptr a, b, a2, b2, g, s;
        mp_size_t sn;

        a = flint_malloc((n + 1) * sizeof(mp_limb_t));
        b = flint_malloc((n + 1) * sizeof(mp_limb_t));
        a2 = flint_malloc((n + 1) * sizeof(mp_limb_t));
        b2 = flint_malloc((n + 1) * sizeof(mp_limb_t));
        g = flint_malloc((n + 1) * sizeof(mp_limb_t));
        s = flint_malloc((n + 1) * sizeof(mp_limb_t));

        flint_mpn_urandomb(a, state, n * FLINT_BITS);
        flint_mpn_urandomb(b, state, n * FLINT_BITS);
        a[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
        b[n - 1] |= UWORD(1) << (FLINT_BITS - 2);
        a[0] |= 1;

#define COPY (flint_mpn_copyi(a2, a, n), flint_mpn_copyi(b2, b, n))
        TIME((COPY, mpn_gcd(g, a2, n, b2, n)), t1);
        TIME((COPY, _flint_mpn_gcd_inplace(g, a2, n, b2, n)), t2);
        TIME((COPY, mpn_gcdext(g, s, &sn, a2, n, b2, n)), t3);
        TIME((COPY, _flint_mpn_gcdext_inplace(g, s, &sn, a2, n, b2, n)), t4);
#undef COPY

        flint_printf("%8wd %10.3e %10.3e %7.3f %10.3e %10.3e %7.3f\n",
            n, t1, t2, t1 / t2, t3, t4, t3 / t4);

        flint_free(a);
        flint_free(b);
        flint_free(a2);
        flint_free(b2);
        flint_free(g);
        flint_free(s);
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
