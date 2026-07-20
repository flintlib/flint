/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "longlong.h"
#include "mpn_extras.h"

/*
    flint_mpn_mulmid_n(rp, a, b, n) is the balanced exact middle product of
    {a, 2n-1} and {b, n}, writing n + 2 limbs.  It is the sum of the diagonals
    n-1 <= p+q < 2n-1 shifted down by n-1; the high n limbs (rp[2..n+2)) are
    exact and the low two are guard limbs.  Reference: accumulate exactly that
    band and compare the high n limbs.
*/
TEST_FUNCTION_START(flint_mpn_mulmid_n, state)
{
    slong ix;

#if !FLINT_HAVE_NATIVE_mpn_mulmid_n
    TEST_FUNCTION_END_SKIPPED(state);
#else
    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        mp_size_t n, an, p, q;
        mp_ptr a, b, rp, acc;

        n = 2 + n_randint(state, 40);
        if (n_randint(state, 500) == 0)
            n = 2 + n_randint(state, 400);
        an = 2 * n - 1;

        a = flint_malloc(sizeof(mp_limb_t) * an);
        b = flint_malloc(sizeof(mp_limb_t) * n);
        rp = flint_malloc(sizeof(mp_limb_t) * (n + 2));
        acc = flint_calloc(n + 3, sizeof(mp_limb_t));

        flint_mpn_rrandom(a, state, an);
        flint_mpn_rrandom(b, state, n);
        flint_mpn_rrandom(rp, state, n + 2);        /* poison */

        flint_mpn_mulmid_n(rp, a, b, n);

        for (p = 0; p < an; p++)
            for (q = 0; q < n; q++)
            {
                mp_size_t s = p + q, off;
                mp_limb_t hi, lo;

                if (s < n - 1 || s >= 2 * n - 1)
                    continue;

                umul_ppmm(hi, lo, a[p], b[q]);
                off = s - (n - 1);
                mpn_add_1(acc + off, acc + off, (n + 3) - off, lo);
                mpn_add_1(acc + off + 1, acc + off + 1, (n + 3) - off - 1, hi);
            }

        if (mpn_cmp(rp + 2, acc + 2, n) != 0)
            TEST_FUNCTION_FAIL(
                    "ix = %wd, n = %wd\n"
                    "a = %{ulong*}\nb = %{ulong*}\n"
                    "exact high limbs = %{ulong*}\n"
                    "got              = %{ulong*}\n",
                    ix, n, a, an, b, n, acc + 2, n, rp + 2, n);

        flint_free(a);
        flint_free(b);
        flint_free(rp);
        flint_free(acc);
    }

    TEST_FUNCTION_END(state);
#endif
}
