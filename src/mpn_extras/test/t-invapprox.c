/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_invapprox, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        mp_ptr x, q, q2, num, r, t;
        mp_size_t xn, n, qn, i;
        int ok;

        switch (n_randint(state, 5))
        {
            case 0:
                xn = 1 + n_randint(state, 4);
                n = xn + n_randint(state, 20);
                break;
            case 1:
                xn = 1 + n_randint(state, 10);
                n = xn + n_randint(state, 3000);
                break;
            case 2:
                xn = 1 + n_randint(state, 400);
                n = xn + n_randint(state, xn + 5);
                break;
            case 3:
                /* around and above the Newton cutoffs */
                xn = FLINT_MPN_INV_NEWTON_LONG_CUTOFF + n_randint(state, 400);
                n = xn + FLINT_MPN_INV_NEWTON_LONG_CUTOFF + n_randint(state, 1500);
                break;
            default:
                xn = 1 + n_randint(state, 1200);
                n = xn + n_randint(state, 1200);
                break;
        }
        qn = n - xn + 2;

        x = flint_malloc(xn * sizeof(mp_limb_t));
        q = flint_malloc(qn * sizeof(mp_limb_t));
        q2 = flint_malloc(qn * sizeof(mp_limb_t));
        num = flint_malloc((n + 1) * sizeof(mp_limb_t));
        r = flint_malloc(xn * sizeof(mp_limb_t));
        t = flint_malloc(qn * sizeof(mp_limb_t));

        flint_mpn_rrandom(x, state, xn);
        switch (n_randint(state, 8))
        {
            case 0:
                /* a power of B, or of 2 */
                flint_mpn_zero(x, xn);
                x[xn - 1] = UWORD(1) << n_randint(state, FLINT_BITS);
                break;
            case 1:
                /* B^xn - 1 */
                for (i = 0; i < xn; i++)
                    x[i] = UWORD_MAX;
                break;
            case 2:
                /* power of 2 + small */
                flint_mpn_zero(x, xn);
                x[xn - 1] = UWORD(1) << n_randint(state, FLINT_BITS);
                x[0] += n_randint(state, 3);
                break;
            case 3:
                x[xn - 1] |= UWORD(1) << (FLINT_BITS - 1);
                break;
            default:
                break;
        }
        if (x[xn - 1] == 0)
            x[xn - 1] = 1;

        /* reference: floor(B^n / x) */
        flint_mpn_zero(num, n);
        num[n] = 1;
        q2[qn - 1] = 0;
        mpn_tdiv_qr(q2, r, 0, num, n + 1, x, xn);

        flint_mpn_invapprox(q, x, xn, n);

        /* q = q2 or q2 + 1 */
        ok = (mpn_cmp(q, q2, qn) == 0);
        if (!ok)
        {
            flint_mpn_copyi(t, q2, qn);
            ok = (mpn_add_1(t, t, qn, 1) == 0) && (mpn_cmp(q, t, qn) == 0);
        }

        if (!ok)
            TEST_FUNCTION_FAIL("xn = %wd, n = %wd\nx = %{ulong*}\nq = %{ulong*}\nq2 = %{ulong*}\n",
                xn, n, x, xn, q, qn, q2, qn);

        flint_free(x);
        flint_free(q);
        flint_free(q2);
        flint_free(num);
        flint_free(r);
        flint_free(t);
    }

    TEST_FUNCTION_END(state);
}
