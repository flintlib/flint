/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "padic_radix.h"
#include "../impl.h"

#define INF (2 * PADIC_RADIX_PACK_BIG)

/* sequences of various structures, with +infinity entries */
static void
_randseq(slong * y, slong len, flint_rand_t state)
{
    slong i, kind = n_randint(state, 9), R = 1 + n_randint(state, n_randint(state, 2) ? 50 : 100000);
    slong m = 1 + n_randint(state, 10), slope = (slong) n_randint(state, 21) - 10, lev = 0;
    slong pinf = n_randint(state, 3) ? 0 : 1 + n_randint(state, 4);

    for (i = 0; i < len; i++)
    {
        switch (kind)
        {
            case 0: y[i] = n_randint(state, R); break;
            case 1: y[i] = slope * i + n_randint(state, 3); break;
            case 2: y[i] = (i % m == 0) ? 0 : R; break;
            case 3: if (i % m == 0) lev = n_randint(state, R); y[i] = lev + n_randint(state, 3); break;
            case 4: y[i] = (i - len / 2) * (i - len / 2) / (1 + m); break;
            case 5: y[i] = n_randint(state, 50) ? n_randint(state, 4) : n_randint(state, R); break;
            case 6: y[i] = slope * i + ((i % m == 0) ? R : 0); break;
            case 7: y[i] = 5; break;
            default: y[i] = (slong) n_randint(state, 2 * R + 1) - R; break;
        }
        if (pinf && n_randint(state, pinf + 1) == 0)
            y[i] = INF;
        if (kind == 2 && pinf && i % m == m - 1)
            y[i] = INF;
    }
}

TEST_FUNCTION_START(padic_radix_minplus, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        slong alen, blen, klo, khi, k;
        slong * a, * b, * r1, * r2;
        int exact;

        alen = 1 + n_randint(state, n_randint(state, 4) ? 50 : 1500);
        blen = 1 + n_randint(state, n_randint(state, 4) ? 50 : 1500);
        khi = 1 + n_randint(state, alen + blen - 1);
        klo = n_randint(state, khi);

        a = flint_malloc(sizeof(slong) * (alen + blen + 2 * (khi - klo)));
        b = a + alen;
        r1 = b + blen;
        r2 = r1 + (khi - klo);

        _randseq(a, alen, state);
        _randseq(b, blen, state);

        _padic_radix_minplus_conv_exact(r1, a, alen, b, blen, klo, khi);

        for (exact = 0; exact <= 1; exact++)
        {
            int claimed = _padic_radix_minplus_conv(r2, a, alen, b, blen, klo, khi, exact);

            for (k = 0; k < khi - klo; k++)
            {
                slong x = (r1[k] >= PADIC_RADIX_PACK_NONE) ? INF : r1[k];
                slong y = (r2[k] >= PADIC_RADIX_PACK_NONE) ? INF : r2[k];

                if (((exact || claimed) && x != y) || y > x || (x == INF && y != INF))
                {
                    flint_printf("FAIL: exact = %d, claimed = %d, alen = %wd, blen = %wd, k = %wd: %wd, %wd\n",
                        exact, claimed, alen, blen, k + klo, x, y);
                    flint_abort();
                }
            }
        }

        flint_free(a);
    }

    TEST_FUNCTION_END(state);
}
