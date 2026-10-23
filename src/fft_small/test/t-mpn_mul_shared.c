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
#include "fft_small.h"

/* z1 = a x and z2 = b x against flint_mpn_mul in both geometries
   (two-prime up to FLINT_MPN_MUL_NP2_MAX_AN limbs, wide beyond), with
   the outputs aliasing the inputs as allowed */
TEST_FUNCTION_START(_fft_small_mpn_mul_shared, state)
{
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        ulong an, bn, xn, size;
        nn_ptr a, b, x, z1, z2, r1, r2, buf;
        int alias;

        size = n_randint(state, 12);
        if (size == 0)
            size = 150000;
        else if (size <= 2)
            size = 60000;
        else if (size <= 4)
            size = 8000;
        else
            size = 1500;

        an = 1 + n_randint(state, size);
        bn = 1 + n_randint(state, size);
        xn = 1 + n_randint(state, size);
        if (n_randint(state, 4) == 0)
            bn = an;

        flint_set_num_threads(1 + n_randint(state, 4));

        /* one buffer holding b, a, x and room for both outputs after
           them, so that either output can be placed over the inputs it
           may overlap */
        buf = flint_malloc((bn + an + xn + an + xn + bn + xn) * sizeof(ulong));
        b = buf;
        a = b + bn;
        x = a + an;
        r1 = flint_malloc((an + xn) * sizeof(ulong));
        r2 = flint_malloc((bn + xn) * sizeof(ulong));

        if (n_randint(state, 2))
        {
            flint_mpn_rrandom(a, state, an);
            flint_mpn_rrandom(b, state, bn);
            flint_mpn_rrandom(x, state, xn);
        }
        else
        {
            flint_mpn_urandomb(a, state, an * FLINT_BITS);
            flint_mpn_urandomb(b, state, bn * FLINT_BITS);
            flint_mpn_urandomb(x, state, xn * FLINT_BITS);
        }
        /* all-ones operands: the largest products the plan must hold */
        if (n_randint(state, 8) == 0)
        {
            memset(a, 0xff, an * sizeof(ulong));
            memset(b, 0xff, bn * sizeof(ulong));
            memset(x, 0xff, xn * sizeof(ulong));
        }

        flint_mpn_mul(r1, FLINT_MAX(an, xn) == an ? a : x, FLINT_MAX(an, xn),
                          FLINT_MAX(an, xn) == an ? x : a, FLINT_MIN(an, xn));
        flint_mpn_mul(r2, FLINT_MAX(bn, xn) == bn ? b : x, FLINT_MAX(bn, xn),
                          FLINT_MAX(bn, xn) == bn ? x : b, FLINT_MIN(bn, xn));

        alias = n_randint(state, 3);
        if (alias == 1)
        {
            /* z1 over a and x */
            z1 = a;
            z2 = x + xn;
        }
        else if (alias == 2)
        {
            /* z2 over b, a and x */
            z2 = buf;
            z1 = x + xn;
        }
        else
        {
            z1 = x + xn;
            z2 = z1 + an + xn;
        }

        _fft_small_mpn_mul_shared(z1, a, an, z2, b, bn, x, xn);

        if (mpn_cmp(z1, r1, an + xn) != 0 || mpn_cmp(z2, r2, bn + xn) != 0)
        {
            flint_printf("FAIL: products\nan = %wu, bn = %wu, xn = %wu, "
                         "alias = %d\n", an, bn, xn, alias);
            fflush(stdout);
            flint_abort();
        }

        flint_free(buf);
        flint_free(r1);
        flint_free(r2);
    }

    TEST_FUNCTION_END(state);
}
