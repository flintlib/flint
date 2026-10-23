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

/* z1 = a x and z2 = b x against flint_mpn_mul, over lengths on both
   sides of the transformed-path rule (from single limbs to beyond the
   two-prime range), operands in any order of lengths */
TEST_FUNCTION_START(flint_mpn_mul_shared, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        mp_size_t an, bn, xn, size;
        mp_ptr a, b, x, z1, z2, r1, r2;

        size = n_randint(state, 20);
        size = (size == 0) ? 60000 : (size <= 3) ? 3000
             : (size <= 9) ? 600 : 40;
        if (iter % 50 == 0 && iter != 0)
            size = 100000;

        an = 1 + n_randint(state, size);
        bn = 1 + n_randint(state, size);
        xn = 1 + n_randint(state, size);
        if (n_randint(state, 4) == 0)
            bn = an;

        if (size >= 3000)
            flint_set_num_threads(1 + n_randint(state, 4));

        a = flint_malloc(an * sizeof(mp_limb_t));
        b = flint_malloc(bn * sizeof(mp_limb_t));
        x = flint_malloc(xn * sizeof(mp_limb_t));
        z1 = flint_malloc((an + xn) * sizeof(mp_limb_t));
        z2 = flint_malloc((bn + xn) * sizeof(mp_limb_t));
        r1 = flint_malloc((an + xn) * sizeof(mp_limb_t));
        r2 = flint_malloc((bn + xn) * sizeof(mp_limb_t));

        flint_mpn_rrandom(a, state, an);
        flint_mpn_rrandom(b, state, bn);
        flint_mpn_rrandom(x, state, xn);
        if (n_randint(state, 8) == 0)
        {
            memset(a, 0xff, an * sizeof(mp_limb_t));
            memset(b, 0xff, bn * sizeof(mp_limb_t));
            memset(x, 0xff, xn * sizeof(mp_limb_t));
        }

        if (an >= xn)
            flint_mpn_mul(r1, a, an, x, xn);
        else
            flint_mpn_mul(r1, x, xn, a, an);
        if (bn >= xn)
            flint_mpn_mul(r2, b, bn, x, xn);
        else
            flint_mpn_mul(r2, x, xn, b, bn);

        flint_mpn_mul_shared(z1, a, an, z2, b, bn, x, xn);

        if (mpn_cmp(z1, r1, an + xn) != 0 || mpn_cmp(z2, r2, bn + xn) != 0)
            TEST_FUNCTION_FAIL("an = %wd, bn = %wd, xn = %wd, transformed = %d\n",
                an, bn, xn, _flint_mpn_mul_shared_use_fft(an, bn, xn));

        flint_free(a);
        flint_free(b);
        flint_free(x);
        flint_free(z1);
        flint_free(z2);
        flint_free(r1);
        flint_free(r2);

        flint_set_num_threads(1);
    }

    TEST_FUNCTION_END(state);
}
