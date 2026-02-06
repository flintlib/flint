/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "radix.h"

FLINT_DLL extern int radix_mulmid_force_threading;

TEST_FUNCTION_START(radix_mulmid_fft_small, state)
{
#if FLINT_HAVE_FFT_SMALL
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c, d;
        slong an, bn, lo, hi;
        int squaring;

        radix_mulmid_force_threading = n_randint(state, 2);
        flint_set_num_threads(1 + n_randint(state, 8));

        radix_init_randtest(radix, state);

        if (n_randint(state, 10) == 0)
        {
            bn = 1 + n_randint(state, 30000);
            an = bn + n_randint(state, 10000);
        }
        else
        {
            bn = 1 + n_randint(state, 30);
            an = bn + n_randint(state, 50);
        }

        squaring = n_randint(state, 2);
        if (squaring)
            an = bn;

        lo = n_randint(state, an + bn);
        hi = lo + 1 + n_randint(state, an + bn - lo);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(bn * sizeof(ulong));
        c = flint_malloc((hi - lo) * sizeof(ulong));
        d = flint_malloc((hi - lo) * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(b, state, bn, radix);
        radix_randtest_limbs(c, state, hi - lo, radix);
        radix_randtest_limbs(d, state, hi - lo, radix);

        if (squaring)
        {
            radix_mulmid_fft_small(c, a, an, a, an, lo, hi, radix);
            radix_mulmid_naive(d, a, an, a, an, lo, hi, radix);
        }
        else
        {
            radix_mulmid_fft_small(c, a, an, b, bn, lo, hi, radix);
            radix_mulmid_naive(d, a, an, b, bn, lo, hi, radix);
        }

        if (mpn_cmp(c, d, hi - lo) != 0)
        {
            flint_printf("FAIL: mul\n");
            flint_printf("an = %wd, bn = %wd, lo = %wd, hi = %wd, squaring = %d\n", an, bn, lo, hi, squaring);
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", b, bn);
            flint_printf("%{ulong*}\n", c, hi - lo);
            flint_printf("%{ulong*}\n", d, hi - lo);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(c);
        flint_free(d);
    }

#endif

    TEST_FUNCTION_END(state);
}
