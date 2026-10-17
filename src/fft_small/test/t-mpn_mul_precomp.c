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

/* products by a precomputed transform against flint_mpn_mul, in both
   geometries: several partners per transform, of all lengths up to
   yn_max, with z over y half of the time and the storage given (of the
   size _size reports, in a block holding nothing else) or allocated */
TEST_FUNCTION_START(_fft_small_mpn_mul_precomp, state)
{
    for (slong iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        ulong xn, ymax, size, sz, yn;
        nn_ptr x, y, z, r;
        double * storage = NULL;
        fft_small_mpn_mul_precomp_t P;
        int j;

        /* both geometries: two-prime up to FLINT_MPN_MUL_NP2_MAX_AN
           limbs, wide beyond */
        size = n_randint(state, 12);
        size = (size == 0) ? 150000 : (size <= 2) ? 60000
             : (size <= 4) ? 8000 : 1500;

        xn = 1 + n_randint(state, size);
        ymax = 1 + n_randint(state, size);

        flint_set_num_threads(1 + n_randint(state, 4));

        x = flint_malloc(xn * sizeof(ulong));
        y = flint_malloc((xn + ymax + 1) * sizeof(ulong));
        z = flint_malloc((xn + ymax + 1) * sizeof(ulong));
        r = flint_malloc((xn + ymax + 1) * sizeof(ulong));

        if (n_randint(state, 2))
            flint_mpn_rrandom(x, state, xn);
        else
            flint_mpn_urandomb(x, state, xn * FLINT_BITS);
        if (n_randint(state, 8) == 0)
            memset(x, 0xff, xn * sizeof(ulong));

        sz = _fft_small_mpn_mul_precomp_size(xn, ymax);
        if (sz == 0)
        {
            flint_printf("FAIL: size\nxn = %wu, ymax = %wu\n", xn, ymax);
            fflush(stdout);
            flint_abort();
        }
        if (n_randint(state, 2))
            storage = flint_aligned_alloc(64, sz * sizeof(double));

        _fft_small_mpn_mul_precomp_init(P, x, xn, ymax, storage);

        for (j = 0; j < 4; j++)
        {
            int alias = n_randint(state, 2);
            nn_ptr zz;

            yn = (j == 0) ? ymax : 1 + n_randint(state, ymax);

            if (n_randint(state, 2))
                flint_mpn_rrandom(y, state, yn);
            else
                flint_mpn_urandomb(y, state, yn * FLINT_BITS);
            if (n_randint(state, 8) == 0)
                memset(y, 0xff, yn * sizeof(ulong));

            if (yn >= xn)
                flint_mpn_mul(r, y, yn, x, xn);
            else
                flint_mpn_mul(r, x, xn, y, yn);

            zz = alias ? y : z;
            _fft_small_mpn_mul_precomp(zz, P, y, yn);
            if (mpn_cmp(zz, r, xn + yn) != 0)
            {
                flint_printf("FAIL: product\nxn = %wu, ymax = %wu, yn = %wu, "
                             "alias = %d\n", xn, ymax, yn, alias);
                fflush(stdout);
                flint_abort();
            }
        }

        _fft_small_mpn_mul_precomp_clear(P);
        if (storage != NULL)
            flint_aligned_free(storage);
        flint_free(x);
        flint_free(y);
        flint_free(z);
        flint_free(r);
    }

    TEST_FUNCTION_END(state);
}
