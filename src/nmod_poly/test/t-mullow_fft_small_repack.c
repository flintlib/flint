/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_mullow_fft_small_repack, state)
{
    slong an, bn, zn;
    nn_ptr a, b, z, w;
    nmod_t mod;

    slong i, j;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_init(&mod, n_randtest_not_zero(state));

        if (n_randint(state, 2))
        {
            an = 1 + n_randtest(state) % 10000;
            bn = 1 + n_randtest(state) % 10000;
            zn = 1 + n_randtest(state) % 10000;
        }
        else
        {
            an = 1 + n_randtest(state) % 1000;
            bn = 1 + n_randtest(state) % 1000;
            zn = 1 + n_randtest(state) % 1000;
        }

        zn = FLINT_MIN(zn, an + bn - 1);

        a = flint_malloc(sizeof(ulong) * an);
        b = flint_malloc(sizeof(ulong) * bn);
        z = flint_malloc(sizeof(ulong) * (zn + 1));
        w = flint_malloc(sizeof(ulong) * zn);

        _nmod_vec_randtest(a, state, an, mod);
        _nmod_vec_randtest(b, state, bn, mod);
        _nmod_vec_randtest(z, state, zn, mod);

        if (n_randint(state, 2))
        {
            for (j = 0; j < an; j++) a[j] = n_randint(state, mod.n);
            for (j = 0; j < bn; j++) b[j] = n_randint(state, mod.n);
        }

        /* Additional check that the unpacking doesn't write one
           coefficient too much */
        z[zn] = 31337;

        if (_nmod_poly_mullow_fft_small_repack(z, a, an, b, bn, zn, mod))
        {
            if (z[zn] != 31337)
                flint_abort();

            if (an >= bn)
                _nmod_poly_mullow_KS(w, a, an, b, bn, 0, zn, mod);
            else
                _nmod_poly_mullow_KS(w, b, bn, a, an, 0, zn, mod);

            if (!_nmod_vec_equal(z, w, zn))
            {
                flint_printf("a = %{ulong*}\n\n", a, an);
                flint_printf("b = %{ulong*}\n\n", b, bn);
                flint_printf("z = %{ulong*}\n\n", z, zn);
                flint_printf("w = %{ulong*}\n\n", w, zn);
                flint_abort();
            }
        }

        flint_free(a);
        flint_free(b);
        flint_free(z);
        flint_free(w);
    }

    TEST_FUNCTION_END(state);
}
