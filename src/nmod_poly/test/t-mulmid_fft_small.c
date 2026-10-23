/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

/* general windows [znlo, zn) of products, operands in either order and
   with coefficients that no output uses, against the full product */
TEST_FUNCTION_START(nmod_poly_mulmid_fft_small, state)
{
#if FLINT_HAVE_FFT_SMALL
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        ulong m;
        slong an, bn, full, znlo, zn, L;
        nn_ptr a, b, z, w;

        switch (n_randint(state, 5))
        {
            case 0:
                m = n_randtest_not_zero(state);
                break;
            case 1:
                m = get_default_mpn_ctx()->ffts[n_randint(state, MPN_CTX_NCRTS)].mod.n;
                break;
            case 2:
                /* primes for which the transforms can be done directly */
                m = n_randint(state, 2) ? UWORD(998244353) : UWORD(7340033);
                break;
            default:
                m = n_randprime(state, 2 + n_randint(state, FLINT_BITS - 1), 1);
        }
        m = FLINT_MAX(m, UWORD(2));
        nmod_init(&mod, m);

        L = n_randint(state, 8) ? 400 : 4000;
        an = 1 + n_randint(state, L);
        bn = 1 + n_randint(state, L);
        full = an + bn - 1;

        switch (n_randint(state, 3))
        {
            case 0:
                /* a middle product of the longer by the shorter operand,
                   possibly with unused low coefficients */
                znlo = FLINT_MIN(an, bn) - 1 + n_randint(state, 3) * n_randint(state, 20);
                znlo = FLINT_MIN(znlo, full - 1);
                zn = znlo + 1 + n_randint(state, FLINT_MAX(an, bn) - FLINT_MIN(an, bn) + 1);
                break;
            case 1:
                /* a short window */
                znlo = n_randint(state, full);
                zn = znlo + 1 + n_randint(state, 30);
                break;
            default:
                znlo = n_randint(state, full);
                zn = znlo + 1 + n_randint(state, full - znlo);
        }
        zn = FLINT_MIN(zn, full);

        a = _nmod_vec_init(an);
        b = _nmod_vec_init(bn);
        z = _nmod_vec_init(zn - znlo + 1);
        w = _nmod_vec_init(full);

        _nmod_vec_randtest(a, state, an, mod);
        _nmod_vec_randtest(b, state, bn, mod);

        /* the output must not be written past zn - znlo */
        z[zn - znlo] = 31337;

        _nmod_poly_mulmid_fft_small(z, a, an, b, bn, znlo, zn, mod);

        if (an >= bn)
            _nmod_poly_mul_KS(w, a, an, b, bn, mod);
        else
            _nmod_poly_mul_KS(w, b, bn, a, an, mod);

        if (!_nmod_vec_equal(z, w + znlo, zn - znlo) || z[zn - znlo] != 31337)
        {
            flint_printf("FAIL: n = %wu, an = %wd, bn = %wd, znlo = %wd, zn = %wd\n",
                         m, an, bn, znlo, zn);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(z);
        _nmod_vec_clear(w);
    }
#endif

    TEST_FUNCTION_END(state);
}
