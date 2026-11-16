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
#include "nmod_poly/impl.h"

/* a modulus, often at the boundaries of the kernel's cases: up to 2^52,
   low halves alone and u52_red up to 2^26, REDC for odd n above; above
   2^52, REDC for odd n, NMOD_RED3 for even n; a mask for powers of 2.
   On 32-bit builds (no kernel: the classical code), those below 2^32. */
static ulong
_ifma_modulus(flint_rand_t state)
{
    static const ulong special[] = {
        1, 2, 3, 1023, 1024, 1025, 1027,
        (UWORD(1) << 26) - 5, UWORD(1) << 26, (UWORD(1) << 26) + 1,
#if FLINT_BITS == 64
        (UWORD(1) << 32) - 5, UWORD(1) << 40,
        (UWORD(1) << 51) + 1, (UWORD(1) << 52) - 1, UWORD(1) << 52,
        (UWORD(1) << 52) + 1, (UWORD(1) << 52) + 2, (UWORD(1) << 53) - 1,
        UWORD(1) << 53,
#endif
        UWORD(1) << (FLINT_BITS - 1), (UWORD(1) << (FLINT_BITS - 1)) + 2,
        UWORD_MAX - 58, UWORD_MAX - 1, UWORD_MAX };
    ulong n;

    if (n_randint(state, 4) == 0)
        return special[n_randint(state, sizeof(special) / sizeof(ulong))];

    n = n_randtest_bits(state, 1 + n_randint(state, FLINT_BITS));
    return FLINT_MAX(n, 1);
}

/* random residues, many n - 1 */
static void
_ifma_randvec(nn_ptr v, slong len, flint_rand_t state, ulong n)
{
    slong i;
    int all_max = (n_randint(state, 4) == 0);

    for (i = 0; i < len; i++)
        v[i] = (all_max || n_randint(state, 4) == 0) ? n - 1 : n_randint(state, n);
}

TEST_FUNCTION_START(nmod_poly_mulmid_ifma, state)
{
    slong iter;

    /* any window, against the classical algorithm, with the CRT variant
       for n > 2^52 one time out of two; the coefficients after the window
       are not written */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        nn_ptr a, b, c, d;
        slong la, lb, len, nlo, nhi, i;
        int crt;

        nmod_init(&mod, _ifma_modulus(state));

        la = 1 + n_randint(state, n_randint(state, 4) ? 40 : 511);
        lb = la + n_randint(state, n_randint(state, 8) ? 300 : 3000);
        if (n_randint(state, 2))
            FLINT_SWAP(slong, la, lb);
        len = la + lb - 1;

        switch (n_randint(state, 5))
        {
            case 0: nlo = 0; nhi = len; break;
            case 1: nlo = 0; nhi = 1 + n_randint(state, len); break;
            case 2: nlo = FLINT_MIN(la, lb) - 1; nhi = FLINT_MAX(la, lb); break;
            case 3: nlo = n_randint(state, len); nhi = len; break;
            default:
                nlo = n_randint(state, len);
                nhi = nlo + 1 + n_randint(state, len - nlo);
        }

        a = _nmod_vec_init(la);
        b = _nmod_vec_init(lb);
        c = _nmod_vec_init(nhi - nlo + 1);
        d = _nmod_vec_init(nhi - nlo);

        _ifma_randvec(a, la, state, mod.n);
        _ifma_randvec(b, lb, state, mod.n);
        c[nhi - nlo] = 12345;

        crt = n_randint(state, 2);
        if (crt)
            _nmod_poly_mulmid_ifma_crt(c, a, la, b, lb, nlo, nhi, mod);
        else
            _nmod_poly_mulmid_ifma(c, a, la, b, lb, nlo, nhi, mod);
        _nmod_poly_mulmid_classical(d, a, la, b, lb, nlo, nhi, mod);

        if (!_nmod_vec_equal(c, d, nhi - nlo) || c[nhi - nlo] != 12345)
        {
            for (i = 0; i < nhi - nlo && c[i] == d[i]; i++) ;
            flint_printf("FAIL (window):\n");
            flint_printf("n = %wu, la = %wd, lb = %wd, nlo = %wd, nhi = %wd, crt = %d, first difference at %wd\n",
                         mod.n, la, lb, nlo, nhi, crt, i);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(c);
        _nmod_vec_clear(d);
    }

#if NMOD_POLY_HAVE_IFMA && !FLINT_WANT_ASSERT
    /* words >= n (up to 2^64 - 1) in the coefficients t and above of both
       operands of a full product: the outputs below t, which do not read
       them, are unchanged (the kernel is lane-local); not with asserts,
       which reject such words in the multiplication of the shorter operand
       by 2^104 mod n (odd n > 2^26), and not without the kernel (the
       classical code takes reduced inputs) */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        nn_ptr a, b, c, d;
        slong la, lb, t, i;
        int crt;

        nmod_init(&mod, _ifma_modulus(state));

        la = 1 + n_randint(state, 511);
        lb = la + n_randint(state, 600);
        if (n_randint(state, 2))
            FLINT_SWAP(slong, la, lb);
        t = 1 + n_randint(state, FLINT_MIN(la, lb));

        a = _nmod_vec_init(la);
        b = _nmod_vec_init(lb);
        c = _nmod_vec_init(la + lb - 1);
        d = _nmod_vec_init(t);

        _ifma_randvec(a, la, state, mod.n);
        _ifma_randvec(b, lb, state, mod.n);
        _nmod_poly_mulmid_classical(d, a, t, b, t, 0, t, mod);

        for (i = t; i < la; i++)
            a[i] = n_randint(state, 2) ? UWORD_MAX - n_randint(state, 5) : n_randtest(state);
        for (i = t; i < lb; i++)
            b[i] = n_randint(state, 2) ? UWORD_MAX - n_randint(state, 5) : n_randtest(state);

        crt = n_randint(state, 2);
        if (crt)
            _nmod_poly_mulmid_ifma_crt(c, a, la, b, lb, 0, la + lb - 1, mod);
        else
            _nmod_poly_mulmid_ifma(c, a, la, b, lb, 0, la + lb - 1, mod);

        if (!_nmod_vec_equal(c, d, t))
        {
            flint_printf("FAIL (words >= n):\n");
            flint_printf("n = %wu, la = %wd, lb = %wd, t = %wd, crt = %d\n", mod.n, la, lb, t, crt);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(c);
        _nmod_vec_clear(d);
    }
#endif

    TEST_FUNCTION_END(state);
}
