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
#include "nmod.h"

/* nmod_mul, _nmod_mul_fullword, _nmod_mul_halflimb, nmod_addmul, nmod_fmma
   and the bounded variants against n_mulmod2 and n_addmod */
TEST_FUNCTION_START(nmod_mul, state)
{
    slong i;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        ulong n, a, b, c, d, e, r, s, t, u;
        flint_bitcnt_t nu;
        flint_bitcnt_t bits = 1 + n_randint(state, FLINT_BITS);
        int halflimb;

        /* n of the given size, sometimes 2^bits - 1, 2^(bits - 1) or
           2^(FLINT_BITS/2) */
        switch (n_randint(state, 8))
        {
            case 0:
                n = (bits == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << bits) - 1;
                break;
            case 1:
                n = UWORD(1) << (FLINT_BITS / 2);
                break;
            case 2:
                n = UWORD(1) << (bits - 1);
                break;
            default:
                n = n_randbits(state, bits);
        }

        nmod_init(&mod, n);

        halflimb = (n <= (UWORD(1) << (FLINT_BITS / 2)));
        if (!NMOD_IS_HALFLIMB(mod) != !halflimb)
            TEST_FUNCTION_FAIL("NMOD_IS_HALFLIMB\nn = %wu\n", n);

        /* operands, sometimes n - 1 */
        a = n_randint(state, 4) ? n_randtest(state) % n : n - 1;
        b = n_randint(state, 4) ? n_randtest(state) % n : n - 1;

        s = n_mulmod2(a, b, n);

        r = nmod_mul(a, b, mod);
        if (r != s)
            TEST_FUNCTION_FAIL("nmod_mul\nn = %wu\na = %wu\nb = %wu\nr = %wu\ns = %wu\n", n, a, b, r, s);

        if (mod.norm == 0)
        {
            r = _nmod_mul_fullword(a, b, mod);
            if (r != s)
                TEST_FUNCTION_FAIL("_nmod_mul_fullword\nn = %wu\na = %wu\nb = %wu\nr = %wu\ns = %wu\n", n, a, b, r, s);
        }

        if (halflimb)
        {
            r = _nmod_mul_halflimb(a, b, mod);
            if (r != s)
                TEST_FUNCTION_FAIL("_nmod_mul_halflimb\nn = %wu\na = %wu\nb = %wu\nr = %wu\ns = %wu\n", n, a, b, r, s);
        }

        /* s + a b and a b + c d */
        c = n_randint(state, 4) ? n_randtest(state) % n : n - 1;
        d = n_randint(state, 4) ? n_randtest(state) % n : n - 1;
        e = n_randint(state, 4) ? n_randtest(state) % n : n - 1;
        t = n_addmod(e, s, n);
        u = n_addmod(s, n_mulmod2(c, d, n), n);

        r = nmod_addmul(e, a, b, mod);
        if (r != t)
            TEST_FUNCTION_FAIL("nmod_addmul\nn = %wu\ne = %wu\na = %wu\nb = %wu\nr = %wu\nt = %wu\n", n, e, a, b, r, t);

        r = nmod_fmma(a, b, c, d, mod);
        if (r != u)
            TEST_FUNCTION_FAIL("nmod_fmma\nn = %wu\na = %wu\nb = %wu\nc = %wu\nd = %wu\nr = %wu\nu = %wu\n", n, a, b, c, d, r, u);

        nu = FLINT_BIT_COUNT(n);
        if (!NMOD_CAN_USE_MUL_BOUNDED(mod) != !(nu >= 3 && nu <= FLINT_BITS - 2))
            TEST_FUNCTION_FAIL("NMOD_CAN_USE_MUL_BOUNDED\nn = %wu\n", n);
        if (!NMOD_CAN_USE_FMMA_BOUNDED(mod) != !(nu >= 3 && nu <= FLINT_BITS - 3))
            TEST_FUNCTION_FAIL("NMOD_CAN_USE_FMMA_BOUNDED\nn = %wu\n", n);

        if (NMOD_CAN_USE_MUL_BOUNDED(mod))
        {
            r = _nmod_mul_bounded(a, b, mod);
            if (r != s)
                TEST_FUNCTION_FAIL("_nmod_mul_bounded\nn = %wu\na = %wu\nb = %wu\nr = %wu\ns = %wu\n", n, a, b, r, s);

            r = _nmod_addmul_bounded(e, a, b, mod);
            if (r != t)
                TEST_FUNCTION_FAIL("_nmod_addmul_bounded\nn = %wu\ne = %wu\na = %wu\nb = %wu\nr = %wu\nt = %wu\n", n, e, a, b, r, t);
        }

        if (NMOD_CAN_USE_FMMA_BOUNDED(mod))
        {
            r = _nmod_fmma_bounded(a, b, c, d, mod);
            if (r != u)
                TEST_FUNCTION_FAIL("_nmod_fmma_bounded\nn = %wu\na = %wu\nb = %wu\nc = %wu\nd = %wu\nr = %wu\nu = %wu\n", n, a, b, c, d, r, u);
        }
    }

    TEST_FUNCTION_END(state);
}
