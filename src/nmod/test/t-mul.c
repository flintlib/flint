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

/* nmod_mul, _nmod_mul_fullword and _nmod_mul_halflimb against n_mulmod2 */
TEST_FUNCTION_START(nmod_mul, state)
{
    slong i;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        ulong n, a, b, r, s;
        flint_bitcnt_t bits = 1 + n_randint(state, FLINT_BITS);
        int halflimb;

        /* n of the given size, sometimes 2^bits - 1 or 2^(FLINT_BITS/2) */
        switch (n_randint(state, 8))
        {
            case 0:
                n = (bits == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << bits) - 1;
                break;
            case 1:
                n = UWORD(1) << (FLINT_BITS / 2);
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
    }

    TEST_FUNCTION_END(state);
}
