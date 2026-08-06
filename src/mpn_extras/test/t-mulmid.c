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

/*
    flint_mpn_mulmid returns the window [zlo, zhi) of a*b as a lower
    approximation: computed = exact - D, where the deficit D is a single carry
    from the products dropped below zlo, bounded by min(an,bn,zlo)*2^64 (a couple
    of backends, e.g. via_mulhigh, use a different but comparably small D).  Since
    D is subtracted modulo 2^(64*zn), a narrow window can wrap so that the result
    exceeds the exact window numerically -- that is allowed.  What we can check is
    that the deficit d = (exact - computed) mod 2^(64*zn) is confined to the low
    two limbs (it fits there whenever the window has room), and that with
    zlo == 0, where nothing is dropped, the window is exact.
*/
TEST_FUNCTION_START(flint_mpn_mulmid, state)
{
    slong ix;

    for (ix = 0; ix < 30000 * flint_test_multiplier(); ix++)
    {
        mp_size_t an, bn, zn, zlo, zhi, k;
        mp_ptr a, b, z, full, d;

        an = 1 + n_randint(state, 30);
        bn = 1 + n_randint(state, 30);
        if (n_randint(state, 1000) == 0)
        {
            an = 1 + n_randint(state, 600);
            bn = 1 + n_randint(state, 600);
        }

        zlo = n_randint(state, an + bn);
        zhi = zlo + 1 + n_randint(state, an + bn - zlo);
        zn = zhi - zlo;

        a = flint_malloc(sizeof(mp_limb_t) * an);
        b = flint_malloc(sizeof(mp_limb_t) * bn);
        z = flint_malloc(sizeof(mp_limb_t) * zn);
        d = flint_malloc(sizeof(mp_limb_t) * zn);
        full = flint_malloc(sizeof(mp_limb_t) * (an + bn));

        flint_mpn_rrandom(a, state, an);
        flint_mpn_rrandom(b, state, bn);
        flint_mpn_rrandom(z, state, zn);        /* poison */

        flint_mpn_mulmid(z, a, an, b, bn, zlo, zhi);

        if (an >= bn)
            mpn_mul(full, a, an, b, bn);
        else
            mpn_mul(full, b, bn, a, an);

        /* d = (exact - computed) mod 2^(64*zn) = deficit (mod 2^(64*zn)) */
        mpn_sub_n(d, full + zlo, z, zn);

        if (zn >= 2)
        {
            for (k = 2; k < zn; k++)
                if (d[k] != 0)
                    TEST_FUNCTION_FAIL(
                            "deficit exceeds the low two window limbs\n"
                            "ix = %wd, (an, bn) = (%wd, %wd), (zlo, zhi) = (%wd, %wd)\n"
                            "a = %{ulong*}\nb = %{ulong*}\n"
                            "exact = %{ulong*}\ngot   = %{ulong*}\n",
                            ix, an, bn, zlo, zhi, a, an, b, bn,
                            full + zlo, zn, z, zn);

            if (d[1] > (mp_limb_t) (an + bn))
                TEST_FUNCTION_FAIL(
                        "deficit too large (second limb %wu > an + bn)\n"
                        "ix = %wd, (an, bn) = (%wd, %wd), (zlo, zhi) = (%wd, %wd)\n"
                        "a = %{ulong*}\nb = %{ulong*}\n"
                        "exact = %{ulong*}\ngot   = %{ulong*}\n",
                        d[1], ix, an, bn, zlo, zhi, a, an, b, bn,
                        full + zlo, zn, z, zn);
        }

        if (zlo == 0)
            for (k = 0; k < zn; k++)
                if (d[k] != 0)
                    TEST_FUNCTION_FAIL(
                            "zlo == 0 window is not exact\n"
                            "ix = %wd, (an, bn) = (%wd, %wd), zhi = %wd\n"
                            "a = %{ulong*}\nb = %{ulong*}\n"
                            "exact = %{ulong*}\ngot   = %{ulong*}\n",
                            ix, an, bn, zhi, a, an, b, bn,
                            full, zhi, z, zhi);

        flint_free(a);
        flint_free(b);
        flint_free(z);
        flint_free(d);
        flint_free(full);
    }

    TEST_FUNCTION_END(state);
}
