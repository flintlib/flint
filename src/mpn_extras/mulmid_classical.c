/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* Propagate a single-limb carry cc into z[k..), extending the valid prefix
   *top by writing fresh limbs as needed, and dropping any carry that would
   land at index >= zn (i.e. at a limb position >= the top of the window). */
static void
_flint_mpn_mulmid_carry(mp_ptr z, mp_size_t * top, mp_size_t zn,
                        mp_size_t k, mp_limb_t cc)
{
    while (cc != 0 && k < zn)
    {
        if (k < *top)
        {
            mp_limb_t t = z[k] + cc;
            cc = (t < z[k]);
            z[k] = t;
            k++;
        }
        else
        {
            z[k] = cc;
            *top = k + 1;
            cc = 0;
        }
    }
}

/*
    Set (z, zhi - zlo) to the limb window [zlo, zhi) of the integer product
    a * b, as a lower approximation: partial products a[p] * b[q] landing
    entirely below limb zlo are dropped, so the carry they would propagate
    into limb zlo is not recovered.  This matches the contract of
    radix_mulmid_classical with radix B = 2^64.

    The low-end deficit (exact window minus computed window) is bounded by
    min(an, bn, zlo) * 2^64; limbs sufficiently above zlo are exact, up to the
    truncation at zhi.  With zlo == 0 the result is exact.

    Requires an >= 1, bn >= 1 and 0 <= zlo < zhi <= an + bn.

    Schoolbook: each row a * b[i] is accumulated with mpn_mul_1 / mpn_addmul_1
    over only the limbs that fall inside the window.  No scratch space is used
    and z is never pre-zeroed: the first write to each output limb is a
    mpn_mul_1 or a carry store, and later rows accumulate on top.
*/
void
flint_mpn_mulmid_classical(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                 mp_size_t zlo, mp_size_t zhi)
{
    mp_size_t i, o, s, e, len, zn, i0, top;
    mp_limb_t cc;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(zlo >= 0);
    FLINT_ASSERT(zhi > zlo);
    FLINT_ASSERT(zhi <= an + bn);

    if (an < bn)
    {
        FLINT_SWAP(mp_srcptr, a, b);
        FLINT_SWAP(mp_size_t, an, bn);
    }

    /* limbs of either input at or above zhi cannot reach the window */
    an = FLINT_MIN(an, zhi);
    bn = FLINT_MIN(bn, zhi);

    zn = zhi - zlo;

    /*
        Phase 1: the rows i in [0, min(zlo, bn-1)] all land at output offset 0
        (their low parts a[0 .. zlo-i) are below zlo and dropped).  These rows
        have non-decreasing length in i, so the largest i is the longest; do it
        first with mpn_mul_1 to seed the prefix, then add the shorter ones on
        top (whose carries propagate within the interior).
    */
    i0 = FLINT_MIN(zlo, bn - 1);

    s = zlo - i0;
    e = FLINT_MIN(an, zhi - i0);
    len = e - s;
    if (len > 0)
    {
        cc = mpn_mul_1(z, a + s, len, b[i0]);
        top = len;
        if (len < zn)
        {
            z[len] = cc;
            top = len + 1;
        }
    }
    else
    {
        /* the longest phase-1 row is already empty, so every product lands
           below zlo: the window starts out all zero */
        z[0] = 0;
        top = 1;
    }

    for (i = i0 - 1; i >= 0; i--)
    {
        s = zlo - i;
        e = FLINT_MIN(an, zhi - i);
        len = e - s;
        if (len <= 0)
            break;                  /* this and all smaller i are below window */
        cc = mpn_addmul_1(z, a + s, len, b[i]);
        _flint_mpn_mulmid_carry(z, &top, zn, len, cc);
    }

    /*
        Phase 2: rows i in (zlo, bn) land at offset i - zlo, sliding right by one
        limb per row.  Each row's region is interior up to the current frontier
        and extends it by at most one limb (written via the carry), exactly as
        in a plain mpn multiplication.
    */
    for (i = zlo + 1; i < bn; i++)
    {
        o = i - zlo;
        len = FLINT_MIN(an, zhi - i);       /* s == 0 */
        if (len <= 0)
            break;                          /* zhi <= i: no further rows reach */

        if (o + len <= top)
        {
            /* the whole row sits inside the current prefix */
            cc = mpn_addmul_1(z + o, a, len, b[i]);
            _flint_mpn_mulmid_carry(z, &top, zn, o + len, cc);
        }
        else
        {
            /* [o, top) overlaps the prefix; [top, o+len) is fresh */
            mp_size_t ov = top - o;         /* overlap length, >= 0 */
            mp_size_t right = o + len;
            mp_limb_t ccf = mpn_mul_1(z + top, a + ov, len - ov, b[i]);
            mp_size_t topf = right;

            if (right < zn)
            {
                z[right] = ccf;
                topf = right + 1;
            }

            if (ov > 0)
            {
                cc = mpn_addmul_1(z + o, a, ov, b[i]);
                top = topf;
                _flint_mpn_mulmid_carry(z, &top, zn, o + ov, cc);
            }
            else
            {
                top = topf;
            }
        }
    }

    /* Any high limbs never reached by a carry are exact zeros: with zhi ==
       an + bn the very top limb can be one of these, so fill explicitly
       rather than asserting the frontier reached the top. */
    if (top < zn)
        flint_mpn_zero(z + top, zn - top);
}
