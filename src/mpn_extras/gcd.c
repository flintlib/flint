/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* {up, *un} = {xp, xn} >> b for 0 <= b < FLINT_BITS, normalised */
static void
_rshift_norm(mp_ptr up, mp_size_t * un, mp_srcptr xp, mp_size_t xn, unsigned int b)
{
    if (b != 0)
    {
        mpn_rshift(up, xp, xn, b);
        xn -= (up[xn - 1] == 0);
    }
    else
        flint_mpn_copyi(up, xp, xn);
    *un = xn;
}

mp_size_t
flint_mpn_gcd(mp_ptr gp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    mp_size_t az, bz, gz, un, vn, gn;
    unsigned int ab, bb, gb;
    mp_ptr up, vp, rp;
    TMP_INIT;

    FLINT_ASSERT(an > 0 && ap[an - 1] != 0);
    FLINT_ASSERT(bn > 0 && bp[bn - 1] != 0);

    if (an == 1)
    {
        gp[0] = mpn_gcd_1(bp, bn, ap[0]);
        return 1;
    }

    if (bn == 1)
    {
        gp[0] = mpn_gcd_1(ap, an, bp[0]);
        return 1;
    }

    /* the common power of two 2^(gz FLINT_BITS + gb) is set aside, and the
       gcd of the odd parts is computed */
    for (az = 0; ap[az] == 0; az++) ;
    for (bz = 0; bp[bz] == 0; bz++) ;
    ab = flint_ctz(ap[az]);
    bb = flint_ctz(bp[bz]);

    if (az != bz)
    {
        gz = FLINT_MIN(az, bz);
        gb = (az < bz) ? ab : bb;
    }
    else
    {
        gz = az;
        gb = FLINT_MIN(ab, bb);
    }

    if (an - az <= 2 && bn - bz <= 2)
    {
        /* both odd parts fit in two limbs */
        mp_limb_t u1, u0, v1, v0, g[2];

        u0 = ap[az];
        u1 = (an - az == 2) ? ap[az + 1] : 0;
        v0 = bp[bz];
        v1 = (bn - bz == 2) ? bp[bz + 1] : 0;

        if (ab != 0)
        {
            u0 = (u0 >> ab) | (u1 << (FLINT_BITS - ab));
            u1 >>= ab;
        }
        if (bb != 0)
        {
            v0 = (v0 >> bb) | (v1 << (FLINT_BITS - bb));
            v1 >>= bb;
        }

        if (u1 == 0 && v1 == 0)
        {
            g[0] = mpn_gcd_11(u0, v0);
            g[1] = 0;
        }
        else if (u1 == 0 || v1 == 0)
        {
            mp_limb_t w[2];
            if (u1 == 0)
            {
                w[0] = v0; w[1] = v1;
                g[0] = mpn_gcd_1(w, 2, u0);
            }
            else
            {
                w[0] = u0; w[1] = u1;
                g[0] = mpn_gcd_1(w, 2, v0);
            }
            g[1] = 0;
        }
        else
        {
            _flint_mpn_gcd_22(g, u1, u0, v1, v0);
        }

        /* g << (gz FLINT_BITS + gb), which fits in min(an, bn) limbs */
        flint_mpn_zero(gp, gz);
        if (gb != 0)
        {
            g[1] = (g[1] << gb) | (g[0] >> (FLINT_BITS - gb));
            g[0] <<= gb;
        }
        gp[gz] = g[0];
        gn = gz + 1;
        if (g[1] != 0)
            gp[gn++] = g[1];
        return gn;
    }

    TMP_START;

    up = TMP_ALLOC((an - az + bn - bz) * sizeof(mp_limb_t));
    vp = up + (an - az);
    _rshift_norm(up, &un, ap + az, an - az, ab);
    _rshift_norm(vp, &vn, bp + bz, bn - bz, bb);

    /* the first operand must be at least as large as the second (by size,
       and for equal sizes by the top limb, which saves a reduction step);
       the gcd is written over the second */
    if (un < vn || (un == vn && up[un - 1] < vp[vn - 1]))
    {
        FLINT_SWAP(mp_ptr, up, vp);
        FLINT_SWAP(mp_size_t, un, vn);
    }

    rp = vp;
    gn = _flint_mpn_gcd_inplace(rp, up, un, vp, vn);

    /* gp = {rp, gn} << (gz FLINT_BITS + gb); this is written last, as gp may
       be ap or bp */
    flint_mpn_zero(gp, gz);
    if (gb != 0)
    {
        mp_limb_t cy = mpn_lshift(gp + gz, rp, gn, gb);
        if (cy != 0)
            gp[gz + gn++] = cy;
    }
    else
    {
        flint_mpn_copyi(gp + gz, rp, gn);
    }

    TMP_END;
    return gz + gn;
}
