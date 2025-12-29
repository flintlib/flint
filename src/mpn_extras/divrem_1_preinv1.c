/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* mpn_divrem_1 -- mpn by limb division.

Copyright 1991, 1993, 1994, 1996, 1998-2000, 2002, 2003 Free Software
Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */


#include "ulong_extras.h"
#include "mpn_extras.h"

mp_limb_t
flint_mpn_divrem_2_1_preinv_norm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv)
{
    mp_limb_t r;
    r = n_divrem_norm(&qp[1], up[1], d);
    udiv_qrnnd_preinv(qp[0], r, r, up[0], d, dinv);
    return r;
}

mp_limb_t
flint_mpn_divrem_3_1_preinv_norm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv)
{
    mp_limb_t r;
    r = n_divrem_norm(&qp[2], up[2], d);
    udiv_qrnnd_preinv(qp[1], r, r, up[1], d, dinv);
    udiv_qrnnd_preinv(qp[0], r, r, up[0], d, dinv);
    return r;
}

mp_limb_t
flint_mpn_divrem_2_1_preinv_unnorm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv, mp_limb_t norm)
{
    mp_limb_t u0, u1, r;

    FLINT_ASSERT(norm >= 1);

    u1 = up[1];
    u0 = up[0];
    if (u1 < d)
    {
        d <<= norm;
        qp[1] = 0;
        r = (u1 << norm) | (u0 >> (FLINT_BITS - norm));
    }
    else
    {
        d <<= norm;
        r = (u1 >> (FLINT_BITS - norm));
        udiv_qrnnd_preinv(qp[1], r, r, (u1 << norm) | (u0 >> (FLINT_BITS - norm)), d, dinv);
    }

    udiv_qrnnd_preinv(qp[0], r, r, u0 << norm, d, dinv);
    return r >> norm;
}

mp_limb_t
flint_mpn_divrem_3_1_preinv_unnorm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv, mp_limb_t norm)
{
    mp_limb_t u0, u1, u2, r;

    FLINT_ASSERT(norm >= 1);

    u2 = up[2];
    u1 = up[1];
    if (u2 < d)
    {
        d <<= norm;
        qp[2] = 0;
        r = (u2 << norm) | (u1 >> (FLINT_BITS - norm));
    }
    else
    {
        d <<= norm;
        r = (u2 >> (FLINT_BITS - norm));
        udiv_qrnnd_preinv(qp[2], r, r, (u2 << norm) | (u1 >> (FLINT_BITS - norm)), d, dinv);
    }

    u0 = up[0];
    udiv_qrnnd_preinv(qp[1], r, r, (u1 << norm) | (u0 >> (FLINT_BITS - norm)), d, dinv);
    udiv_qrnnd_preinv(qp[0], r, r, u0 << norm, d, dinv);
    return r >> norm;
}

mp_limb_t
flint_mpn_divrem_1_preinv(mp_ptr qp, mp_srcptr up, mp_size_t n, mp_limb_t d, mp_limb_t dinv, mp_limb_t norm)
{
    mp_limb_t u2, u3, r;
    slong i;

    if (norm == 0)
    {
        r = n_divrem_norm(&qp[n - 1], up[n - 1], d);
        for (i = n - 2; i >= 0; i--)
            udiv_qrnnd_preinv(qp[i], r, r, up[i], d, dinv);
        return r;
    }
    else
    {
        if (n == 1)
            return n_divrem_preinv_unnorm(qp, up[0], d, dinv, norm);

        u3 = up[n - 1];
        u2 = up[n - 2];
        if (u3 < d)
        {
            d <<= norm;
            qp[n - 1] = 0;
            r = (u3 << norm) | (u2 >> (FLINT_BITS - norm));
        }
        else
        {
            d <<= norm;
            r = (u3 >> (FLINT_BITS - norm));
            udiv_qrnnd_preinv(qp[n - 1], r, r, (u3 << norm) | (u2 >> (FLINT_BITS - norm)), d, dinv);
        }

        for (i = n - 2; i >= 1; i--)
        {
            u3 = up[i - 1];
            udiv_qrnnd_preinv(qp[i], r, r, (u2 << norm) | (u3 >> (FLINT_BITS - norm)), d, dinv);
            u2 = u3;
        }

        udiv_qrnnd_preinv(qp[0], r, r, u2 << norm, d, dinv);
        return r >> norm;
    }
}

