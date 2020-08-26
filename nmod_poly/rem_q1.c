/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void _nmod_poly_rem_q1(mp_ptr R, 
                       mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                       nmod_t mod)
{
    slong i;
    mp_limb_t invL, t, q0, q1, t1, t0, s1, s0;

    FLINT_ASSERT(lenA == lenB + 1);
    invL = (B[lenB-1] == 1) ? 1 : n_invmod(B[lenB-1], mod.n);

    if (lenB < 2)
        return;

    q1 = n_mulmod2_preinv(A[lenA-1], invL, mod.n, mod.ninv);
    t  = n_mulmod2_preinv(q1, B[lenB-2], mod.n, mod.ninv);
    t  = n_submod(t, A[lenA-2], mod.n);
    q0 = n_mulmod2_preinv(t, invL, mod.n, mod.ninv);
    q1 = nmod_neg(q1, mod);

    /* R = A + (q1*x + q0)*B */
    t = A[0];
    NMOD_ADDMUL(t, q0, B[0], mod);
    R[0] = t;

    if (mod.norm >= (FLINT_BITS + 1)/2 + 1)
    {
        for (i = 1; i < lenB - 1; i++)
        {
            NMOD_RED2(R[i], 0, A[i] + q1*B[i - 1] + q0*B[i], mod);
        }
    }
    else if (mod.norm != 0)
    {
        for (i = 1; i < lenB - 1; i++)
        {
            umul_ppmm(t1, t0, q1, B[i - 1]);
            umul_ppmm(s1, s0, q0, B[i]);
            add_ssaaaa(t1, t0, t1, t0, 0, A[i]);
            add_ssaaaa(t1, t0, t1, t0, s1, s0);
            t1 = FLINT_MIN(t1, t1 - mod.n);
            FLINT_ASSERT(t1 < mod.n);
            NMOD_RED2(R[i], t1, t0, mod);
        }
    }
    else
    {
        for (i = 1; i < lenB - 1; i++)
        {
            t = A[i];
            NMOD_ADDMUL(t, q1, B[i - 1], mod);
            NMOD_ADDMUL(t, q0, B[i], mod);
            R[i] = t;
        }
    }
}

