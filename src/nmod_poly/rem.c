/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_poly.h"

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

void _nmod_poly_rem(mp_ptr R, mp_srcptr A, slong lenA,
                              mp_srcptr B, slong lenB, nmod_t mod)
{
    if (lenA - lenB == 1)
    {
        _nmod_poly_rem_q1(R, A, lenA, B, lenB, mod);
    }
    else if (lenB >= 2)
    {
        mp_ptr Q;
        TMP_INIT;

        TMP_START;
        Q = TMP_ALLOC((lenA - lenB + 1) * sizeof(mp_limb_t));
        _nmod_poly_divrem(Q, R, A, lenA, B, lenB, mod);
        TMP_END;
    }
}

void nmod_poly_rem(nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    nmod_poly_t tR;
    mp_ptr r;

    if (lenB == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_rem). Division by zero.\n");
    }
    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        return;
    }

    if (R == B)
    {
        nmod_poly_init2_preinv(tR, B->mod.n, B->mod.ninv, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_rem(r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (R == B)
    {
        nmod_poly_swap(R, tR);
        nmod_poly_clear(tR);
    }

    R->length = lenB - 1;
    _nmod_poly_normalise(R);
}
