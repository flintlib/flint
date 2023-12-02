/*
    Copyright (C) 2010, 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

static
slong NMOD_DIVREM_BC_ITCH(slong lenA, slong lenB, nmod_t mod)
{
    const flint_bitcnt_t bits =
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);

    if (bits <= FLINT_BITS)
        return lenA;
    else if (bits <= 2 * FLINT_BITS)
        return 2*(lenA + lenB - 1);
    else
        return 3*(lenA + lenB - 1);
}


void _nmod_poly_divrem_q0_preinv1(mp_ptr Q, mp_ptr R,
                          mp_srcptr A, mp_srcptr B, slong lenA, mp_limb_t invL, nmod_t mod)
{
    if (lenA == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invL, mod);
    }
    else
    {
        Q[0] = nmod_mul(A[lenA-1], invL, mod);

        if (R == A)
        {
            _nmod_vec_scalar_addmul_nmod(R, B, lenA - 1, nmod_neg(Q[0], mod), mod);
        }
        else
        {
            _nmod_vec_scalar_mul_nmod(R, B, lenA - 1, Q[0], mod);
            _nmod_vec_sub(R, A, R, lenA - 1, mod);
        }
    }
}

void _nmod_poly_divrem_q1_preinv1(mp_ptr Q, mp_ptr R,
                          mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                          mp_limb_t invL, nmod_t mod)
{
    if (lenB == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invL, mod);
    }
    else
    {
        mp_limb_t q0, q1, t, t0, t1, s0, s1;
        slong i;

        q1 = nmod_mul(A[lenA-1], invL, mod);
        t = nmod_mul(q1, B[lenB-2], mod);
        t = nmod_sub(t, A[lenA-2], mod);
        q0 = nmod_mul(t, invL, mod);

        R[0] = nmod_addmul(A[0], q0, B[0], mod);

        Q[0] = nmod_neg(q0, mod);
        Q[1] = q1;
        q1 = nmod_neg(q1, mod);

        if (mod.norm >= (FLINT_BITS + 1) / 2 + 1)
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
}

void
_nmod_poly_divrem_basecase_preinv1_1(mp_ptr Q, mp_ptr R, mp_ptr W,
                             mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                             mp_limb_t invL,
                             nmod_t mod)
{
    slong iR;
    mp_ptr ptrQ = Q - lenB + 1;
    mp_ptr R1 = W;

    flint_mpn_copyi(R1, A, lenA);

    for (iR = lenA - 1; iR >= lenB - 1; iR--)
    {
        if (R1[iR] == 0)
        {
            ptrQ[iR] = WORD(0);
        }
        else
        {
            ptrQ[iR] = n_mulmod2_preinv(R1[iR], invL, mod.n, mod.ninv);

            if (lenB > 1)
            {
                const mp_limb_t c = n_negmod(ptrQ[iR], mod.n);
                mpn_addmul_1(R1 + iR - lenB + 1, B, lenB - 1, c);
            }
        }
    }

    if (lenB > 1)
        _nmod_vec_reduce(R, R1, lenB - 1, mod);
}

void
_nmod_poly_divrem_basecase_preinv1_2(mp_ptr Q, mp_ptr R, mp_ptr W,
                             mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                             mp_limb_t invL,
                             nmod_t mod)
{
    slong iR, i;
    mp_ptr B2 = W, R2 = W + 2*(lenB - 1), ptrQ = Q - lenB + 1;

    for (i = 0; i < lenB - 1; i++)
    {
        B2[2 * i] = B[i];
        B2[2 * i + 1] = 0;
    }
    for (i = 0; i < lenA; i++)
    {
        R2[2 * i] = A[i];
        R2[2 * i + 1] = 0;
    }

    for (iR = lenA - 1; iR >= lenB - 1; )
    {
        mp_limb_t r =
            n_ll_mod_preinv(R2[2 * iR + 1], R2[2 * iR], mod.n, mod.ninv);

        while ((iR + 1 >= lenB) && (r == WORD(0)))
        {
            ptrQ[iR--] = WORD(0);
            if (iR + 1 >= lenB)
                r = n_ll_mod_preinv(R2[2 * iR + 1], R2[2 * iR], mod.n,
                                    mod.ninv);
        }

        if (iR + 1 >= lenB)
        {
            ptrQ[iR] = n_mulmod2_preinv(r, invL, mod.n, mod.ninv);

            if (lenB > 1)
            {
                const mp_limb_t c = n_negmod(ptrQ[iR], mod.n);
                mpn_addmul_1(R2 + 2 * (iR - lenB + 1), B2, 2 * lenB - 2, c);
            }
            iR--;
        }
    }

    for (iR = 0; iR < lenB - 1; iR++)
        R[iR] = n_ll_mod_preinv(R2[2*iR+1], R2[2*iR], mod.n, mod.ninv);
}

void
_nmod_poly_divrem_basecase_preinv1_3(mp_ptr Q, mp_ptr R, mp_ptr W,
                                     mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                                     mp_limb_t invL,
                                     nmod_t mod)
{
    slong iR, i;
    mp_ptr B3 = W, R3 = W + 3*(lenB - 1), ptrQ = Q - lenB + 1;

    for (i = 0; i < lenB - 1; i++)
    {
        B3[3 * i] = B[i];
        B3[3 * i + 1] = 0;
        B3[3 * i + 2] = 0;
    }
    for (i = 0; i < lenA; i++)
    {
        R3[3 * i] = A[i];
        R3[3 * i + 1] = 0;
        R3[3 * i + 2] = 0;
    }

    for (iR = lenA - 1; iR >= lenB - 1; )
    {
        mp_limb_t r =
            n_lll_mod_preinv(R3[3 * iR + 2], R3[3 * iR + 1],
                             R3[3 * iR], mod.n, mod.ninv);

        while ((iR + 1 >= lenB) && (r == WORD(0)))
        {
            ptrQ[iR--] = WORD(0);
            if (iR + 1 >= lenB)
                r = n_lll_mod_preinv(R3[3 * iR + 2], R3[3 * iR + 1],
                                     R3[3 * iR], mod.n, mod.ninv);
        }

        if (iR + 1 >= lenB)
        {
            ptrQ[iR] = n_mulmod2_preinv(r, invL, mod.n, mod.ninv);

            if (lenB > 1)
            {
                const mp_limb_t c = n_negmod(ptrQ[iR], mod.n);
                mpn_addmul_1(R3 + 3 * (iR - lenB + 1), B3, 3 * lenB - 3, c);
            }
            iR--;
        }
    }

    for (iR = 0; iR < lenB - 1; iR++)
        R[iR] = n_lll_mod_preinv(R3[3 * iR + 2], R3[3 * iR + 1],
                                 R3[3 * iR], mod.n, mod.ninv);
}

void
_nmod_poly_divrem_basecase_preinv1(mp_ptr Q, mp_ptr R,
                           mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                           mp_limb_t invB,
                           nmod_t mod)
{
    if (lenA == lenB + 1)
    {
        _nmod_poly_divrem_q1_preinv1(Q, R, A, lenA, B, lenB, invB, mod);
    }
    else if (lenA == lenB)
    {
        _nmod_poly_divrem_q0_preinv1(Q, R, A, B, lenB, invB, mod);
    }
    else if (lenB == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invB, mod);
    }
    else
    {
        mp_ptr W;
        TMP_INIT;
        slong bits = 2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);

        TMP_START;
        W = TMP_ALLOC(NMOD_DIVREM_BC_ITCH(lenA, lenB, mod)*sizeof(mp_limb_t));

        if (bits <= FLINT_BITS)
            _nmod_poly_divrem_basecase_preinv1_1(Q, R, W, A, lenA, B, lenB, invB, mod);
        else if (bits <= 2 * FLINT_BITS)
            _nmod_poly_divrem_basecase_preinv1_2(Q, R, W, A, lenA, B, lenB, invB, mod);
        else
            _nmod_poly_divrem_basecase_preinv1_3(Q, R, W, A, lenA, B, lenB, invB, mod);

        TMP_END;
    }
}

void
_nmod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA,
                                  mp_srcptr B, slong lenB, nmod_t mod)
{
    mp_limb_t invB;

    invB = (B[lenB - 1] == 1) ? 1 : n_invmod(B[lenB - 1], mod.n);
    _nmod_poly_divrem_basecase_preinv1(Q, R, A, lenA, B, lenB, invB, mod);
}

void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R,
                      const nmod_poly_t A, const nmod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    nmod_poly_t tQ, tR;
    mp_ptr q, r;

    if (lenB == 0)
    {
        if (nmod_poly_modulus(B) == 1)
        {
            nmod_poly_set(Q, A);
            nmod_poly_zero(R);
            return;
        } else
        {
            flint_throw(FLINT_ERROR, "Exception (nmod_poly_divrem). Division by zero.");
        }
    }

    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2_preinv(tQ, A->mod.n, A->mod.ninv, lenA - lenB + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
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

    _nmod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, tQ);
        nmod_poly_clear(tQ);
    }
    if (R == B)
    {
        nmod_poly_swap(R, tR);
        nmod_poly_clear(tR);
    }

    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}
