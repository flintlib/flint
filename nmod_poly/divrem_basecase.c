/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010, 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_divrem_basecase_1(mp_ptr Q, mp_ptr R, mp_ptr W,
                             mp_srcptr A, len_t lenA, mp_srcptr B, len_t lenB,
                             nmod_t mod)
{
    const mp_limb_t invL = n_invmod(B[lenB - 1], mod.n);
    len_t iR;
    mp_ptr ptrQ = Q - lenB + 1;
    mp_ptr R1 = W;
    
    flint_mpn_copyi(R1, A, lenA);

    for (iR = lenA - 1; iR >= lenB - 1; iR--)
    {
        if (R1[iR] == 0)
        {
            ptrQ[iR] = 0L;
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
_nmod_poly_divrem_basecase_2(mp_ptr Q, mp_ptr R, mp_ptr W,
                             mp_srcptr A, len_t lenA, mp_srcptr B, len_t lenB,
                             nmod_t mod)
{
    const mp_limb_t invL = n_invmod(B[lenB - 1], mod.n);
    len_t iR, i;
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

        while ((iR + 1 >= lenB) && (r == 0L))
        {
            ptrQ[iR--] = 0L;
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
_nmod_poly_divrem_basecase_3(mp_ptr Q, mp_ptr R, mp_ptr W,
                             mp_srcptr A, len_t lenA, mp_srcptr B, len_t lenB,
                             nmod_t mod)
{
    const mp_limb_t invL = n_invmod(B[lenB - 1], mod.n);
    len_t iR, i;
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

        while ((iR + 1 >= lenB) && (r == 0L))
        {
            ptrQ[iR--] = 0L;
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
_nmod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, mp_ptr W,
                           mp_srcptr A, len_t lenA, mp_srcptr B, len_t lenB,
                           nmod_t mod)
{
    const len_t bits =
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);

    if (bits <= FLINT_BITS)
        _nmod_poly_divrem_basecase_1(Q, R, W, A, lenA, B, lenB, mod);
    else if (bits <= 2 * FLINT_BITS)
        _nmod_poly_divrem_basecase_2(Q, R, W, A, lenA, B, lenB, mod);
    else
        _nmod_poly_divrem_basecase_3(Q, R, W, A, lenA, B, lenB, mod);
}

void
nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A,
                          const nmod_poly_t B)
{
    const len_t lenA = A->length, lenB = B->length;
    mp_ptr Q_coeffs, R_coeffs, W;
    nmod_poly_t t1, t2;

    if (lenB == 0)
    {
        printf("Exception (nmod_poly_divrem). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2_preinv(t1, B->mod.n, B->mod.ninv, lenA - lenB + 1);
        Q_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        Q_coeffs = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_poly_init2_preinv(t2, B->mod.n, B->mod.ninv, lenB - 1);
        R_coeffs = t2->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        R_coeffs = R->coeffs;
    }

    W = _nmod_vec_init(NMOD_DIVREM_BC_ITCH(lenA, lenB, A->mod));
    
    _nmod_poly_divrem_basecase(Q_coeffs, R_coeffs, W, A->coeffs, lenA,
                               B->coeffs, lenB, B->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, t1);
        nmod_poly_clear(t1);
    }
    if (R == A || R == B)
    {
        nmod_poly_swap(R, t2);
        nmod_poly_clear(t2);
    }
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_vec_clear(W);
    _nmod_poly_normalise(R);
}
