/*
    Copyright (C) 2010, 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#  if HAVE_ALLOCA_H
#   include <alloca.h>
#  else
#   if _MSC_VER
#    include <malloc.h>
#    define alloca _alloca
#   else
#    ifdef __DECC
#     define alloca(x) __ALLOCA(x)
#    else
#     ifdef BSD
#      include <stdlib.h>
#     else
#      error Could not find alloca
#     endif
#    endif
#   endif
#  endif
# endif
#endif

#include "gmp.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "flint-impl.h"

void
_nmod_poly_divrem_basecase_1(ulong_ptr Q, ulong_ptr R, ulong_ptr W,
                             ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB,
                             nmod_t mod)
{
    const ulong invL = n_invmod(B[lenB - 1], mod.n);
    slong iR;
    ulong_ptr ptrQ = Q - lenB + 1;
    ulong_ptr R1 = W;
    
    FLINT_MPN_COPYI(R1, A, lenA);

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
                const ulong c = n_negmod(ptrQ[iR], mod.n);
                mpn_addmul_1(R1 + iR - lenB + 1, B, lenB - 1, c);
            }
        }
    }

    if (lenB > 1)
        _nmod_vec_reduce(R, R1, lenB - 1, mod);
}

void
_nmod_poly_divrem_basecase_2(ulong_ptr Q, ulong_ptr R, ulong_ptr W,
                             ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB,
                             nmod_t mod)
{
    const ulong invL = n_invmod(B[lenB - 1], mod.n);
    slong iR, i;
    ulong_ptr B2 = W, R2 = W + 2*(lenB - 1), ptrQ = Q - lenB + 1;

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
        ulong r = 
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
                const ulong c = n_negmod(ptrQ[iR], mod.n);
                mpn_addmul_1(R2 + 2 * (iR - lenB + 1), B2, 2 * lenB - 2, c);
            }
            iR--;
        }
    }

    for (iR = 0; iR < lenB - 1; iR++)
        R[iR] = n_ll_mod_preinv(R2[2*iR+1], R2[2*iR], mod.n, mod.ninv);
}

void
_nmod_poly_divrem_basecase_3(ulong_ptr Q, ulong_ptr R, ulong_ptr W,
                             ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB,
                             nmod_t mod)
{
    const ulong invL = n_invmod(B[lenB - 1], mod.n);
    slong iR, i;
    ulong_ptr B3 = W, R3 = W + 3*(lenB - 1), ptrQ = Q - lenB + 1;

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
        ulong r =
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
                const ulong c = n_negmod(ptrQ[iR], mod.n);
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
_nmod_poly_divrem_basecase(ulong_ptr Q, ulong_ptr R, ulong_ptr W,
                           ulong_srcptr A, slong lenA, ulong_srcptr B, slong lenB,
                           nmod_t mod)
{
    const slong bits =
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
    const slong lenA = A->length, lenB = B->length;
    ulong_ptr Q_coeffs, R_coeffs, W;
    nmod_poly_t t1, t2;
    TMP_INIT;

    if (lenB == 0)
    {
        if (nmod_poly_modulus(B) == 1)
        {
            nmod_poly_set(Q, A);
            nmod_poly_zero(R);
            return;
        }
        else
            flint_throw(FLINT_DIVZERO, "nmod_poly_divrem\n");
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

    TMP_START;
    W = TMP_ALLOC(NMOD_DIVREM_BC_ITCH(lenA, lenB, A->mod)*sizeof(ulong));
    
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

    TMP_END;
    _nmod_poly_normalise(R);
}
