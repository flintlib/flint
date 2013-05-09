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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void _fmpq_poly_xgcd(fmpz *G, fmpz_t denG, 
                     fmpz *S, fmpz_t denS, fmpz *T, fmpz_t denT, 
                     const fmpz *A, const fmpz_t denA, len_t lenA, 
                     const fmpz *B, const fmpz_t denB, len_t lenB)
{
    int alloc = 0;
    fmpz *primA, *primB, *C, *D;
    fmpz_t cA, cB;
    len_t lenG, lenC, lenD;

    fmpz_init(cA);
    fmpz_init(cB);

    _fmpz_vec_content(cA, A, lenA);
    _fmpz_vec_content(cB, B, lenB);

    if (fmpz_is_one(cA))
    {
        if (fmpz_is_one(cB))
        {
            primA = (fmpz *) A;
            primB = (fmpz *) B;
        }
        else
        {
            alloc |= 1;
            primA = (fmpz *) A;
            primB = _fmpz_vec_init(lenB);
            _fmpz_vec_scalar_divexact_fmpz(primB, B, lenB, cB);
        }
    }
    else
    {
        if (fmpz_is_one(cA))
        {
            alloc |= 2;
            primA = _fmpz_vec_init(lenA);
            primB = (fmpz *) B;
            _fmpz_vec_scalar_divexact_fmpz(primA, A, lenA, cA);
        }
        else
        {
            alloc |= 3;
            primA = _fmpz_vec_init(lenA + lenB);
            primB = primA + lenA;
            _fmpz_vec_scalar_divexact_fmpz(primA, A, lenA, cA);
            _fmpz_vec_scalar_divexact_fmpz(primB, B, lenB, cB);
        }
    }

    _fmpz_poly_gcd(G, primA, lenA, primB, lenB);

    for (lenG = lenB - 1; !G[lenG]; lenG--) ;
    lenG++;

    if (lenG > 1)
    {
        alloc |= 4;
        lenC = lenA - lenG + 1;
        lenD = lenB - lenG + 1;
        C = _fmpz_vec_init(lenC + lenD);
        D = C + lenC;
        _fmpz_poly_div(C, primA, lenA, G, lenG);
        _fmpz_poly_div(D, primB, lenB, G, lenG);
    }
    else
    {
        lenC = lenA;
        lenD = lenB;
        C = primA;
        D = primB;
    }

    _fmpz_poly_xgcd(denG, S, T, C, lenC, D, lenD);

    if (!fmpz_is_one(denA))
        _fmpz_vec_scalar_mul_fmpz(S, S, lenD, denA);
    fmpz_mul(cA, cA, denG);
    fmpz_mul(denS, cA, G + (lenG - 1));

    if (!fmpz_is_one(denB))
        _fmpz_vec_scalar_mul_fmpz(T, T, lenC, denB);
    fmpz_mul(cB, cB, denG);
    fmpz_mul(denT, cB, G + (lenG - 1));

    _fmpz_vec_zero(S + lenD, lenB - lenD);
    _fmpz_vec_zero(T + lenC, lenA - lenC);

    _fmpq_poly_canonicalise(S, denS, lenD);
    _fmpq_poly_canonicalise(T, denT, lenC);

    fmpz_set(denG, G + (lenG - 1));

    if ((alloc & 3) == 1)
        _fmpz_vec_clear(primB, lenB);
    else if ((alloc & 3) == 2)
        _fmpz_vec_clear(primA, lenA);
    else if ((alloc & 3) == 3)
        _fmpz_vec_clear(primA, lenA + lenB);

    if ((alloc & 4))
        _fmpz_vec_clear(C, lenC + lenD);

    fmpz_clear(cA);
    fmpz_clear(cB);
}

void fmpq_poly_xgcd(fmpq_poly_t G, fmpq_poly_t S, fmpq_poly_t T, 
                    const fmpq_poly_t A, const fmpq_poly_t B)
{
    if (G == S || G == T || S == T)
    {
        printf("Exception (fmpq_poly_xgcd). Output arguments aliased.\n");
        abort();
    }

    if (A->length < B->length)
    {
       fmpq_poly_xgcd(G, T, S, B, A);
    } else
    {
       len_t lenA = A->length, lenB = B->length, lenG = lenB;

       if (lenA == 0)  /* lenA = lenB = 0 */
       {
           fmpq_poly_zero(G);
           fmpq_poly_zero(S);
           fmpq_poly_zero(T);
       }
       else if (lenB == 0)  /* lenA > lenB = 0 */
       {
           fmpq_poly_make_monic(G, A);
           fmpq_poly_zero(T);
           fmpq_poly_fit_length(S, 1);
           _fmpq_poly_set_length(S, 1);
           if (fmpz_sgn(A->coeffs + (lenA - 1)) > 0)
           {
              fmpz_set(S->coeffs, A->den);
              fmpz_set(S->den, A->coeffs + (lenA - 1));
           }
           else
           {
              fmpz_neg(S->coeffs, A->den);
              fmpz_neg(S->den, A->coeffs + (lenA - 1));
           }
           fmpq_poly_canonicalise(S);
       }
       else if (lenB == 1) /* lenA >= lenB = 1 */
       {
           fmpq_poly_set_ui(G, 1);
           fmpq_poly_zero(S);
           fmpq_poly_fit_length(T, 1);
           _fmpq_poly_set_length(T, 1);
           if (fmpz_sgn(B->coeffs) > 0)
           {
               fmpz_set(T->coeffs, B->den);
               fmpz_set(T->den, B->coeffs);
           }
           else
           {
               fmpz_neg(T->coeffs, B->den);
               fmpz_neg(T->den, B->coeffs);
           }
       }
       else /* lenA >= lenB >= 2 */
       {
           /* Aliasing */
           if (G == A || G == B)
           {
               fmpq_poly_t tG;

               fmpq_poly_init2(tG, lenG);
               fmpq_poly_xgcd(tG, S, T, A, B);
               fmpq_poly_swap(tG, G);
               fmpq_poly_clear(tG);
           } 
           else if (S == A || S == B)
           {
               fmpq_poly_t tS;

               fmpq_poly_init2(tS, lenB);
               fmpq_poly_xgcd(G, tS, T, A, B);
               fmpq_poly_swap(tS, S);
               fmpq_poly_clear(tS);
           }
           else if (T == A || T == B)
           {
               fmpq_poly_t tT;

               fmpq_poly_init2(tT, lenA);
               fmpq_poly_xgcd(G, S, tT, A, B);
               fmpq_poly_swap(tT, T);
               fmpq_poly_clear(tT);
           }
           else /* no aliasing */
           {

               fmpq_poly_fit_length(G, lenG);
               fmpq_poly_fit_length(S, lenB);
               fmpq_poly_fit_length(T, lenA);

               _fmpq_poly_xgcd(G->coeffs, G->den, S->coeffs, S->den, T->coeffs, T->den, 
                    A->coeffs, A->den, lenA, B->coeffs, B->den, lenB);
               
               _fmpq_poly_set_length(G, lenG);
               _fmpq_poly_set_length(S, lenB);
               _fmpq_poly_set_length(T, lenA);

               _fmpq_poly_normalise(G);
               _fmpq_poly_normalise(S);
               _fmpq_poly_normalise(T);
           }
       }
    }
}

