/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_mod_poly.h"
#include "fmpz_vec.h"

#define __set(B, lenB, A, lenA)      \
do {                                 \
    _fmpz_vec_set((B), (A), (lenA)); \
    (lenB) = (lenA);                 \
} while (0)

#define __rem(R, lenR, A, lenA, B, lenB)                              \
do {                                                                  \
    if ((lenA) >= (lenB))                                             \
    {                                                                 \
        fmpz_invmod(invB, B + lenB - 1, mod);                         \
        _fmpz_mod_poly_rem((R), (A), (lenA), (B), (lenB), invB, mod); \
        (lenR) = (lenB) - 1;                                          \
        FMPZ_VEC_NORM((R), (lenR));                                   \
    }                                                                 \
    else                                                              \
    {                                                                 \
        _fmpz_vec_set((R), (A), (lenA));                              \
        (lenR) = (lenA);                                              \
    }                                                                 \
} while (0)

/*
    XXX: Currently supports aliasing between {A,a} and {B,b}.
 */

slong _fmpz_mod_poly_hgcd_res(fmpz **M, slong *lenM, 
                     fmpz *A, slong *lenA, fmpz *B, slong *lenB, 
                     const fmpz *a, slong lena, const fmpz *b, slong lenb, 
                     const fmpz_t mod, fmpz_t r)
{
    const slong lenW = 22 * lena + 16 * (FLINT_CLOG2(lena) + 1);
    slong sgnM;
    fmpz_mod_poly_res_t res;
    fmpz *W;
    
    fmpz_init(res->res);
    fmpz_init(res->lc);
    
    fmpz_set(res->res, r);
    fmpz_set(res->lc, b + lenb - 1);
    res->len0 = lena;
    res->len1 = lenb;
    res->off = 0;

    W = _fmpz_vec_init(lenW);

    if (M == NULL)
    {
        sgnM = _fmpz_mod_poly_hgcd_recursive(NULL, NULL, 
                                         A, lenA, B, lenB, 
                                         a, lena, b, lenb, W, mod, 0, res);
    }
    else
    {
        sgnM = _fmpz_mod_poly_hgcd_recursive(M, lenM, 
                                         A, lenA, B, lenB, 
                                         a, lena, b, lenb, W, mod, 1, res);
    }
    
    if (*lenB < lenb) /* make sure something happened */
    {
       if (*lenB >= 1)
       {
          fmpz_powm_ui(res->lc, res->lc, res->len0 - *lenB, mod);
          fmpz_mul(res->res, res->res, res->lc);
          fmpz_mod(res->res, res->res, mod);
          
          if (((res->len0 | res->len1) & 1) == 0)
             fmpz_negmod(res->res, res->res, mod);
       } else
       {
          if (res->len1 == 1) 
          {
             fmpz_powm_ui(res->lc, res->lc, res->len0 - 1, mod);
             fmpz_mul(res->res, res->res, res->lc);
             fmpz_mod(res->res, res->res, mod);
          } else
             fmpz_zero(res->res);
       }
    }

    fmpz_set(r, res->res);

    fmpz_clear(res->res);
    fmpz_clear(res->lc);
    _fmpz_vec_clear(W, lenW);

    return sgnM;
}

/*
    XXX: Incidentally, this implementation currently supports aliasing.  
    But since this may change in the future, no function other than 
    nmod_poly_resultant_hgcd() should rely on this.
 */

void _fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz *A, slong lenA, 
                               const fmpz *B, slong lenB, const fmpz_t mod)
{
    slong len1 = FLINT_MIN(lenA, lenB), len2 = 2 * lenB;
    fmpz *G = _fmpz_vec_init(len1);
    fmpz *J = _fmpz_vec_init(len2);
    fmpz *R = J + lenB;
    fmpz_t lc, invB;

    slong lenG, lenJ, lenR;

    fmpz_init(invB);
    fmpz_init(lc);

    fmpz_set_ui(res, 1);
    fmpz_set(lc, B + lenB - 1);

    __rem(R, lenR, A, lenA, B, lenB);

    if (lenR == 0)
    {
        if (lenB == 1)
        {
           fmpz_powm_ui(lc, lc, lenA - 1, mod);
           fmpz_mul(res, res, lc);
           fmpz_mod(res, res, mod);
        } else
           fmpz_zero(res);
    }
    else
    {
       fmpz_powm_ui(lc, lc, lenA - lenR, mod);
       fmpz_mul(res, res, lc);
       fmpz_mod(res, res, mod);
       
       if (((lenA | lenB) & 1) == 0)
          fmpz_negmod(res, res, mod);
       
       _fmpz_mod_poly_hgcd_res(NULL, NULL, G, &(lenG), J, &(lenJ), B, lenB, R, lenR, mod, res);

        while (lenJ != 0)
        {
            fmpz_set(lc, J + lenJ - 1);
            
            __rem(R, lenR, G, lenG, J, lenJ);

            if (lenR == 0)
            {
               if (lenJ == 1)
               {
                  fmpz_powm_ui(lc, lc, lenG - 1, mod);
                  fmpz_mul(res, res, lc);
                  fmpz_mod(res, res, mod);
               } else
                  fmpz_zero(res);
              
               break;
            } else
            {
               fmpz_powm_ui(lc, lc, lenG - lenR, mod);
               fmpz_mul(res, res, lc);
               fmpz_mod(res, res, mod);

               if (((lenG | lenJ) & 1) == 0)
                  fmpz_negmod(res, res, mod);
            }

            if (lenJ < FMPZ_MOD_POLY_GCD_CUTOFF)
            {
                fmpz_t r;
                
                fmpz_init(r);
                _fmpz_mod_poly_resultant_euclidean(r, J, lenJ, R, lenR, mod);

                fmpz_mul(res, res, r);
                fmpz_mod(res, res, mod);
                
                fmpz_clear(r);

                break;
            }

            _fmpz_mod_poly_hgcd_res(NULL, NULL, G, &(lenG), J, &(lenJ), J, lenJ, R, lenR, mod, res);
        }
    }

    _fmpz_vec_clear(J, len2);
    _fmpz_vec_clear(G, len1);

    fmpz_clear(lc);
    fmpz_clear(invB);
}

void fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz_mod_poly_t A, 
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    if (A->length == 0 || B->length == 0) 
    {
       fmpz_zero(res);
    }
    else if (A->length < B->length)
    {
        fmpz_mod_poly_resultant_hgcd(res, B, A, ctx);

        if (((A->length | B->length) & 1) == 0)
           fmpz_negmod(res, res, fmpz_mod_ctx_modulus(ctx));
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length;
        
        /* lenA >= lenB >= 1 */
        _fmpz_mod_poly_resultant_hgcd(res, A->coeffs, lenA,
                                   B->coeffs, lenB, fmpz_mod_ctx_modulus(ctx));
    }
}

#undef __set
#undef __rem

