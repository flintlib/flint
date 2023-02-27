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
#include "nmod_poly.h"
#include "mpn_extras.h"

#define __set(B, lenB, A, lenA)      \
do {                                 \
    _nmod_vec_set((B), (A), (lenA)); \
    (lenB) = (lenA);                 \
} while (0)

#define __rem(R, lenR, A, lenA, B, lenB)                    \
do {                                                        \
    if ((lenA) >= (lenB))                                   \
    {                                                       \
        _nmod_poly_rem((R), (A), (lenA), (B), (lenB), mod); \
        (lenR) = (lenB) - 1;                                \
        MPN_NORM((R), (lenR));                              \
    }                                                       \
    else                                                    \
    {                                                       \
        _nmod_vec_set((R), (A), (lenA));                    \
        (lenR) = (lenA);                                    \
    }                                                       \
} while (0)

/*
    XXX: Currently supports aliasing between {A,a} and {B,b}.
 */

slong _nmod_poly_hgcd_res(mp_ptr *M, slong *lenM, 
                     mp_ptr A, slong *lenA, mp_ptr B, slong *lenB, 
                     mp_srcptr a, slong lena, mp_srcptr b, slong lenb, 
                     nmod_t mod, mp_limb_t * r)
{
    const slong lenW = 22 * lena + 16 * (FLINT_CLOG2(lena) + 1);
    slong sgnM;
    nmod_poly_res_t res;
    mp_ptr W;

    res->res = *r;
    res->lc = b[lenb - 1];
    res->len0 = lena;
    res->len1 = lenb;
    res->off = 0;

    W = _nmod_vec_init(lenW);

    if (M == NULL)
    {
        sgnM = _nmod_poly_hgcd_recursive(NULL, NULL, 
                                         A, lenA, B, lenB, 
                                         a, lena, b, lenb, W, mod, 0, res);
    }
    else
    {
        sgnM = _nmod_poly_hgcd_recursive(M, lenM, 
                                         A, lenA, B, lenB, 
                                         a, lena, b, lenb, W, mod, 1, res);
    }
    
    if (*lenB < lenb) /* make sure something happened */
    {
       if (*lenB >= 1)
       {
          res->lc  = n_powmod2_preinv(res->lc, res->len0 - *lenB, mod.n, mod.ninv);
          res->res = n_mulmod2_preinv(res->res, res->lc, mod.n, mod.ninv);

          if (((res->len0 | res->len1) & 1) == 0)
             res->res = nmod_neg(res->res, mod);
       } else
       {
          if (res->len1 == 1) 
          {
             res->lc  = n_powmod2_preinv(res->lc, res->len0 - 1, mod.n, mod.ninv);
             res->res = n_mulmod2_preinv(res->res, res->lc, mod.n, mod.ninv);
          } else
             res->res = 0;
       }
    }

    *r = res->res;

    _nmod_vec_clear(W);

    return sgnM;
}

/*
    XXX: Incidentally, this implementation currently supports aliasing.  
    But since this may change in the future, no function other than 
    nmod_poly_resultant_hgcd() should rely on this.
 */

mp_limb_t _nmod_poly_resultant_hgcd(mp_srcptr A, slong lenA, 
                               mp_srcptr B, slong lenB, nmod_t mod)
{
    const slong cutoff = FLINT_BIT_COUNT(mod.n) <= 8 ? 
                        NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;

    mp_ptr G = _nmod_vec_init(FLINT_MIN(lenA, lenB));
    mp_ptr J = _nmod_vec_init(2 * lenB);
    mp_ptr R = J + lenB;
    mp_limb_t res = 1;

    slong lenG, lenJ, lenR;

    mp_limb_t lc = B[lenB - 1];

    __rem(R, lenR, A, lenA, B, lenB);

    if (lenR == 0)
    {
        if (lenB == 1)
        {
           lc  = n_powmod2_preinv(lc, lenA - 1, mod.n, mod.ninv);
           res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);
        } else
           res = 0;
    }
    else
    {
       lc  = n_powmod2_preinv(lc, lenA - lenR, mod.n, mod.ninv);
       res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);

       if (((lenA | lenB) & 1) == 0)
          res = nmod_neg(res, mod);
       
       _nmod_poly_hgcd_res(NULL, NULL, G, &(lenG), J, &(lenJ), B, lenB, R, lenR, mod, &res);

        while (lenJ != 0)
        {
            lc = J[lenJ - 1];
            
            __rem(R, lenR, G, lenG, J, lenJ);

            if (lenR == 0)
            {
               if (lenJ == 1)
               {
                  lc  = n_powmod2_preinv(lc, lenG - 1, mod.n, mod.ninv);
                  res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);
               } else
                  res = 0;
              
               break;
            } else
            {
               lc  = n_powmod2_preinv(lc, lenG - lenR, mod.n, mod.ninv);
               res = n_mulmod2_preinv(res, lc, mod.n, mod.ninv);

               if (((lenG | lenJ) & 1) == 0)
                  res = nmod_neg(res, mod);
            }

            if (lenJ < cutoff)
            {
                mp_limb_t r = _nmod_poly_resultant_euclidean(J, lenJ, R, lenR, mod);

                res = n_mulmod2_preinv(res, r, mod.n, mod.ninv);

                break;
            }

            _nmod_poly_hgcd_res(NULL, NULL, G, &(lenG), J, &(lenJ), J, lenJ, R, lenR, mod, &res);
        }
    }

    _nmod_vec_clear(J);
    _nmod_vec_clear(G);

    return res;
}

mp_limb_t nmod_poly_resultant_hgcd(const nmod_poly_t A, const nmod_poly_t B)
{
    mp_limb_t res;

    if (A->length == 0 || B->length == 0) 
    {
       return 0;
    } 
           
    if (A->length < B->length)
    {
        res = nmod_poly_resultant_hgcd(B, A);

        if (((A->length | B->length) & 1) == 0)
           res = nmod_neg(res, A->mod);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length;
        
        /* lenA >= lenB >= 1 */
        res = _nmod_poly_resultant_hgcd(A->coeffs, lenA,
                                               B->coeffs, lenB, A->mod);

    }

    return res;
}

#undef __set
#undef __rem

