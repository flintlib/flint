/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart

******************************************************************************/

#include "nf.h"

void nf_init(nf_t nf, const fmpq_poly_t pol)
{
    slong i, j;
    slong len = pol->length, deg = len - 1;

    fmpz * pow, * lead = fmpq_poly_numref(pol) + len - 1;

    fmpq_poly_init(nf->pol);
    fmpq_poly_set(nf->pol, pol);

    /**** Set up precomputed inverse of leading coeff of f(x) ****/

    if (fmpz_is_one(fmpq_poly_denref(pol)) /* denominator is one and numerator is monic */
     && fmpz_is_one(lead))
       nf->flag = NF_MONIC;
    else
    {
       fmpz_preinvn_init(nf->pinv.qq, lead);
       nf->flag = NF_GENERIC;
    }

    /**** Set up precomputed powers x^i mod f(x) ****/

    if (len < 2)
    {
       flint_printf("Exception (nf_init). Degree must be at least 1.\n");
       abort();
    } else if (len == 2) /* linear case */
       nf->flag |= NF_LINEAR;
    else if (len == 3) /* quadratic case */
    {
       nf->flag |= NF_QUADRATIC;
       if (fmpz_is_one(pol->coeffs + 0) && fmpz_is_zero(pol->coeffs + 1) &&
            fmpz_is_one(pol->coeffs + 2) && fmpz_is_one(pol->den))
          nf->flag |= NF_GAUSSIAN;
    }
    else if (len <= NF_POWERS_CUTOFF) /* compute powers of generator mod pol */
    {
       if (nf->flag & NF_MONIC)
       {
          nf->powers.zz->powers = _fmpz_poly_powers_precompute(fmpq_poly_numref(pol), 
                                       len);
          nf->powers.zz->len = len;
       }
       else
       {
          nf->powers.qq->powers = _fmpq_poly_powers_precompute(fmpq_poly_numref(pol), 
                                       fmpq_poly_denref(pol), len);
          nf->powers.qq->len = len;
       }
   }

   /**** Set up precomputed traces S_k = \sum _i theta_i^k for roots theta_i of f(x) ****/

   /* 
      Uses the recursive formula from pp. 163 of "A Course in Computational Algebraic 
      Number Theory" by Henri Cohen
   */

  fmpq_poly_init2(nf->traces, deg);
   pow = fmpq_poly_denref(nf->traces);

   for (i = 1; i < deg; i++)
   {
      fmpz_mul_si(fmpq_poly_numref(nf->traces) + i, 
         fmpq_poly_numref(pol) + deg - i, i); 
      
      for (j = i - 1; j >= 1; j--)
      {
         fmpz_mul(fmpq_poly_numref(nf->traces) + i, 
            fmpq_poly_numref(nf->traces) + i, lead);
         fmpz_addmul(fmpq_poly_numref(nf->traces) + i, 
            fmpq_poly_numref(pol) + deg - j, 
            fmpq_poly_numref(nf->traces) + i - j);
      }
      
      fmpz_neg(fmpq_poly_numref(nf->traces) + i, 
         fmpq_poly_numref(nf->traces) + i);
   }

   for (i = 1; i < deg; i++)
   {
      fmpz_mul(fmpq_poly_numref(nf->traces) + deg - i,
         fmpq_poly_numref(nf->traces) + deg - i, pow);
      fmpz_mul(pow, pow, lead);
   }

   fmpz_mul_si(fmpq_poly_numref(nf->traces), pow, deg);
}

