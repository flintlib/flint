/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 William Hart

******************************************************************************/

#include "nf_elem.h"

void _nf_elem_reduce(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      /* nothing to be done */
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      fmpz * const aden = QNF_ELEM_DENREF(a);
      
      if (!fmpz_is_zero(anum + 2))
      {
         fmpz * pnum = fmpq_poly_numref(nf->pol);
           
         if (nf->flag & NF_MONIC)
         {
            fmpz_submul(anum + 1, anum + 2, pnum + 1);
            fmpz_submul(anum, anum + 2, pnum);
         } else
         {
            fmpz * prod = _fmpz_vec_init(3);

            _fmpq_poly_scalar_mul_fmpq(prod, prod + 2, 
                                  pnum, pnum + 2, 2, anum + 2, aden);
            _fmpq_poly_sub_can(anum, aden, anum, aden, 2, prod, prod + 2, 2, 0);

            _fmpz_vec_clear(prod, 3);
         }

         fmpz_zero(anum + 2);
      }
   } else /* generic nf_elem */
   {
      const slong len = nf->pol->length;
      slong plen = NF_ELEM(a)->length;
      
      if (plen >= len)
      {
         if (nf->flag & NF_MONIC)
         {
            if (len <= NF_POWERS_CUTOFF)
            {
               _fmpz_poly_rem_powers_precomp(NF_ELEM_NUMREF(a), plen,
                  fmpq_poly_numref(nf->pol), len, nf->powers.zz->powers);

               _fmpq_poly_set_length(NF_ELEM(a), len - 1);
               _fmpq_poly_normalise(NF_ELEM(a));
              
            } else
            {
               fmpz * q = _fmpz_vec_init(plen - len + 1);
               fmpz * r = _fmpz_vec_init(plen);
               slong i;
               
               _fmpz_vec_set(r, NF_ELEM_NUMREF(a), plen);

               _fmpz_poly_divrem(q, NF_ELEM_NUMREF(a), r, plen, 
                  fmpq_poly_numref(nf->pol), len, 0);

               _fmpz_vec_clear(r, plen);
               _fmpz_vec_clear(q, plen - len + 1);
          
               for (i = len - 2; i >= 0 && fmpz_is_zero(NF_ELEM_NUMREF(a) + i); i--);
               NF_ELEM(a)->length = i + 1;      
            }
         }
         else
         {
            fmpq_poly_t t;
        
            if (len <= NF_POWERS_CUTOFF)
            {
               _fmpq_poly_rem_powers_precomp(NF_ELEM_NUMREF(a), 
                  fmpq_poly_denref(NF_ELEM(a)), plen,
                  fmpq_poly_numref(nf->pol), fmpq_poly_denref(nf->pol), 
                  len, nf->powers.qq->powers);

               _fmpq_poly_set_length(NF_ELEM(a), len - 1);
               _fmpq_poly_normalise(NF_ELEM(a));
            } else
            {
               fmpq_poly_init2(t, 2*len - 3);

               _fmpq_poly_rem(t->coeffs, t->den,
                  NF_ELEM(a)->coeffs, NF_ELEM(a)->den, plen, 
                  nf->pol->coeffs, nf->pol->den, len, nf->pinv.qq); 
           
               _fmpq_poly_set_length(t, len - 1);
               _fmpq_poly_normalise(t);
        
               fmpq_poly_swap(t, NF_ELEM(a));
               fmpq_poly_clear(t);
            }
         }
      }
   }
}

void nf_elem_reduce(nf_elem_t a, const nf_t nf)
{
   if (!(nf->flag & NF_LINEAR))
      _nf_elem_reduce(a, nf);
   
   nf_elem_canonicalise(a, nf);
}
