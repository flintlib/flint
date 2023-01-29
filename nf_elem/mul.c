/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 William Hart

******************************************************************************/

#include "nf_elem.h"

void _nf_elem_mul_gaussian(fmpz * anum, fmpz * aden,
        const fmpz * bnum, const fmpz * bden,
        const fmpz * cnum, const fmpz * cden)
{
   fmpz_t t;
   fmpz_init(t);

   if (anum == bnum || anum == cnum)  /* aliasing */
   {
      if (bnum == cnum && bden == cden) /* squaring */
      {
         fmpz_fmms(t, bnum + 0, bnum + 0, bnum + 1, bnum + 1);
         fmpz_mul(anum + 1, bnum + 0, bnum + 1);
         fmpz_mul_2exp(anum + 1, anum + 1, 1);
      }
      else
      {
         fmpz_fmms(t, bnum + 0, cnum + 0, bnum + 1, cnum + 1);
         fmpz_fmma(anum + 1, bnum + 0, cnum + 1, bnum + 1, cnum + 0);
      }
      fmpz_swap(anum + 0, t);
   }
   else
   {
      if (bnum == cnum && bden == cden) /* squaring */
      {
         fmpz_fmms(anum + 0, bnum + 0, bnum + 0, bnum + 1, bnum + 1);
         fmpz_mul(anum + 1, bnum + 0, bnum + 1);
         fmpz_mul_2exp(anum + 1, anum + 1, 1);
      }
      else
      {
         fmpz_fmms(anum + 0, bnum + 0, cnum + 0, bnum + 1, cnum + 1);
         fmpz_fmma(anum + 1, bnum + 0, cnum + 1, bnum + 1, cnum + 0);
      }
   }
   fmpz_zero(anum + 2);
   fmpz_mul(aden, bden, cden);
   if (!fmpz_is_one(aden))
   {
#if (                                                           \
        (__FLINT_VERSION > 2)                                   \
        ||                                                      \
        (__FLINT_VERSION >= 2 && __FLINT_VERSION_MINOR >= 8)    \
    )
      fmpz_gcd3(t, anum + 0, anum + 1, aden);
#else
      fmpz_gcd(t, anum + 0, anum + 1);
      fmpz_gcd(t, t, aden);
#endif
      if (!fmpz_is_one(t))
      {
         fmpz_divexact(anum + 0, anum + 0, t);
         fmpz_divexact(anum + 1, anum + 1, t);
         fmpz_divexact(aden, aden, t);
      }
   }
   fmpz_clear(t);
}

void _nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, 
                                     const nf_elem_t c, const nf_t nf, int red)
{
   if (nf->flag & NF_LINEAR)
   {
      const fmpz * const bnum = LNF_ELEM_NUMREF(b);
      const fmpz * const bden = LNF_ELEM_DENREF(b);
      const fmpz * const cnum = LNF_ELEM_NUMREF(c);
      const fmpz * const cden = LNF_ELEM_DENREF(c);
      fmpz * const anum = LNF_ELEM_NUMREF(a);
      fmpz * const aden = LNF_ELEM_DENREF(a);
      
      fmpz_mul(anum, bnum, cnum);
      fmpz_mul(aden, bden, cden);
   }
   else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      const fmpz * const bden = QNF_ELEM_DENREF(b);
      const fmpz * const cnum = QNF_ELEM_NUMREF(c);
      const fmpz * const cden = QNF_ELEM_DENREF(c);
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      fmpz * const aden = QNF_ELEM_DENREF(a);

      fmpz_mul(anum, bnum, cnum);
      fmpz_fmma(anum + 1, bnum, cnum + 1, bnum + 1, cnum);
      fmpz_mul(anum + 2, bnum + 1, cnum + 1);

      fmpz_mul(aden, bden, cden);

      if (red && !fmpz_is_zero(anum + 2))
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
      const slong len1 = NF_ELEM(b)->length;
      const slong len2 = NF_ELEM(c)->length;
      const slong len = nf->pol->length;
      slong plen = len1 + len2 - 1;
      
      if (len1 == 0 || len2 == 0)
      {
         nf_elem_zero(a, nf);

         return;
      }

      fmpq_poly_fit_length(NF_ELEM(a), plen);
      if (len1 >= len2)
      {
         _fmpz_poly_mul(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(b), len1,
            NF_ELEM_NUMREF(c), len2);
      }
      else
      {
          _fmpz_poly_mul(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(c), len2,
             NF_ELEM_NUMREF(b), len1);
      }

      fmpz_mul(fmpq_poly_denref(NF_ELEM(a)), fmpq_poly_denref(NF_ELEM(b)),
         fmpq_poly_denref(NF_ELEM(c)));

      _fmpq_poly_set_length(NF_ELEM(a), plen);

      if (red && plen >= len)
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

void _nf_elem_mul(nf_elem_t a, const nf_elem_t b, 
                                     const nf_elem_t c, const nf_t nf)
{
   _nf_elem_mul_red(a, b, c, nf, 1);
}

void nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, 
                                     const nf_elem_t c, const nf_t nf, int red)
{
   nf_elem_t t;
   
   if (nf->flag & NF_LINEAR)
   {
      _fmpq_mul(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a), 
                LNF_ELEM_NUMREF(b), LNF_ELEM_DENREF(b),
                LNF_ELEM_NUMREF(c), LNF_ELEM_DENREF(c));
   }
   else if ((nf->flag & NF_GAUSSIAN) &&
      fmpz_is_zero(QNF_ELEM_NUMREF(b) + 2) &&
      fmpz_is_zero(QNF_ELEM_NUMREF(c) + 2))
   {
     _nf_elem_mul_gaussian(QNF_ELEM_NUMREF(a), QNF_ELEM_DENREF(a),
        QNF_ELEM_NUMREF(b), QNF_ELEM_DENREF(b),
        QNF_ELEM_NUMREF(c), QNF_ELEM_DENREF(c));
   }
   else
   {
      if (a == b || a == c)
      {
         nf_elem_init(t, nf);

         _nf_elem_mul_red(t, b, c, nf, red);
         nf_elem_swap(t, a, nf);

         nf_elem_clear(t, nf);
      }
      else
         _nf_elem_mul_red(a, b, c, nf, red);

      nf_elem_canonicalise(a, nf);
   }
}

void nf_elem_mul(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
   nf_elem_mul_red(a, b, c, nf, 1);
}
