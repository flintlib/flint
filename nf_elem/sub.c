/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart

******************************************************************************/

#include "nf_elem.h"

void _nf_elem_sub_lf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can)
{
   const fmpz * const p = LNF_ELEM_NUMREF(b);
   const fmpz * const q = LNF_ELEM_DENREF(b);
   const fmpz * const r = LNF_ELEM_NUMREF(c);
   const fmpz * const s = LNF_ELEM_DENREF(c);
   fmpz * const rnum = LNF_ELEM_NUMREF(a);
   fmpz * const rden = LNF_ELEM_DENREF(a);
   fmpz_t t;

   if (can)
      _fmpq_sub(rnum, rden, p, q, r, s);
   else
   {
      /* Same denominator */
      if (fmpz_equal(q, s))
      {
         fmpz_sub(rnum, p, r);
         fmpz_set(rden, q);

         return;
      }

      /* p/q is an integer */
      if (fmpz_is_one(q))
      {
         fmpz_init(t);

         fmpz_mul(t, p, s);
         fmpz_sub(rnum, t, r);
         fmpz_set(rden, s);

         fmpz_clear(t);

         return;
      }

      /* r/s is an integer */
      if (fmpz_is_one(s))
      {
         fmpz_init(t);

         fmpz_mul(t, r, q);
         fmpz_sub(rnum, t, p);
         fmpz_set(rden, q);

         fmpz_clear(t);

         return;
      }

      /*
         We want to compute p/q - r/s which is (p*s - q*r, q*s).
      */

      fmpz_init(t);

      fmpz_mul(t, q, r);
      fmpz_mul(rnum, p, s);
      fmpz_sub(rnum, rnum, t);
      fmpz_mul(rden, q, s);

      fmpz_clear(t);
   }
}

void _nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can)
{
   fmpz_t d;

   const fmpz * const bnum = QNF_ELEM_NUMREF(b);
   const fmpz * const bden = QNF_ELEM_DENREF(b);
   
   const fmpz * const cnum = QNF_ELEM_NUMREF(c);
   const fmpz * const cden = QNF_ELEM_DENREF(c);
   
   fmpz * const anum = QNF_ELEM_NUMREF(a);
   fmpz * const aden = QNF_ELEM_DENREF(a);

   fmpz_init(d);
   fmpz_one(d);

   if (fmpz_equal(bden, cden))
   {
      fmpz_sub(anum, bnum, cnum);
      fmpz_sub(anum + 1, bnum + 1, cnum + 1);
      fmpz_sub(anum + 2, bnum + 2, cnum + 2);
      fmpz_set(aden, bden);

      if (can && !fmpz_is_one(aden))
      {
#if (                                                           \
        (__FLINT_VERSION > 2)                                   \
        ||                                                      \
        (__FLINT_VERSION >= 2 && __FLINT_VERSION_MINOR >= 8)    \
    )
         fmpz_gcd3(d, anum + 0, anum + 1, anum + 2);
#else
         fmpz_gcd(d, anum + 0, anum + 1);
         fmpz_gcd(d, d, anum + 2);
#endif
         if (!fmpz_is_one(d))
         {
            fmpz_gcd(d, d, aden);

            if (!fmpz_is_one(d))
            {
               fmpz_divexact(anum, anum, d);
               fmpz_divexact(anum + 1, anum + 1, d);
               fmpz_divexact(anum + 2, anum + 2, d);
               fmpz_divexact(aden, aden, d);
            }
         }
      }

      fmpz_clear(d);

      return;
   }

   if (!fmpz_is_one(bden) && !fmpz_is_one(cden))
      fmpz_gcd(d, bden, cden);

   if (fmpz_is_one(d))
   {
      fmpz_mul(anum, bnum, cden);
      fmpz_mul(anum + 1, bnum + 1, cden);
      fmpz_mul(anum + 2, bnum + 2, cden);
      fmpz_submul(anum, cnum, bden);
      fmpz_submul(anum + 1, cnum + 1, bden);
      fmpz_submul(anum + 2, cnum + 2, bden);
      fmpz_mul(aden, bden, cden);
   } else
   {
      fmpz_t bden1;
      fmpz_t cden1;
      
      fmpz_init(bden1);
      fmpz_init(cden1);
      
      fmpz_divexact(bden1, bden, d);
      fmpz_divexact(cden1, cden, d);
        
      fmpz_mul(anum, bnum, cden1);
      fmpz_mul(anum + 1, bnum + 1, cden1);
      fmpz_mul(anum + 2, bnum + 2, cden1);
      fmpz_submul(anum, cnum, bden1);
      fmpz_submul(anum + 1, cnum + 1, bden1);
      fmpz_submul(anum + 2, cnum + 2, bden1);
      
      if (fmpz_is_zero(anum) && fmpz_is_zero(anum + 1) && fmpz_is_zero(anum + 2))
         fmpz_one(aden);
      else
      {
         if (can)
         {
            fmpz_t e;
            
            fmpz_init(e);
              
#if (                                                           \
        (__FLINT_VERSION > 2)                                   \
        ||                                                      \
        (__FLINT_VERSION >= 2 && __FLINT_VERSION_MINOR >= 8)    \
    )
            fmpz_gcd3(e, anum + 0, anum + 1, anum + 2);
#else
            fmpz_gcd(e, anum + 0, anum + 1);
            fmpz_gcd(e, e, anum + 2);
#endif
            if (!fmpz_is_one(e))
               fmpz_gcd(e, e, d);
            
            if (fmpz_is_one(e))
               fmpz_mul(aden, bden, cden1);
            else
            {
                fmpz_divexact(anum, anum, e);
                fmpz_divexact(anum + 1, anum + 1, e);
                fmpz_divexact(anum + 2, anum + 2, e);
                fmpz_divexact(bden1, bden, e);
                fmpz_mul(aden, bden1, cden1);
            }
            
            fmpz_clear(e);
         } else
            fmpz_mul(aden, bden, cden1);
      }

      fmpz_clear(bden1);
      fmpz_clear(cden1);
   }

   fmpz_clear(d);
}

void nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
   if (a == c)
   {
      nf_elem_t t;

      nf_elem_init(t, nf);

      _nf_elem_sub_qf(t, b, c, nf, 1);
      nf_elem_swap(t, a, nf);

      nf_elem_clear(t, nf);
   } else
      _nf_elem_sub_qf(a, b, c, nf, 1);
}
