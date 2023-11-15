/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "nf_elem.h"

void _nf_elem_trace(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf)
{
   slong i;

   if (nf->flag & NF_LINEAR)
   {
      const fmpz * const anum = LNF_ELEM_NUMREF(a);
      const fmpz * const aden = LNF_ELEM_DENREF(a);

      fmpz_set(rnum, anum);
      fmpz_set(rden, aden);
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const anum = QNF_ELEM_NUMREF(a);
      const fmpz * const aden = QNF_ELEM_DENREF(a);
      const fmpz * const tnum = fmpq_poly_numref(nf->traces);
      const fmpz * const tden = fmpq_poly_denref(nf->traces);

      slong alen = 2;
      while (alen > 0 && fmpz_is_zero(anum + alen - 1))
         alen--;

      if (alen == 0)
      {
         fmpz_zero(rnum);
         fmpz_one(rden);
      } else
      {
         fmpz_mul(rnum, anum, tnum);

         if (alen == 2)
            fmpz_addmul(rnum, anum + 1, tnum + 1);

         fmpz_mul(rden, aden, tden);

         _fmpq_canonicalise(rnum, rden);
      }
   } else /* generic nf_elem */
   {
      const fmpz * const anum = NF_ELEM_NUMREF(a);
      const fmpz * const aden = NF_ELEM_DENREF(a);
      const fmpz * const tnum = fmpq_poly_numref(nf->traces);
      const fmpz * const tden = fmpq_poly_denref(nf->traces);

      slong alen = NF_ELEM(a)->length;

      if (alen == 0)
      {
         fmpz_zero(rnum);
         fmpz_one(rden);
      } else
      {
         fmpz_mul(rnum, anum, tnum);

         for (i = 1; i < alen; i++)
            fmpz_addmul(rnum, anum + i, tnum + i);

         fmpz_mul(rden, aden, tden);

         _fmpq_canonicalise(rnum, rden);
      }
   }
}

void nf_elem_trace(fmpq_t res, const nf_elem_t a, const nf_t nf)
{
   _nf_elem_trace(fmpq_numref(res), fmpq_denref(res), a, nf);
}
