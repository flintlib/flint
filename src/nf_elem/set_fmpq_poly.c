/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nf_elem.h"

void nf_elem_set_fmpq_poly(nf_elem_t a, const fmpq_poly_t pol, const nf_t nf)
{
   if (fmpq_poly_length(pol) >= fmpq_poly_length(nf->pol))
   {
       fmpq_poly_t r;
       fmpq_poly_init(r);
       fmpq_poly_rem(r, pol, nf->pol);
       nf_elem_set_fmpq_poly(a, r, nf);
       fmpq_poly_clear(r);
       return;
   }

   if (nf->flag & NF_LINEAR)
   {
      if (pol->length == 0)
	  {
	     fmpz_zero(LNF_ELEM_NUMREF(a));
		 fmpz_one(LNF_ELEM_DENREF(a));
	  } else if (pol->length == 1)
	  {
	     fmpz_set(LNF_ELEM_NUMREF(a), fmpq_poly_numref(pol));
         fmpz_set(LNF_ELEM_DENREF(a), fmpq_poly_denref(pol));
	  }
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);

      if (pol->length == 0)
	  {
	     fmpz_zero(anum);
		 fmpz_zero(anum + 1);
		 fmpz_one(QNF_ELEM_DENREF(a));
	  } else if (pol->length == 1)
	  {
		 fmpz_zero(anum + 1);
		 fmpz_set(anum, fmpq_poly_numref(pol));
		 fmpz_set(QNF_ELEM_DENREF(a), fmpq_poly_denref(pol));
	  } else
	  {
	     fmpz_set(anum, fmpq_poly_numref(pol));
	     fmpz_set(anum + 1, fmpq_poly_numref(pol) + 1);
         fmpz_set(QNF_ELEM_DENREF(a), fmpq_poly_denref(pol));
	  }
   } else
      fmpq_poly_set(NF_ELEM(a), pol);
}
