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

void nf_elem_fmpz_sub(nf_elem_t a, const fmpz_t c, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz * den = LNF_ELEM_DENREF(a);
	  fmpz * num = LNF_ELEM_NUMREF(a);
	  const fmpz * const den2 = LNF_ELEM_DENREF(b);
	  const fmpz * const num2 = LNF_ELEM_NUMREF(b);
	  
      _fmpq_sub_fmpz(num, den, num2, den2, c);
	  fmpz_neg(num, num);
   }
   else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * den = QNF_ELEM_DENREF(a);
      fmpz * num = QNF_ELEM_NUMREF(a);
	  
      nf_elem_neg(a, b, nf);
      fmpz_addmul(num, den, c);
      _fmpq_poly_canonicalise(num, den, 2);
   } else
   {
      fmpq_poly_fmpz_sub(NF_ELEM(a), c, NF_ELEM(b));
   }
}
