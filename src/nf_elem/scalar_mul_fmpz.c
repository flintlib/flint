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

void nf_elem_scalar_mul_fmpz(nf_elem_t a, const nf_elem_t b, 
                                                       const fmpz_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz * den = LNF_ELEM_DENREF(a);
	  fmpz * num = LNF_ELEM_NUMREF(a);
	  const fmpz * const den2 = LNF_ELEM_DENREF(b);
	  const fmpz * const num2 = LNF_ELEM_NUMREF(b);
	  
      fmpz_mul(num, num2, c);
	  fmpz_set(den, den2);
	  _fmpq_canonicalise(num, den);
   }
   else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * den = QNF_ELEM_DENREF(a);
	  fmpz * num = QNF_ELEM_NUMREF(a);
	  const fmpz * const den2 = QNF_ELEM_DENREF(b);
	  const fmpz * const num2 = QNF_ELEM_NUMREF(b);
	  
	  _fmpz_vec_scalar_mul_fmpz(num, num2, 2, c);
	  fmpz_set(den, den2);
	  _fmpq_poly_canonicalise(num, den, 2);
   } else
   {
      fmpq_poly_scalar_mul_fmpz(NF_ELEM(a), NF_ELEM(b), c);
   }

}
