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

void nf_elem_fmpq_sub(nf_elem_t a, const fmpq_t c, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz * den = LNF_ELEM_DENREF(a);
	  fmpz * num = LNF_ELEM_NUMREF(a);
	  const fmpz * const den2 = LNF_ELEM_DENREF(b);
	  const fmpz * const num2 = LNF_ELEM_NUMREF(b);
	  
      _fmpq_sub(num, den, fmpq_numref(c), fmpq_denref(c), num2, den2);
   }
   else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * den = QNF_ELEM_DENREF(a);
	  fmpz * num = QNF_ELEM_NUMREF(a);
	  const fmpz * const den2 = QNF_ELEM_DENREF(b);
	  const fmpz * const num2 = QNF_ELEM_NUMREF(b);
	  slong len = 2;
	  
	  while (len != 0 && fmpz_is_zero(num2 + len - 1))
	     len--;
	  
      if (len == 0)
	  {
	     fmpz_set(num, fmpq_numref(c));
		 fmpz_set(den, fmpq_denref(c));
	  } else if (len == 1)
	     _fmpq_sub(num, den, fmpq_numref(c), fmpq_denref(c), num2, den2);
	  else
	  {
	     /* fast path */
		 if (fmpz_equal(fmpq_denref(c), den2))
		 {
    		fmpz_sub(num, fmpq_numref(c), num2);
            fmpz_neg(num + 1, num2 + 1);
			fmpz_set(den, den2);
		 } else /* slow path */
		 {
		    fmpz_t d1, d2, g;
			
			fmpz_init(d1);
			fmpz_init(d2);
			fmpz_init(g);
			
			nf_elem_set(a, b, nf);

	        fmpz_gcd(g, fmpq_denref(c), den);
			fmpz_divexact(d1, fmpq_denref(c), g);
			fmpz_divexact(d2, den, g);

			fmpz_mul(num + 1, num + 1, d1);
			fmpz_mul(num, num, d1);
			fmpz_mul(den, den, d1);
			
			fmpz_submul(num, d2, fmpq_numref(c));
		    fmpz_neg(num, num);
		    fmpz_neg(num + 1, num + 1);
			
			fmpz_clear(g);
			fmpz_clear(d1);
			fmpz_clear(d2);
		 }
		 
		 _fmpq_poly_canonicalise(num, den, 2);
	  }
   } else
   {
      fmpq_poly_fmpq_sub(NF_ELEM(a), c, NF_ELEM(b));
   }
}
