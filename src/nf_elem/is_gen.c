/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 William Hart

******************************************************************************/

#include "nf_elem.h"

int nf_elem_is_gen(const nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_t t1, t2;
	  int is_gen;
	  
	  /* fast path */
	  if (fmpz_equal(LNF_ELEM_DENREF(a), nf->pol->coeffs + 1))
	    return fmpz_cmpabs(LNF_ELEM_DENREF(a), nf->pol->coeffs) == 0
            && fmpz_sgn(LNF_ELEM_DENREF(a)) == -fmpz_sgn(nf->pol->coeffs);	
			
	  /* slow path */
	  fmpz_init(t1);
	  fmpz_init(t2);
	  
	  fmpz_mul(t1, LNF_ELEM_NUMREF(a), nf->pol->coeffs + 1);
	  fmpz_mul(t2, LNF_ELEM_DENREF(a), nf->pol->coeffs);
	  fmpz_neg(t1, t1);
	  
	  is_gen = fmpz_equal(t1, t2);
	  
	  fmpz_clear(t1);
	  fmpz_clear(t2);
	  
	  return is_gen;
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      return fmpz_equal(anum + 1, QNF_ELEM_DENREF(a)) 
	      && fmpz_is_zero(anum);
   } else
      return fmpq_poly_length(NF_ELEM(a)) == 2
	      && fmpz_equal(NF_ELEM(a)->coeffs + 1, NF_ELEM(a)->den)
		  && fmpz_is_zero(NF_ELEM(a)->coeffs);
}
