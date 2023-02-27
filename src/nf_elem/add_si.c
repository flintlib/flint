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

void nf_elem_add_si(nf_elem_t a, const nf_elem_t b, slong c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz * den = LNF_ELEM_DENREF(a);
	  fmpz * num = LNF_ELEM_NUMREF(a);
	  
      nf_elem_set(a, b, nf);
	  
	  if (c >= 0)
	     fmpz_addmul_ui(num, den, c);
	  else
	     fmpz_submul_ui(num, den, -c);
	  _fmpq_canonicalise(num, den);
   }
   else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * den = QNF_ELEM_DENREF(a);
	  fmpz * num = QNF_ELEM_NUMREF(a);
	  slong len = 2;
	  
	  nf_elem_set(a, b, nf);
	  
	  while (len != 0 && fmpz_is_zero(num + len - 1))
	     len--;
	  
      if (c >= 0)
	     fmpz_addmul_ui(num, den, c);
	  else
	     fmpz_submul_ui(num, den, -c);
	  _fmpq_poly_canonicalise(num, den, len);
   } else
   {
      fmpq_poly_add_si(NF_ELEM(a), NF_ELEM(b), c);
   }
}
