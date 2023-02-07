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

void nf_elem_get_coeff_fmpz(fmpz_t a, const nf_elem_t b, 
                                                        slong i, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      if (i > 0)
         fmpz_zero(a);
      else
         fmpz_set(a, LNF_ELEM_NUMREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      
      if (i > 2) /* element may be unreduced */
         fmpz_zero(a);
      else
         fmpz_set(a, bnum + i);
   } else
      fmpq_poly_get_coeff_fmpz(a, NF_ELEM(b), i);
}
