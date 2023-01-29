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

char * nf_elem_get_str_pretty(const nf_elem_t a, 
                              const char * var, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      const fmpz * const den = LNF_ELEM_DENREF(a);
      const fmpz * const num = LNF_ELEM_NUMREF(a);
      slong len = 1 - fmpz_is_zero(num);

      return _fmpq_poly_get_str_pretty(num, den, len, var);
   }
   else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const den = QNF_ELEM_DENREF(a);
      const fmpz * const num = QNF_ELEM_NUMREF(a);
      slong len = 3;

      while (len != 0 && fmpz_is_zero(num + len - 1))
         len--;

      return _fmpq_poly_get_str_pretty(num, den, len, var);
   } else
   {
      return fmpq_poly_get_str_pretty(NF_ELEM(a), var);
   }
}
