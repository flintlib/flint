/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nf_elem.h"

void nf_elem_get_coeff_fmpq(fmpq_t a, const nf_elem_t b,
                                                        slong i, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      if (i > 0)
         fmpq_zero(a);
      else
      {
         fmpz_set(fmpq_numref(a), LNF_ELEM_NUMREF(b));
         fmpz_set(fmpq_denref(a), LNF_ELEM_DENREF(b));
      }
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);

      if (i > 2) /* element may be unreduced */
         fmpq_zero(a);
      else
      {
         fmpz_set(fmpq_numref(a), bnum + i);
         fmpz_set(fmpq_denref(a), QNF_ELEM_DENREF(b));
      }

      fmpq_canonicalise(a);
   } else
      fmpq_poly_get_coeff_fmpq(a, NF_ELEM(b), i);
}
