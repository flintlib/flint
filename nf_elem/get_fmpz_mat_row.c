/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart
    Copyright (C) 2015 Claus Fieker

******************************************************************************/

#include "nf_elem.h"

void nf_elem_get_fmpz_mat_row(fmpz_mat_t M, const slong i, fmpz_t den, 
                                              const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(fmpz_mat_entry(M, i, 0), LNF_ELEM_NUMREF(b));
      fmpz_set(den, LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      fmpz_set(fmpz_mat_entry(M, i, 0), bnum);
      fmpz_set(fmpz_mat_entry(M, i, 1), bnum + 1);
      fmpz_set(den, QNF_ELEM_DENREF(b));
   } else
   {
      slong j;
      for (j = 0; j < NF_ELEM(b)->length; j++) 
      {
         fmpz_set(fmpz_mat_entry(M, i, j), NF_ELEM_NUMREF(b) + j);
      }
      for ( ; j < nf->pol->length - 1; j++)
      {
         fmpz_zero(fmpz_mat_entry(M, i, j));
      }
      fmpz_set(den, NF_ELEM_DENREF(b));
   }
}
