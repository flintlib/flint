/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart
    Copyright (C) 2015 Claus Fieker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "nf_elem.h"
#include "fmpq_poly.h"

void nf_elem_set_fmpz_mat_row(nf_elem_t b, const fmpz_mat_t M,
                                      const slong i, fmpz_t den, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(LNF_ELEM_NUMREF(b), fmpz_mat_entry(M, i, 0));
      fmpz_set(LNF_ELEM_DENREF(b), den);
      _fmpq_canonicalise(LNF_ELEM_NUMREF(b), LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz_t d;
      fmpz * const bnum = QNF_ELEM_NUMREF(b);

      fmpz_init(d);

      fmpz_set(bnum, fmpz_mat_entry(M, i, 0));
      fmpz_set(bnum + 1, fmpz_mat_entry(M, i, 1));
      fmpz_set(QNF_ELEM_DENREF(b), den);

      fmpz_gcd(d, bnum, bnum + 1);

      if (!fmpz_is_one(d))
      {
          fmpz_gcd(d, d, QNF_ELEM_DENREF(b));

          if (!fmpz_is_one(d))
          {
              fmpz_divexact(bnum, bnum, d);
              fmpz_divexact(bnum + 1, bnum + 1, d);
              fmpz_divexact(QNF_ELEM_DENREF(b), QNF_ELEM_DENREF(b), d);
          }
       }

       fmpz_clear(d);
   } else
   {
      slong j;
      for (j = nf->pol->length - 2; j >= 0; j--)
         if (!fmpz_is_zero(fmpz_mat_entry(M, i, j)))
            break;
      _fmpq_poly_set_length(NF_ELEM(b), j + 1);
      for (; j >= 0; j--)
      fmpq_poly_set_coeff_fmpz(NF_ELEM(b), j, fmpz_mat_entry(M, i, j));
      fmpz_set(NF_ELEM_DENREF(b), den);
      fmpq_poly_canonicalise(NF_ELEM(b));
   }
}
