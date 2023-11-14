/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "nf_elem.h"

void _nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      if (a == b)
         fmpz_swap(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a));
      else
      {
         fmpz_set(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(b));
         fmpz_set(LNF_ELEM_DENREF(a), LNF_ELEM_NUMREF(b));
      }
      _fmpq_canonicalise(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      fmpz * const aden = QNF_ELEM_DENREF(a);
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      const fmpz * const bden = QNF_ELEM_DENREF(b);
      fmpz * t = _fmpz_vec_init(6);
      slong len = 2;

      while (len > 0 && fmpz_is_zero(bnum + len - 1))
         len--;

      _fmpq_poly_xgcd(t + 3, t + 5, t, t + 2, anum, aden,
             fmpq_poly_numref(nf->pol), fmpq_poly_denref(nf->pol), 3, bnum, bden, len);

      _fmpz_vec_clear(t, 6);
   } else
   {
      fmpq_poly_t g, t;

      fmpq_poly_init(g);
      fmpq_poly_init(t);

      fmpq_poly_xgcd(g, NF_ELEM(a), t, NF_ELEM(b), nf->pol);

      fmpq_poly_clear(t);
      fmpq_poly_clear(g);
   }
}

void nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   nf_elem_t t;

   if (a == b)
   {
      nf_elem_init(t, nf);

      _nf_elem_inv(t, b, nf);
      nf_elem_swap(t, a, nf);

      nf_elem_clear(t, nf);
   }
   else
      _nf_elem_inv(a, b, nf);
}
