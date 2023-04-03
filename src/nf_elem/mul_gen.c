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

void nf_elem_mul_gen(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
  if (nf->flag & NF_LINEAR)
  {
      fmpz * den = LNF_ELEM_DENREF(a);
	    fmpz * num = LNF_ELEM_NUMREF(a);
      /* _fmpq_mul assumes a positive denominator */
      if (fmpz_sgn(fmpq_poly_numref(nf->pol) + 1) < 0)
      {
          fmpz_t t;
          fmpz_init(t);
          fmpz_neg(t, fmpq_poly_numref(nf->pol) + 1);
          _fmpq_mul(num, den, LNF_ELEM_NUMREF(b), LNF_ELEM_DENREF(b), fmpq_poly_numref(nf->pol), t);
          fmpz_clear(t);
      }
      else
      {
          _fmpq_mul(num, den, LNF_ELEM_NUMREF(b), LNF_ELEM_DENREF(b), fmpq_poly_numref(nf->pol), fmpq_poly_numref(nf->pol) + 1);
          fmpz_neg(num, num);
      }

      _fmpq_canonicalise(num, den);
  }
  else if (nf->flag & NF_QUADRATIC)
  {
      fmpz * anum = QNF_ELEM_NUMREF(a);
      fmpz const * bnum = QNF_ELEM_NUMREF(b);

      fmpz_set(anum + 2, bnum + 1);
      fmpz_set(anum + 1, bnum);
      fmpz_zero(anum);
      fmpz_set(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));

      nf_elem_reduce(a, nf);
      nf_elem_canonicalise(a, nf);
  }
  else
  {
      fmpq_poly_shift_left(NF_ELEM(a), NF_ELEM(b), 1);
      nf_elem_reduce(a, nf);
      nf_elem_canonicalise(a, nf);
  }
}
