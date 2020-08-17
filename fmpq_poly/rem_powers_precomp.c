/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void _fmpq_poly_rem_powers_precomp(fmpz * A, fmpz_t denA, slong m, 
                                 const fmpz * B, const fmpz_t denB, slong n, 
                                 fmpq_poly_struct * const powers)
{
   slong i;
   fmpq_poly_t prod;
   fmpz_t den;

   if (m >= 2*n) /* TODO: make this case use the precomputed squares */
   {
      fmpz * R = _fmpz_vec_init(m);
      fmpz_init(den);

      _fmpz_vec_set(R, A, m);
      fmpz_set(den, denA);

      _fmpq_poly_rem(A, denA, R, den, m, B, denB, n, NULL);

      _fmpz_vec_clear(R, m);
      fmpz_clear(den);

      return;
   }

   if (m < n)
      return;

   fmpz_init(den);
   
   fmpq_poly_init2(prod, n - 1);
   fmpz_set(den, denA);

   for (i = n - 1; i < m; i++)
   {
      _fmpz_vec_scalar_mul_fmpz(fmpq_poly_numref(prod), 
         fmpq_poly_numref(powers + i), powers[i].length, A + i);
      fmpz_mul(fmpq_poly_denref(prod), fmpq_poly_denref(powers + i), den);
      _fmpq_poly_add_can(A, denA, A, denA, n - 1, fmpq_poly_numref(prod), 
         fmpq_poly_denref(prod), powers[i].length, 0);
   }

   fmpq_poly_clear(prod);
   fmpz_clear(den);
}

void 
fmpq_poly_rem_powers_precomp(fmpq_poly_t R, const fmpq_poly_t A, 
                    const fmpq_poly_t B, const fmpq_poly_powers_precomp_t B_inv)
{
    fmpq_poly_t tR;
    fmpz * r, * d;
    slong len1 = A->length, len2 = B->length;
    
    if (len1 < len2)
    {
        fmpq_poly_set(R, A);
        return;
    }

    if (R == B)
    {
        fmpq_poly_init2(tR, len1);
        r = fmpq_poly_numref(tR);
        d = fmpq_poly_denref(tR);
    }
    else
    {
        fmpq_poly_fit_length(R, len1);
        r = fmpq_poly_numref(R);
        d = fmpq_poly_denref(R);
    }

    if (R == B || R != A)
    {
       _fmpz_vec_set(r, fmpq_poly_numref(A), len1);
       fmpz_set(d, fmpq_poly_denref(A));
    }

    _fmpq_poly_rem_powers_precomp(r, d, len1, fmpq_poly_numref(B),
                        fmpq_poly_denref(B), len2, B_inv->powers);

    if (R == B)
    {
        _fmpq_poly_set_length(tR, len2 - 1);
        fmpq_poly_swap(tR, R);
        fmpq_poly_clear(tR);
    }
    else
        _fmpq_poly_set_length(R, len2 - 1);

    _fmpq_poly_normalise(R);
}
