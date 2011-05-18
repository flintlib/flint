/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 William Hart
   
******************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

int _fmpz_poly_divides(fmpz * q, const fmpz * a, 
                       long len1, const fmpz * b, long len2)
{
   fmpz * r = _fmpz_vec_init(len1);
   
   _fmpz_poly_divrem(q, r, a, len1, b, len2);
   while ((len1) && r[len1 - 1] == 0) len1--;
   
   _fmpz_vec_clear(r, len1);

   return (len1 == 0);
}

int fmpz_poly_divides(fmpz_poly_t q, const fmpz_poly_t a, const fmpz_poly_t b)
{
   long len_out;
   int res;

   if (b->length == 0)
   {
       printf("Exception: division by zero in fmpz_poly_divides\n");
       abort();
   }

   if (a->length == 0)
   {
      fmpz_poly_zero(q);

      return 1;
   }

   if (a->length < b->length)
      return 0;

   len_out = a->length - b->length + 1;
      
   if (q == a || q == b)
   {
      fmpz_poly_t temp;
      
      fmpz_poly_init2(temp, len_out);
      res = _fmpz_poly_divides(temp->coeffs, a->coeffs, a->length, b->coeffs, b->length);
      temp->length = len_out;
      _fmpz_poly_normalise(temp);
      fmpz_poly_swap(q, temp);
      fmpz_poly_clear(temp);
   } else
   {
      fmpz_poly_fit_length(q, len_out);
      res = _fmpz_poly_divides(q->coeffs, a->coeffs, a->length, b->coeffs, b->length);
      q->length = len_out;
      _fmpz_poly_normalise(q);
   }

   return res;
}