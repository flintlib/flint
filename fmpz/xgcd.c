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

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
{
   fmpz_t t1, t2;
   fmpz * f1, * g1;

   if (fmpz_is_zero(f))
   {
      fmpz_set(d, g);
      fmpz_set_ui(a, 0);
      fmpz_set_ui(b, 1);
   } else if (fmpz_cmp(f, g) == 0)
   {
      fmpz_set(d, f);
      fmpz_set_ui(a, 1);
      fmpz_set_si(b, 0);
   } else
   {
      int sign1 = fmpz_sgn(f);
      int sign2 = fmpz_sgn(g);

      fmpz_init(t1);
      fmpz_init(t2);

      /* support aliasing */
      if (d == f || a == f || sign1 < 0)
      {
         f1 = t1;
         if (sign1 < 0) fmpz_neg(f1, f);
         else fmpz_set(f1, f);
      } else
         f1 = (fmpz *) f;
      
      if (d == g || a == g || sign2 < 0)
      {
         g1 = t2;
         if (sign2 < 0) fmpz_neg(g1, g);
         else fmpz_set(g1, g);
      } else
         g1 = (fmpz *) g;
      
      if (fmpz_cmp(f1, g1) < 0)
      {
         fmpz_gcdinv(d, a, f1, g1);
         fmpz_mul(t1, a, f1);
         fmpz_sub(t1, d, t1);
         fmpz_divexact(b, t1, g1);
      } else /* g < f */
      {
         fmpz_gcdinv(d, b, g1, f1);
         fmpz_mul(t2, b, g1);
         fmpz_sub(t2, d, t2);
         fmpz_divexact(a, t2, f1);
      }

      if (sign1 < 0) fmpz_neg(a, a);
      if (sign2 < 0) fmpz_neg(b, b);

      fmpz_clear(t1);
      fmpz_clear(t2);
   } 
}
