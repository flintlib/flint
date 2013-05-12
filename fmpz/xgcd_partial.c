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

void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, 
                                    fmpz_t r2, fmpz_t r1, const fmpz_t L)
{
   fmpz_t q, r;
   long aa2, aa1, bb2, bb1, rr1, rr2, qq, bb, t1, t2, t3, i;
   long bits, bits1, bits2;

   fmpz_init(q); fmpz_init(r);
   
   fmpz_zero(co2);
   fmpz_set_si(co1, -1);
  
   while (!fmpz_is_zero(r1) && fmpz_cmp(r1, L) > 0)
   {
      bits2 = fmpz_bits(r2);
      bits1 = fmpz_bits(r1);
      bits = FLINT_MAX(bits2, bits1) - FLINT_BITS + 1;
      if (bits < 0) bits = 0;
      
      fmpz_tdiv_q_2exp(r, r2, bits);
      rr2 = fmpz_get_ui(r);
      fmpz_tdiv_q_2exp(r, r1, bits);
      rr1 = fmpz_get_ui(r);
      fmpz_tdiv_q_2exp(r, L, bits);
      bb = fmpz_get_ui(r);

      aa2 = 0; aa1 = 1;
      bb2 = 1; bb1 = 0;

      for (i = 0; rr1 != 0 && rr1 > bb; i++)
      {
         qq = rr2 / rr1;

         t1 = rr2 - qq*rr1; 
         t2 = aa2 - qq*aa1; 
         t3 = bb2 - qq*bb1; 

         if (i & 1)
         {
            if (t1 < -t3 || rr1 - t1 < t2 - aa1) break;
         } else 
         {
            if (t1 < -t2 || rr1 - t1 < t3 - bb1) break;
         }

         rr2 = rr1; rr1 = t1;
         aa2 = aa1; aa1 = t2;
         bb2 = bb1; bb1 = t3;
      }

      if (i == 0)
      {
         fmpz_fdiv_qr(q, r2, r2, r1);
         fmpz_swap(r2, r1);

         fmpz_submul(co2, co1, q);
         fmpz_swap(co2, co1);
      } else
      {
         fmpz_mul_si(r, r2, bb2);
         if (aa2 >= 0)
            fmpz_addmul_ui(r, r1, aa2);
         else
            fmpz_submul_ui(r, r1, -aa2);
         fmpz_mul_si(r1, r1, aa1);
         if (bb1 >= 0)
            fmpz_addmul_ui(r1, r2, bb1);
         else
            fmpz_submul_ui(r1, r2, -bb1);
         fmpz_set(r2, r);

         fmpz_mul_si(r, co2, bb2);
         if (aa2 >= 0)
            fmpz_addmul_ui(r, co1, aa2);
         else
            fmpz_submul_ui(r, co1, -aa2);
         fmpz_mul_si(co1, co1, aa1);
         if (bb1 >= 0)
            fmpz_addmul_ui(co1, co2, bb1);
         else
            fmpz_submul_ui(co1, co2, -bb1);
         fmpz_set(co2, r);

         if (fmpz_sgn(r1) < 0) { fmpz_neg(co1, co1); fmpz_neg(r1, r1); }
         if (fmpz_sgn(r2) < 0) { fmpz_neg(co2, co2); fmpz_neg(r2, r2); }
      }
   }

   if (fmpz_sgn(r2) < 0)
   { 
      fmpz_neg(co2, co2); fmpz_neg(co1, co1);
      fmpz_neg(r2, r2);
   }

   fmpz_clear(q); fmpz_clear(r);
}
