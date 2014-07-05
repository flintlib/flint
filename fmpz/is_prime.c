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

    Copyright (C) 2014 William Hart
   
******************************************************************************/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int fmpz_is_prime(const fmpz_t n)
{
   int res = 0;
   
   if (fmpz_cmp_ui(n, 1) <= 0)
      return 0;
   
   if (fmpz_is_even(n))
      return (fmpz_cmp_ui(n, 2) == 0);

   if (fmpz_is_square(n))
      return 0;

   {
      double logd = log(fmpz_get_d(n));
      ulong limit = (ulong) (logd*logd*logd/10.0) + 2;

      fmpz_t F, F2, F3;

      fmpz_init(F);
      fmpz_init(F2);
      fmpz_init(F3);

      res = fmpz_is_prime_pocklington(F, n, limit);

      if (res == 1)
      {
         fmpz_mul(F2, F, F);
         if (fmpz_cmp(F2, n) < 0)
         {
            fmpz_mul(F3, F2, F);
            if (fmpz_cmp(F3, n) >= 0) /* Brillhart, Lehmer, Selfridge test */
            {
               fmpz_t n1, c2, c1;

               fmpz_init(n1);
               fmpz_init(c2);
               fmpz_init(c1);

               fmpz_sub_ui(n1, n, 1); /* n is 1 mod F */
               fmpz_tdiv_q(n1, n1, F);

               fmpz_tdiv_qr(c2, c1, n1, F); /* Let n = c2*F^2 + c1*F + 1 */

               fmpz_mul(c1, c1, c1); /* check if c1^2 - 4*c2 is a square */
               fmpz_submul_ui(c1, c2, 4);

               if (fmpz_is_square(c1))
                  res = 0;
               /* else n is prime (res == 1) */

               fmpz_clear(n1);
               fmpz_clear(c2);
               fmpz_clear(c1);
            } else /* p + 1 test */
            {
               fmpz_t FF;
         
               fmpz_init(FF);
         
               res = fmpz_is_prime_morrison(FF, n, limit);

               if (res == 1)
               {
                  fmpz_sub_ui(FF, FF, 1); /* need F - 1 > sqrt(n) */
                  fmpz_mul(F2, FF, FF);
                  if (fmpz_cmp(F2, n) <= 0)
                     res = -1;
               }

               fmpz_clear(FF);
            }
         }
      } 

      fmpz_clear(F);
      fmpz_clear(F2);
      fmpz_clear(F3);
   }

   return res;
      
}
