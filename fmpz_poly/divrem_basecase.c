/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2008, 2009 William Hart
   
*****************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A, long A_len, 
								const fmpz * B, long B_len)
{
   if (B_len == 0)
   {
      printf("Exception : Divide by zero in _fmpz_poly_divrem_basecase.\n");
      abort();
   }
   
   long coeff = A_len;
   
   fmpz * qB;
   
   const fmpz * B_lead = B + B_len - 1; 
   fmpz * coeff_Q = Q + A_len - B_len;
   
   int want_rem = (R == NULL) ? 0 : 1;
   
   while (coeff >= B_len)
   {
      if (fmpz_cmpabs(A + coeff - 1, B_lead) >= 0) break;
      else 
	  {
		 fmpz_zero(coeff_Q);
		 coeff_Q--;
		 coeff--; 
	  }
   }
   
   if (want_rem) _fmpz_vec_copy(R, A, A_len);
   
   if (coeff < B_len)
      return;
    
   if (!want_rem) 
   {
      R = _fmpz_vec_init(coeff);
	  mpn_copyi(R + B_len - 1, A + B_len - 1, coeff - B_len + 1);
   }
   
   long B1;
   if (want_rem) B1 = B_len;
   else B1 = B_len - 1;
   long B2 = B_len;
   
   long B1_orig = B1;
   qB = _fmpz_vec_init(B1);
   
   while (coeff >= B_len)
   {
      if (fmpz_cmpabs(R + coeff - 1, B_lead) < 0) fmpz_zero(coeff_Q);
      else
      {
         fmpz_fdiv_q(coeff_Q, R + coeff - 1, B_lead);
         
         _fmpz_vec_scalar_mul_fmpz(qB, B, B1, coeff_Q); 
         
         fmpz * R_sub = R + coeff - B2;
         _fmpz_vec_sub(R_sub, R_sub, qB, B1);
      }
      
      if ((!want_rem) && (B1 >= coeff - B_len + 1))
      {
         B++;
         B1--;
         B2--;
      }
      
      coeff--;
	  coeff_Q--;
   }
         
   _fmpz_vec_clear(qB, B1_orig);
}
