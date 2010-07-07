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

void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, ulong * d, const fmpz * A, ulong A_len, 
								const fmpz * B, ulong B_len)
{
   if (B_len == 0)
   {
      printf("Exception : Divide by zero in _fmpz_poly_pseudo_divrem_basecase.\n");
      abort();
   }
   
   ulong coeff = A_len;
   
   fmpz * qB;
   *d = 0;
   
   const fmpz * B_lead = B + B_len - 1; 
   fmpz * coeff_Q = Q + A_len - B_len;
   fmpz_t rem;
   int scale;
   ulong i;

   _fmpz_vec_copy(R, A, A_len);
   
   if (coeff < B_len)
      return;
   
   fmpz_init(rem);

   ulong B1 = B_len - 1;
   ulong B1_orig = B1;
   qB = _fmpz_vec_init(B1);

   for (i = 0; i < A_len - B_len + 1; i++)
      fmpz_zero(Q + i);
   
   while (coeff >= B_len)
   {
      fmpz_fdiv_qr(coeff_Q, rem, R + coeff - 1, B_lead);
      if (fmpz_is_zero(rem)) scale = 0;
      else
	  {
		 _fmpz_vec_scalar_mul_fmpz(Q, Q, A_len - B_len + 1, B_lead);
		 fmpz_set(coeff_Q, R + coeff - 1);
         scale = 1;
         (*d)++;
      }
      
      if (B_len > 1) _fmpz_vec_scalar_mul_fmpz(qB, B, B1, coeff_Q); 
         
	  if (scale) _fmpz_vec_scalar_mul_fmpz(R, R, A_len, B_lead);
      
      if (B_len > 1)
         _fmpz_vec_sub(R + coeff - B_len, R + coeff - B_len, qB, B_len - 1);
      
	  fmpz_zero(R + coeff - 1);
      
      coeff--;
	  coeff_Q--;
   }
         
   fmpz_clear(rem);
   _fmpz_vec_clear(qB, B1_orig);
}

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
								 ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (B->length == 0)
   {
	   printf("Exception: division by zero in fmpz_poly_pseudo_divrem_basecase.\n");
       abort();
   }

   if (A->length < B->length)
   {
      fmpz_poly_zero(Q);
	  fmpz_poly_set(R, A);
	  *d = 0;
	  return;
   }
   
   fmpz_poly_t t1, t2;
   fmpz * Q_coeffs, * R_coeffs;
   
   if ((Q == A) || (Q == B))
   {
	   fmpz_poly_init2(t1, A->length - B->length + 1);
	   Q_coeffs = t1->coeffs;
   } else 
   {
	   fmpz_poly_fit_length(Q, A->length - B->length + 1);
       Q_coeffs = Q->coeffs;
   }
	  
   if ((R == A) || (R == B))
   {
	   fmpz_poly_init2(t2, A->length);
	   R_coeffs = t2->coeffs;
   } else    
   {
	   fmpz_poly_fit_length(R, A->length);
       R_coeffs = R->coeffs;
   }

   _fmpz_poly_pseudo_divrem_basecase(Q_coeffs, R_coeffs, d, A->coeffs, A->length, B->coeffs, B->length);
	
   if ((Q == A) || (Q == B))
   {
	   _fmpz_poly_set_length(t1, A->length - B->length + 1);
       fmpz_poly_swap(t1, Q);
	   fmpz_poly_clear(t1);
   } else
       _fmpz_poly_set_length(Q, A->length - B->length + 1);

   if ((R == A) || (R == B))
   {
	   _fmpz_poly_set_length(t2, A->length);
       fmpz_poly_swap(t2, R);
	   fmpz_poly_clear(t2);
   } else
	   _fmpz_poly_set_length(R, A->length);
   
   _fmpz_poly_normalise(Q);	  
   _fmpz_poly_normalise(R);	  
}
