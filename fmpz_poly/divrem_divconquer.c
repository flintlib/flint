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

void _fmpz_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
				const fmpz * A, long A_len, const fmpz * B, long B_len)
{
   if (A_len < 2*B_len - 1)
   {
      /* Convert unbalanced division into a 2*q - 1 by q division */
      const fmpz * t_A, * t_B, * t_B2;
	  fmpz * d1q1, * d2q1;
      
      long q = A_len - B_len + 1;
      long q2 = B_len - q;
      
	  t_A = A + q2;
      t_B = B + q2;
      t_B2 = B;
      
      d1q1 = R + q2;
	  _fmpz_poly_divrem_divconquer_recursive(Q, d1q1, t_A, t_B, q); 
      
      /*
         Compute d2q1 = Q*t_B2
         It is of length q2+q-1
      */
      
      d2q1 = _fmpz_vec_init(q2 + q - 1);
	  if (q >= q2)
		 _fmpz_poly_mul(d2q1, Q, q, t_B2, q2);
	  else
		 _fmpz_poly_mul(d2q1, t_B2, q2, Q, q);

      /*
         Compute BQ = d1q1*x^n1 + d2q1
         It has length n1+n2-1
		 Then compute R = A-BQ
      */
      
      _fmpz_vec_copy(R, d2q1, q2);
	  _fmpz_vec_add(R + q2, R + q2, d2q1 + q2, q - 1);
	  _fmpz_vec_sub(R, A, R, A_len);
	  
	  _fmpz_vec_clear(d2q1, q2 + q - 1);
            
      return;   
   } 

   if (A_len > 2*B_len - 1)
   {
      const fmpz * p1;
	  fmpz * d1q1, * dq1, * q1, * q2;
	  
	  // We shift A right until it is length 2*B->length -1
      // We call this polynomial p1
      
      long shift = A_len - 2*B_len + 1;
      p1 = A + shift;
	  
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being length B->length
         d1q1 = d1*q1 is length 2*B->length-1
      */
      
      dq1 = _fmpz_vec_init(A_len);
	  d1q1 = dq1 + shift;
	  q1 = Q + shift;
	  
      _fmpz_poly_divrem_divconquer_recursive(q1, d1q1, p1, B, B_len); 
       
      /* 
         We have dq1 = d1*q1*x^shift
         dq1 is of lengthA->length
 
         Compute R = A - dq1 
         The first B->length coefficients represent
		 remainder terms (zero if division is exact), 
		 leaving A->length - B->length significant 
		 terms which we use in the division
      */
   
      _fmpz_vec_copy(dq1, A, shift);
	  _fmpz_vec_sub(dq1 + shift, A + shift, dq1 + shift, B_len - 1);
      _fmpz_vec_sub(R + A_len - B_len, A + A_len - B_len, dq1 + A_len - B_len, B_len);

      /*
         Compute q2 = trunc(R) div B
         It is a smaller division than the original 
         since trunc(R)->length = A->length - B->length
      */
   
      q2 = Q;
	  _fmpz_poly_divrem_divconquer(q2, R, dq1, A_len - B_len, B, B_len); 
      
      /*
         We have Q = q1*x^shift + q2
         Q has length B->length+shift
         Note q2 has length shift since 
         it is an A->length-B->length 
         by B->length division
      
         We've also written the 
		 remainder in place
      */
      
      _fmpz_vec_clear(dq1, A_len);
      
      return;
   } 

   /* A_len = 2*B_len - 1 */
   fmpz * QB = _fmpz_vec_init(A_len);
   
   _fmpz_poly_divrem_divconquer_recursive(Q, QB, A, B, B_len);
   _fmpz_vec_sub(R, A, QB, A_len);
   
   _fmpz_vec_clear(QB, A_len);
}

void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, 
								 const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (B->length == 0)
   {
	   printf("Exception: division by zero in fmpz_poly_divrem_divconquer.\n");
       abort();
   }

   if (A->length < B->length)
   {
      fmpz_poly_zero(Q);
	  fmpz_poly_set(R, A);
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

   _fmpz_poly_divrem_divconquer(Q_coeffs, R_coeffs, A->coeffs, A->length, B->coeffs, B->length);
	
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
