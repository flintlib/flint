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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

// Assumes poly1 and poly2 are not length 0
void _fmpz_poly_mul_classical(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   ulong len1 = poly1->length;
   ulong len2 = poly2->length;
   ulong len_out = len1 + len2 - 1;
   
   if ((len1 == 1) && (len2 == 1)) // Special case if the length of both inputs is 1
   {
      fmpz_mul(res->coeffs, poly1->coeffs, poly2->coeffs);      
   } else // Ordinary case
   {
      long i, j;
      
      // Set res[i] = poly1[i]*poly2[0] 
      if (poly2->coeffs[0])
		for (i = 0; i < len1; i++)
           fmpz_mul(res->coeffs + i, poly1->coeffs + i, poly2->coeffs);
	  else 
		for (i = 0; i < len1; i++)
           fmpz_zero(res->coeffs + i);
   
      // Set res[i+len1-1] = in1[len1-1]*in2[i]
      if (poly1->coeffs[len1 - 1])
	     for (i = 1; i < len2; i++)
            fmpz_mul(res->coeffs + i + len1 - 1, poly1->coeffs + len1 - 1, poly2->coeffs + i);  
	  else 
         for (i = 1; i < len2; i++)
            fmpz_zero(res->coeffs + i + len1 - 1);
      
      // out[i+j] += in1[i]*in2[j] 
      for (i = 0; i < len1 - 1; i++)
      {      
         fmpz c = poly1->coeffs[i];
		 if (c)
		 {
		    if (!COEFF_IS_MPZ(c))
		    {
		       if (c < 0L) 
			      for (j = 1; j < len2; j++)
                     fmpz_submul_ui(res->coeffs + i + j, poly2->coeffs + j, -c);
			   else
                  for (j = 1; j < len2; j++)
				     fmpz_addmul_ui(res->coeffs + i + j, poly2->coeffs + j, c);
		    } else
		       for (j = 1; j < len2; j++)
                  fmpz_addmul(res->coeffs + i + j, poly1->coeffs + i, poly2->coeffs + j);
	     }
      }
   } 
   
   _fmpz_poly_set_length(res, len_out);
   _fmpz_poly_normalise(res);
}

void fmpz_poly_mul_classical(fmpz_poly_t res, 
                         const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   if ((poly1->length == 0) || (poly2->length == 0))  
   {
      fmpz_poly_zero(res);
      return;
   }

   if (res == poly1 || res == poly2)
   {
	   fmpz_poly_t temp;
	   fmpz_poly_init(temp);
	   fmpz_poly_fit_length(temp, poly1->length + poly2->length - 1);
      _fmpz_poly_mul_classical(temp, poly1, poly2);
	   fmpz_poly_swap(res, temp);
	   fmpz_poly_clear(temp);
   } else
   {
	   fmpz_poly_fit_length(res, poly1->length + poly2->length - 1);
      _fmpz_poly_mul_classical(res, poly1, poly2);
   }
}
