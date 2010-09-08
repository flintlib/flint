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

   Copyright (C) 2010 William Hart
   Copyright (C) 2010 Daniel Woodhouse
   
*****************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mpoly.h"

ulong nmod_mpoly_get_coeff_of_product(nmod_mpoly_t poly1, nmod_mpoly_t poly2, ulong mon){

   long i = 0;
   long j = poly2->length -1;
   ulong bits1, bits2;
   ulong a1 = 0;
   ulong a2 = 0;
   ulong a3 = 0;
   ulong n1 = 0;
   ulong n2 = 0;
   ulong cy = 0;
   //Has it overflowed?
   int b = 0;
    
   bits1 = nmod_mpoly_max_bits(poly1);
   bits2 = nmod_mpoly_max_bits(poly2);

   b = bits1 + bits2 + FLINT_BIT_COUNT(poly2->length); //?? is the FLINT_BIT_COUNT(len2) neccessary? -YES
   if (b <= FLINT_BITS) b = 0; // output will fit in a mp_limb_t
   else b = 2; // unsigned <= 3 limbs

   while(i < poly1->length && j >= 0){
       
      if(poly1->exps[i] + poly2->exps[j] == mon){
         
         //if it fits in a single limb
         if(b == 0){
         	a1 += poly1->coeffs[i] * poly2->coeffs[j];
         }
         //if it does not fit in a single limb
         else if(b == 2){
            umul_ppmm(n2, n1, poly1->coeffs[i], poly2->coeffs[j]);
            add_ssaaaa(cy, a1, (mp_limb_t) 0, a1, (mp_limb_t) 0, n1);
            add_ssaaaa(cy, a2, (mp_limb_t) 0, a2, (mp_limb_t) 0, n2 + cy); // cannot overflow  
            if (cy) // has overflowed
	       {
		b = 3;
		a3 = cy;
	       }
         } else // general case
            {
            umul_ppmm(n2, n1, poly1->coeffs[i], poly2->coeffs[j]);
                        
            add_ssaaaa(cy, a1, (mp_limb_t) 0, a1, (mp_limb_t) 0, n1);
            add_ssaaaa(a3, a2, a3, a2, (mp_limb_t) 0, n2 + cy); // cannot overflow
         }
         

         i++;
         j = j-1; 
      } else if(poly1->exps[i] + poly2->exps[j] > mon){
         j = j-1;
         }
      else{
         i++;
         }
   }
      //now return modded value
   if (b == 0) // fits into a small fmpz
      {
         
      NMOD_RED(a1, a1, poly1->mod);
      return a1; 
   } else if (b == 2) // fits into two limbs unsigned
      {
		
      NMOD2_RED2(a1, a2, a1, poly1->mod);
      return a1;
   } else
   {
      NMOD_RED3(a1, a3, a2, a1, poly1->mod); 
      return a1;  
   }

   return 0;
}
