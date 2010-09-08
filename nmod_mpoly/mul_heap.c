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

   This implementation was inspired by the reference implementation of
   Roman Pearce.
   
*****************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mpoly.h"
 
#include "longlong.h"
#include "nmod_vec.h"


void nmod_mpoly_mul_heap(nmod_mpoly_t res, nmod_mpoly_t poly1, nmod_mpoly_t poly2)
{
   
   ulong len1 = poly1->length;
   ulong len2 = poly2->length;
   ulong i, j, bits1 = 0, bits2 = 0, b, b_old;

   nmod_mpoly_heap_t * heap = (nmod_mpoly_heap_t *) malloc((len2 + 1)*sizeof(nmod_mpoly_heap_t));
   nmod_mpoly_entry_t * entries = (nmod_mpoly_entry_t *) malloc(len2*sizeof(nmod_mpoly_entry_t));
   
   ulong n; // index of final entry in heap (0 means nothing in heap yet as indices start at 1)
   ulong len_out = 0; // length of product polynomial
   ulong res_alloc = res->alloc;
   
   //fmpz_init(temp); Now redundant??

   /* compute the maximum bitsize for poly1 and poly2 */
   bits1 = nmod_mpoly_max_bits(poly1);
   bits2 = nmod_mpoly_max_bits(poly2);

   b = bits1 + bits2 + FLINT_BIT_COUNT(len2); //?? is the FLINT_BIT_COUNT(len2) neccessary? -YES
   if (b <= FLINT_BITS) b = 0; // output will fit in a mp_limb_t
   else b = 2; // unsigned <= 3 limbs
  
   b_old = b;
   
   /* place all products poly1[0]*poly2[j] into our array of entries */
  
   for (j = 0; j < len2; j++)
   {
      
      entries[j].next = NULL;
      entries[j].i1 = 0;
      entries[j].i2 = j;
      
   }  

   /* 
      insert only the first of these entries into the heap - we inject the others later 
	  we do this manually as the heap_insert can't handle an empty heap
   */

   heap[1].entry = entries; //N.B. :provides a pointer to the first entry 
   heap[1].exp = poly1->exps[0] + poly2->exps[0];
   n = 1;

   /*
      this is where we perform the heap multiplication
   */
   while (n)
   {
          
          nmod_mpoly_entry_t * chain = heap[1].entry; // chain starting from top of heap
	  nmod_mpoly_entry_t * next = chain; // to iterate along chain
          nmod_mpoly_entry_t * prev;
          mp_limb_t exp = heap[1].exp; // current exponents
	  mp_limb_t c1, c2; // coefficients
          mp_limb_t accum0 = 0, accum1 = 0, accum2 = 0, cy;
	  
	  //int sign;
	  ulong size;
          mp_limb_t n1, n2, n3; 

	  b = b_old;
      
	  /* make sure there's enough space for the new term */
	  if (len_out >= res_alloc)
	  {
	         nmod_mpoly_fit_length(res, len_out + 1);
		 res_alloc = res->alloc;
		 
	  }
          
	  /* chain all heap elements with the same exp together and remove from heap */
	  while (1)
	  {
	     prev = next;
	     next = next->next;
	     if (next) continue;
		 
		 // end of a chain
		 nmod_mpoly_reheapify(heap, &n); // remove from heap
	     
		 if ((n) && (exp == heap[1].exp)) // next element the same exp?
		 {
			prev->next = heap[1].entry;
			next = prev->next;
	     } else break; // chain is complete
	  }
          
          next = chain; // go back to start of chain
	  
	  /* add together all the terms in the chain.*/
	  do
	  {      
                 
                 prev = next;
		 
		 //get the coefficients of the corresponding terms
		 
		 c1 = poly1->coeffs[next->i1];
		 
		 c2 = poly2->coeffs[next->i2];
                 
		 /* addmul */
		 if (b == 0) // fits in a single limb
		 {
		    accum0 += (c1 * c2);
		 } else if (b == 2) // input coeffs fit in a small and are unsigned
		 {
		    umul_ppmm(n2, n1, c1, c2);
                    add_ssaaaa(cy, accum0, (mp_limb_t) 0, accum0, (mp_limb_t) 0, n1);
                    add_ssaaaa(cy, accum1, (mp_limb_t) 0, accum1, (mp_limb_t) 0, n2 + cy); // cannot overflow  
		    if (cy) // has overflowed
		    {
			b = 3;
			accum2 = cy;
		    }
		 } else // general case
		 {
			
		    umul_ppmm(n2, n1, c1, c2);
                        
                    add_ssaaaa(cy, accum0, (mp_limb_t) 0, accum0, (mp_limb_t) 0, n1);
                    add_ssaaaa(accum2, accum1, accum2, accum1, (mp_limb_t) 0, n2 + cy); // cannot overflow
 		 }
             
	     if (next->i1 == 0 && next->i2 < len2 - 1) // inject next entry poly1[0]*poly2[j+1]
	        nmod_mpoly_heap_insert(heap, &n, next + 1, poly1->exps[0] + poly2->exps[next->i2 + 1]);
		
		          
             next = next->next;
		  
             /* insert poly1[i+1]*poly2[j] */
	     prev->i1++;
		 
	     if (prev->i1 < len1)
		 nmod_mpoly_heap_insert(heap, &n, prev, poly1->exps[prev->i1] + poly2->exps[prev->i2]);
	  } while (next); 
	  //reached the end of the chain


	  /* store term */
	  if (b == 0) // fits into a small fmpz
	  {
	         
             NMOD_RED(res->coeffs[len_out], accum0, res->mod);
		 
	  } else if (b == 2) // fits into two limbs unsigned
	  {		
             NMOD2_RED2(res->coeffs[len_out], accum1, accum0, res->mod);
	  } else
	  {
             NMOD_RED(accum2, accum2, res->mod); //is this neccessary?
	     NMOD_RED3(res->coeffs[len_out], accum2, accum1, accum0, res->mod);    
	  }

	  res->exps[len_out] = exp;
	  if (res->coeffs[len_out] != 0) //So that zero terms are ignored
	     len_out++;
   }
   //End of multiplication

   free(heap);
   free(entries);

   res->length = len_out;
}
