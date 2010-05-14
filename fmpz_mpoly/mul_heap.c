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

   This implementation was inspired by the reference implementation of
   Roman Pearce.
   
*****************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_mul_heap(fmpz_mpoly_t res, fmpz_mpoly_t poly1, fmpz_mpoly_t poly2)
{
   ulong len1 = poly1->length;
   ulong len2 = poly2->length;
   ulong i, j;

   fmpz_mpoly_heap_t * heap = (fmpz_mpoly_heap_t *) malloc((len2 + 1)*sizeof(fmpz_mpoly_heap_t));
   fmpz_mpoly_entry_t * entries = (fmpz_mpoly_entry_t *) malloc(len2*sizeof(fmpz_mpoly_entry_t));
   
   ulong n; // index of final entry in heap (0 means nothing in heap yet as indices start at 1)
   ulong len_out = 0; // length of product polynomial
   ulong res_alloc = res->alloc;

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
   heap[1].entry = entries; 
   heap[1].exp = poly1->exps[0] + poly2->exps[0];
   n = 1;

   /*
      this is where we perform the heap multiplication
   */
   while (n)
   {
      fmpz_mpoly_entry_t * chain = heap[1].entry; // chain starting from top of heap
	  fmpz_mpoly_entry_t * next = chain; // to iterate along chain
      fmpz_mpoly_entry_t * prev;
      fmpz exp = heap[1].exp; // current exponents
	  mp_limb_t accum[3] = { 0, 0, 0 };
	  mp_limb_t cy;
	  int sign;
	  ulong size;

	  /* make sure there's enough space for the new term */
	  if (len_out >= res_alloc)
	  {
	     fmpz_mpoly_fit_length(res, len_out + 1);
		 res_alloc = res->alloc;
	  }

	  /* chain all heap elements with the same exp together and remove from heap */
	  while (1)
	  {
	     prev = next;
		 next = next->next;
	     if (next) continue;
		 
		 // end of a chain
		 fmpz_mpoly_reheapify(heap, &n); // remove from heap
	     
		 if ((n) && (exp == heap[1].exp)) // next element the same exp?
		 {
			prev->next = heap[1].entry;
			next = prev->next;
	     } else break; // chain is complete
	  }
      
      next = chain; // go back to start of chain
	  
	  do
	  {
         prev = next;
		 fmpz c1 = poly1->coeffs[next->i1];
		 fmpz c2 = poly2->coeffs[next->i2];
         mp_limb_signed_t n1, n2, n3;

		 /* addmul */
		 if (!COEFF_IS_MPZ(c1) && !COEFF_IS_MPZ(c2))
		 {
		    smul_ppmm(n2, n1, c1, c2);
            n3 = (n2>>(FLINT_BITS - 1)); //sign extend
            add_ssaaaa(cy, accum[0], (mp_limb_t) 0, accum[0], (mp_limb_t) 0, n1);
            add_ssaaaa(accum[2], accum[1], accum[2], accum[1], (mp_limb_t) 0, cy);
            add_ssaaaa(accum[2], accum[1], accum[2], accum[1], n3, n2);
		 } else
		 {
		    printf("not implemented yet\n");
			abort();
		 }

	     if (next->i1 == 0 && next->i2 < len2 - 1) // inject next entry poly1[1]*poly2[j+1]
		    fmpz_mpoly_heap_insert(heap, &n, next + 1, poly1->exps[0] + poly2->exps[next->i2 + 1]);
         
		 next = next->next;
		  
		 /* insert poly1[i+1]*poly2[j] */
		 prev->i1++;
         if (prev->i1 < len1)
		    fmpz_mpoly_heap_insert(heap, &n, prev, poly1->exps[prev->i1] + poly2->exps[prev->i2]);
	  } while (next);

	  /* store term */
	  if (accum[2] < (mp_limb_signed_t) 0) // negate
	  {
		 sign = 1;
		 add_ssaaaa(cy, accum[0], (mp_limb_t) 0, ~accum[0], (mp_limb_t) 0, (mp_limb_t) 1);
         add_ssaaaa(accum[2], accum[1], ~accum[2], ~accum[1], (mp_limb_t) 0, cy);
	  } else
	     sign = 0;

	  size = 3;
	  if (accum[2] == (mp_limb_t) 0) // normalise
	  {
		 size = 2;
		 if (accum[1] == (mp_limb_t) 0)
	     {
	        size = 1;
			if (accum[0] == (mp_limb_t) 0)
			   size = 0;
		 }
	  }  

	  if (size == 0) 
	     continue;
	  
	  if (size <= 1)
	  {
	     fmpz_set_ui(res->coeffs + len_out, accum[0]);
		 if (sign) fmpz_neg(res->coeffs + len_out, res->coeffs + len_out);
	  } else
	  {
		 __mpz_struct * ptr = _fmpz_promote(res->coeffs + len_out);
	     mpz_realloc(ptr, size);
		 for (i = 0; i < size; i++)
			ptr->_mp_d[i] = accum[i];
		 ptr->_mp_size = size;
         if (sign) ptr->_mp_size = -ptr->_mp_size;
	  }

	  fmpz_set_ui(res->exps + len_out, exp);
	  
	  len_out++;
   }

   free(heap);
   free(entries);

   res->length = len_out;
}
