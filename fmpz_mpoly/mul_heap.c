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
   ulong i, j, bits1 = 0, bits2 = 0, b, b_old;
   int signed_c = 0;
   fmpz_t temp;

   fmpz_mpoly_heap_t * heap = (fmpz_mpoly_heap_t *) malloc((len2 + 1)*sizeof(fmpz_mpoly_heap_t));
   fmpz_mpoly_entry_t * entries = (fmpz_mpoly_entry_t *) malloc(len2*sizeof(fmpz_mpoly_entry_t));
   
   ulong n; // index of final entry in heap (0 means nothing in heap yet as indices start at 1)
   ulong len_out = 0; // length of product polynomial
   ulong res_alloc = res->alloc;
   
   fmpz_init(temp);

   /* compute the maximum bitsize for poly1 and poly2 */
   for (i = 0; i < len1; i++)
   {
      b = fmpz_bits(poly1->coeffs + i);
	  if (b > bits1) bits1 = b;
	  if (fmpz_sgn(poly1->coeffs + i) < 0) signed_c = 1;
   }

   for (i = 0; i < len2; i++)
   {
      b = fmpz_bits(poly2->coeffs + i);
	  if (b > bits2) bits2 = b;
	  if (fmpz_sgn(poly2->coeffs + i) < 0) signed_c = 1;
   }

   b = bits1 + bits2 + FLINT_BIT_COUNT(len2);
   if (b < FLINT_BITS - 2) b = 0; // output will fit in a small fmpz
   else if (!signed_c && bits1 <= FLINT_BITS - 2 && bits2 <= FLINT_BITS - 2) b = 2; // unsigned <= 3 limbs
   else b = 3;

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
	  fmpz c1, c2; // coefficients
      mp_limb_t accum0 = 0, accum1 = 0, accum2 = 0, * ptr, cy;
	  __mpz_struct * coeff;
	  int sign;
	  ulong size;
      mp_limb_signed_t n1, n2, n3;

	  b = b_old;
      
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
		 c1 = poly1->coeffs[next->i1];
		 c2 = poly2->coeffs[next->i2];
         
		 /* addmul */
		 if (b == 0) // fits in a small fmpz
		 {
		    accum0 += ((mp_limb_signed_t) c1 * (mp_limb_signed_t) c2);
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
			if (!COEFF_IS_MPZ(c1) && !COEFF_IS_MPZ(c2))
		    {
		       smul_ppmm(n2, n1, c1, c2);
               n3 = (n2>>(FLINT_BITS - 1)); //sign extend
               add_ssaaaa(cy, accum0, (mp_limb_t) 0, accum0, (mp_limb_t) 0, n1);
               add_ssaaaa(accum2, accum1, accum2, accum1, (mp_limb_t) 0, cy);
               add_ssaaaa(accum2, accum1, accum2, accum1, n3, n2);
		    } else
		       fmpz_addmul(temp, poly1->coeffs + next->i1, poly2->coeffs + next->i2);
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
	  if (b == 0) // fits into a small fmpz
	  {
	     _fmpz_demote(res->coeffs + len_out);
		 res->coeffs[len_out] = accum0;
	  } else if (b == 2) // fits into two limbs unsigned
	  {		
		 if (accum1 != 0 || accum0 > COEFF_MAX)
		 {
			if (!COEFF_IS_MPZ(res->coeffs[len_out])) 
			   coeff = _fmpz_promote(res->coeffs + len_out);
			if (coeff->_mp_alloc < 2) mpz_realloc(coeff, 2);
			coeff->_mp_d[0] = accum0;
			coeff->_mp_d[1] = accum1;
			coeff->_mp_size = 2 - (accum1 == 0);
		 } else 
		 {
			_fmpz_demote(res->coeffs + len_out);
			if (accum0 == 0) continue;
			res->coeffs[len_out] = accum0;
		 }
	  } else
	  {
		 if (sign = (accum2 < (mp_limb_signed_t) 0)) // if negative, negate
	     {
		    add_ssaaaa(cy, accum0, (mp_limb_t) 0, ~accum0, (mp_limb_t) 0, (mp_limb_t) 1);
            add_ssaaaa(accum2, accum1, ~accum2, ~accum1, (mp_limb_t) 0, cy);
	     } 

	     if (accum2 != (mp_limb_t) 0) size = 3; // normalise
		 else if (accum1 != (mp_limb_t) 0) size = 2;
	     else 
		 {
			 _fmpz_demote(res->coeffs + len_out);
			 if (accum0 != (mp_limb_t) 0) size = 1;
		     else 
			 {
				fmpz_set(res->coeffs + len_out, temp); // add rest of coeff
        		fmpz_zero(temp); // zero temp for next coefficient
  		        if (fmpz_is_zero(res->coeffs + len_out)) continue;
				else
				{
	               fmpz_set_ui(res->exps + len_out, exp);
	               len_out++;
				   continue;
				}
			 }
	     }
	      
		 if (size == 1 && accum0 <= COEFF_MAX) // fits in a small fmpz
		 {
			if (sign) accum0 = -accum0;
			res->coeffs[len_out] = accum0;
		 } else // general case
         {
		 	coeff = _fmpz_promote(res->coeffs + len_out);
			ptr = mpz_realloc(coeff, 3);
			if (sign) coeff->_mp_size = -size;
			else coeff->_mp_size = size;
			ptr[0] = accum0;
			ptr[1] = accum1;
			ptr[2] = accum2;
	     }

		 fmpz_add(res->coeffs + len_out, res->coeffs + len_out, temp); // add rest of coeff
		 fmpz_zero(temp); // zero temp for next coefficient
	  }

	  fmpz_set_ui(res->exps + len_out, exp);
	  
	  len_out++;
   }

   if (b == 3) fmpz_clear(temp);
   free(heap);
   free(entries);

   res->length = len_out;
}
