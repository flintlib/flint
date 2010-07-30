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
   
*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_mul_KS(mp_ptr output, mp_srcptr input1, long length1, 
					   mp_srcptr input2, long length2, ulong bits_in, nmod_t mod)
{   
   long len_out = length1 + length2 - 1;    
   long bits; 

   if (!bits_in) 
   {
      long bits1 = _nmod_vec_max_bits(input1, length1);
      long bits2 = (input1 == input2) ? bits1 : _nmod_vec_max_bits(input2, length2);
   
      long log_length = FLINT_BIT_COUNT(length2);
      
	  bits = bits1 + bits2 + log_length;
   } else
	  bits = (long) bits_in;
   
   mp_ptr mpn1, mpn2, res;

   long limbs1 = (length1 * bits - 1) / FLINT_BITS + 1;
   long limbs2 = (length2 * bits - 1) / FLINT_BITS + 1;
      
   mpn1 = (mp_ptr) malloc(sizeof(mp_limb_t)*limbs1);
   mpn2 = (input1 == input2) ? mpn1 : (mp_ptr) malloc(sizeof(mp_limb_t)*limbs2);

   _nmod_poly_bit_pack(mpn1, input1, length1, bits);
   
   if (input1 != input2)
      _nmod_poly_bit_pack(mpn2, input2, length2, bits);
   
   res = (mp_ptr) malloc(sizeof(mp_limb_t)*(limbs1 + limbs2));
   
   if (input1 != input2) mpn_mul(res, mpn1, limbs1, mpn2, limbs2);
   else mpn_mul_n(res, mpn1, mpn1, limbs1);
   
   free(mpn2);
   if(input1 != input2)
      free(mpn1);
  
   _nmod_poly_bit_unpack(output, res, len_out, bits, mod); 
   
   free(res);
}

void nmod_poly_mul_KS(nmod_poly_t res, 
                         const nmod_poly_t poly1, const nmod_poly_t poly2, ulong bits_in)
{
   long len_out;
   
   if ((poly1->length == 0) || (poly2->length == 0))  
   {
      nmod_poly_zero(res);
      return;
   }

   len_out = poly1->length + poly2->length - 1;

   if (res == poly1 || res == poly2)
   {
	   nmod_poly_t temp;
	   nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
	   if (poly1->length >= poly2->length)
		  _nmod_poly_mul_KS(temp->coeffs, poly1->coeffs, poly1->length, 
		              poly2->coeffs, poly2->length, bits_in, poly1->mod);
	   else
		  _nmod_poly_mul_KS(temp->coeffs, poly2->coeffs, poly2->length, 
		              poly1->coeffs, poly1->length, bits_in, poly1->mod);
	   nmod_poly_swap(res, temp);
	   nmod_poly_clear(temp);
   } else
   {
	   nmod_poly_fit_length(res, len_out);
	   if (poly1->length >= poly2->length)
		  _nmod_poly_mul_KS(res->coeffs, poly1->coeffs, poly1->length, 
		               poly2->coeffs, poly2->length, bits_in, poly1->mod);
	   else
		  _nmod_poly_mul_KS(res->coeffs, poly2->coeffs, poly2->length, 
		               poly1->coeffs, poly1->length, bits_in, poly1->mod);
   }

   res->length = len_out;
   _nmod_poly_normalise(res);
}
