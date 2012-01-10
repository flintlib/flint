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

    Copyright (C) 2008-2011 William Hart
    
******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fft.h"

void _fmpz_poly_mul_SS(fmpz * output, const fmpz * input1, long length1, 
                       const fmpz * input2, long length2, const long bits_in)
{
   long len_out = length1 + length2 - 1;
   long output_bits, n, limbs, size, i;
	long log_length2 = 0, log_length = 0;
   mp_limb_t * ptr, * t1, * t2, * tt, * s1, ** ii, ** jj;
   long bits1, bits2;
   int sign = 0;

   while ((1L<<log_length) < len_out) log_length++;
   
	if (!bits_in)
	{
		ulong size1 = _fmpz_vec_max_limbs(input1, length1); 
      ulong size2 = _fmpz_vec_max_limbs(input2, length2);
   
      while ((1<<log_length2) < length2) log_length2++; 
      
      /* Start with an upper bound on the number of bits needed */
      output_bits = FLINT_BITS * (size1 + size2) + log_length2 + 1; 
	} else
		output_bits = FLINT_ABS(bits_in);
   
	/* round up for sqrt2 trick */
   output_bits = (((output_bits - 1) >> (log_length - 2)) + 1) << (log_length - 2);
   
   limbs = (output_bits - 1) / FLINT_BITS + 1; /* initial size of FFT coeffs */
   size = limbs + 1;
   n = (1L<<(log_length - 2));

   ii = malloc((4*(n + n*size) + 5*size)*sizeof(mp_limb_t));
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
      ii[i] = ptr;
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;
   tt = s1 + size;
   
   if (input1 != input2)
   {
      jj = malloc(4*(n + n*size)*sizeof(mp_limb_t));
      for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
         jj[i] = ptr;
   } else jj = ii;
   
   /* put coefficients into FFT polys */
   bits1 = _fmpz_vec_get_fft(ii, input1, limbs, length1, bits_in);
   for (i = length1; i < 4*n; i++)
      mpn_zero(ii[i], limbs + 1);

   if (input1 != input2) 
   {
      bits2 = _fmpz_vec_get_fft(jj, input2, limbs, length2, bits_in);
      for (i = length2; i < 4*n; i++)
         mpn_zero(jj[i], limbs + 1);
   }
   else bits2 = bits1;
   
   if (!bits_in)
	{
		if (bits1 < 0L || bits2 < 0L) 
      {
         sign = 1;  
         bits1 = FLINT_ABS(bits1);
         bits2 = FLINT_ABS(bits2);
      }
   
      /* Recompute the number of limbs now that we know how large everything really is */
      output_bits = bits1 + bits2 + log_length2 + sign;
      
      /* round up output bits */
      output_bits = (((output_bits - 1) >> (log_length - 2)) + 1) << (log_length - 2);
      
      limbs = (output_bits - 1) / FLINT_BITS + 1;
	} else if (bits_in < 0L) sign = 1;
           
   fft_convolution(ii, jj, log_length - 2, limbs, len_out, &t1, &t2, &s1, tt); 
    
   _fmpz_vec_set_fft(output, len_out, ii, limbs, sign); /* write output */
   
   free(ii); 
   if (input1 != input2) 
      free(jj);
}

void
fmpz_poly_mul_SS(fmpz_poly_t res,
                 const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    long rlen;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (len1 == 1 || len2 == 1)
    {
        fmpz_poly_mul_classical(res, poly1, poly2);
        return;
    }

    rlen = len1 + len2 - 1;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        fmpz_poly_mul_SS(t, poly1, poly2);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    fmpz_poly_fit_length(res, rlen);
    if (len1 >= len2)
        _fmpz_poly_mul_SS(res->coeffs, poly1->coeffs, len1,
                          poly2->coeffs, len2, 0);
    else
        _fmpz_poly_mul_SS(res->coeffs, poly2->coeffs, len2,
                          poly1->coeffs, len1, 0);
    _fmpz_poly_set_length(res, rlen);
}
