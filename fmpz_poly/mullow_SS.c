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
#include "fmpz_poly.h"
#include "fft.h"
#include "fft_tuning.h"

void _fmpz_poly_mullow_SS(fmpz * output, const fmpz * input1, len_t len1, 
               const fmpz * input2, len_t len2, len_t trunc)
{
    len_t len_out = len1 + len2 - 1;
    len_t loglen  = FLINT_CLOG2(len_out);
    len_t loglen2 = FLINT_CLOG2(len2);
    len_t n = (1L << (loglen - 2));

    len_t output_bits, limbs, size, i;
    mp_limb_t * ptr, * t1, * t2, * tt, * s1, ** ii, ** jj;
    len_t bits1, bits2;
    int sign = 0;

    ulong size1 = _fmpz_vec_max_limbs(input1, len1); 
    ulong size2 = _fmpz_vec_max_limbs(input2, len2);

    /* Start with an upper bound on the number of bits needed */
    output_bits = FLINT_BITS * (size1 + size2) + loglen2 + 1; 
    
    /* round up for sqrt2 trick */
    output_bits = (((output_bits - 1) >> (loglen - 2)) + 1) << (loglen - 2);

    limbs = (output_bits - 1) / FLINT_BITS + 1; /* initial size of FFT coeffs */
    if (limbs > FFT_MULMOD_2EXPP1_CUTOFF) /* can't be worse than next power of 2 limbs */
        limbs = (1L << FLINT_CLOG2(limbs));
    size = limbs + 1;

    /* allocate space for ffts */
    ii = flint_malloc((4*(n + n*size) + 5*size)*sizeof(mp_limb_t));
    for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
        ii[i] = ptr;
    t1 = ptr;
    t2 = t1 + size;
    s1 = t2 + size;
    tt = s1 + size;

    if (input1 != input2)
    {
        jj = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t));
        for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
            jj[i] = ptr;
    } else jj = ii;

    /* put coefficients into FFT vecs */
    bits1 = _fmpz_vec_get_fft(ii, input1, limbs, len1);
    for (i = len1; i < 4*n; i++)
        flint_mpn_zero(ii[i], limbs + 1);

    if (input1 != input2) 
    {
        bits2 = _fmpz_vec_get_fft(jj, input2, limbs, len2);
        for (i = len2; i < 4*n; i++)
            flint_mpn_zero(jj[i], limbs + 1);
    }
    else bits2 = bits1;

    if (bits1 < 0L || bits2 < 0L) 
    {
        sign = 1;  
        bits1 = FLINT_ABS(bits1);
        bits2 = FLINT_ABS(bits2);
    }

    /* Recompute the number of bits/limbs now that we know how large everything is */
    output_bits = bits1 + bits2 + loglen2 + sign;

    /* round up output bits for sqrt2 */
    output_bits = (((output_bits - 1) >> (loglen - 2)) + 1) << (loglen - 2);

    limbs = (output_bits - 1) / FLINT_BITS + 1;
    limbs = fft_adjust_limbs(limbs); /* round up limbs for Nussbaumer */
    
    fft_convolution(ii, jj, loglen - 2, limbs, len_out, &t1, &t2, &s1, tt); 

    _fmpz_vec_set_fft(output, trunc, ii, limbs, sign); /* write output */

    flint_free(ii); 
    if (input1 != input2) 
        flint_free(jj);
}

void
fmpz_poly_mullow_SS(fmpz_poly_t res,
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t n)
{
    const len_t len1 = poly1->length;
    const len_t len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (len1 == 1 || len2 == 1)
    {
        fmpz_poly_mullow_classical(res, poly1, poly2, n);
        return;
    }

    fmpz_poly_fit_length(res, n);

    if (len1 >= len2)
        _fmpz_poly_mullow_SS(res->coeffs, poly1->coeffs, len1,
                                          poly2->coeffs, len2, n);
    else
        _fmpz_poly_mullow_SS(res->coeffs, poly2->coeffs, len2,
                                          poly1->coeffs, len1, n);

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
