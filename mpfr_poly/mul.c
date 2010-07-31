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

#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include <math.h>
#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_poly.h"

/*
   Remove any terms resulting from the convolution which are
   below the precision bound, i.e. with exponent < 2*s - prec
   + filter_bits.
*/
void _mpfr_poly_filter(mpfr * poly, long len, double s, mpfr_prec_t prec, mpfr_prec_t filter_bits)
{
   long i;
   long cutoff = ((long)(s/log(2)) - prec) + filter_bits;
   
   for (i = 0; i < len; i++)
   {
	  if (mpfr_get_exp(poly + i) <= cutoff)
         mpfr_set_ui(poly + i, 0, GMP_RNDN);
   }
}

/*
   We multiply poly1 of length len1 by poly2 of length len2 and store 
   the result in poly1. We assume that poly1 and poly2 have space for a 
   full convolution of power of 2 length.

   Inputs are destroyed.

   We assume that all coefficients p1_i, p2_i satisfy 
     s - prec <= log|p_i| < s and len1, len2 > 0
*/
void _mpfr_poly_mul_inplace(mpfr * poly1, long len1, 
							mpfr * poly2, long len2, mpfr_prec_t prec, double s, mpfr_prec_t fb)
{
   if (len1*prec < MUL_INPLACE_CUTOFF || len2*prec < MUL_INPLACE_CUTOFF) // use classical
   {
      mpfr * temp = _mpfr_vec_init(len1, prec);
      
	  _mpfr_vec_copy(temp, poly1, len1);
      
	  if (len1 >= len2)
		 _mpfr_poly_mul_classical(poly1, temp, len1, poly2, len2, prec);
	  else
         _mpfr_poly_mul_classical(poly1, poly2, len2, temp, len1, prec);

	  _mpfr_poly_filter(poly1, len1 + len2 - 1, s, prec, fb);

	  _mpfr_vec_clear(temp, len1);
   } else // use FFT
   {
      long log_len = FLINT_BIT_COUNT(len1 + len2 - 1); // ceil(log_2())
	  long length = (1L << log_len);
      long i;

	  for (i = len1; i < length; i++)
	     mpfr_set_ui(poly1 + i, 0, GMP_RNDN);
	  for (i = len2; i < length; i++)
	     mpfr_set_ui(poly2 + i, 0, GMP_RNDN);

      _mpfr_poly_convolution_FHT(poly1, poly2, log_len, prec);

	  _mpfr_poly_filter(poly1, len1 + len2 - 1, s, prec, fb);
   }
}

/*
   We multiply the coefficients of poly by a scaling factor. Coefficient i
   is multiplied by shift*scale^i. If shift is NULL a value of 1 is used
   for shift.
*/
void _mpfr_poly_shift_scale(mpfr * res, mpfr * poly, long len, 
							        mpfr_t shift, mpfr_t scale, mpfr_prec_t prec)
{
   long i;
   mpfr_t mult;
   mpfr_init2(mult, prec);

   if (shift != NULL)
	   mpfr_set(mult, shift, GMP_RNDN);
   else
	   mpfr_set_ui(mult, 1, GMP_RNDN);

   for (i = 0; i + 1 < len; i++)
   {
      mpfr_mul(res + i, poly + i, mult, GMP_RNDN);
	  mpfr_mul(mult, mult, scale, GMP_RNDN);
   }
   mpfr_mul(res + i, poly + i, mult, GMP_RNDN);

   mpfr_clear(mult);
}

void _mpfr_poly_mul_scale(mpfr * res, mpfr * poly1, long len1, mpfr * poly2, long len2, 
						    double s1, double i1, double i2, long w, mpfr_prec_t prec, mpfr_prec_t fb)
{
	mpfr * t1, * t2;
	long len_out = len1 + len2 - 1;
    long log_len = 0;
	long length, length2;
	long i;
	long w1, w2;
	
	mpfr_t shift, scale;

	mpfr_init2(shift, prec);
	mpfr_init2(scale, prec);
	
	w1 = FLINT_MIN(w, len1);
	w2 = FLINT_MIN(w, len2);
	
	log_len = FLINT_BIT_COUNT(len1 + w2 - 1);
	length = (1L << log_len);

	t1 = _mpfr_vec_init(length, prec);
    t2 = _mpfr_vec_init(length, prec);
	
	mpfr_set_d(scale, -s1, GMP_RNDN);
	mpfr_exp(scale, scale, GMP_RNDN);
	
	_mpfr_poly_shift_scale(t1, poly1, len1, NULL, scale, prec);
	_mpfr_vec_zero(t1 + len1, length - len1);
	
	mpfr_set_d(shift, i1 - i2, GMP_RNDN);
	mpfr_exp(shift, shift, GMP_RNDN);

	_mpfr_poly_shift_scale(t2, poly2, w2, shift, scale, prec);
	_mpfr_vec_zero(t2 + w2, length - w2);
	
	_mpfr_poly_mul_inplace(t1, len1, t2, w2, prec, 2*i1, fb);

	mpfr_ui_div(scale, 1, scale, GMP_RNDN);
    mpfr_ui_div(shift, 1, shift, GMP_RNDN);
    
    _mpfr_poly_shift_scale(res, t1, len1 + w2 - 1, shift, scale, prec);
	_mpfr_vec_zero(res + len1 + w2 - 1, len_out - len1 - w2 + 1);
	
	for (i = w2; i < len2; i += w2)
	{
		long l2 = FLINT_MIN(w2, len2 - i);

		log_len = FLINT_BIT_COUNT(w1 + l2 - 1);
		length2 = (1L << log_len);

		_mpfr_vec_copy(t1, poly1 + len1 - w1, w1);
		_mpfr_vec_zero(t1 + w1, length2 - w1);
	
		_mpfr_vec_copy(t2, poly2 + i, l2);
		_mpfr_vec_zero(t2 + l2, length2 - l2);
	
        _mpfr_poly_mul_inplace(t1, w1, t2, l2, prec, i1 + i2, fb);
 
		_mpfr_vec_add(res + len1 - w1 + i, res + len1 - w1 + i, t1, w1 + l2 - 1);
	}

    _mpfr_vec_clear(t1, length);
    _mpfr_vec_clear(t2, length);

	mpfr_clear(shift);
	mpfr_clear(scale);
}

void mpfr_poly_mul(mpfr_poly_t res, mpfr_poly_t poly1, mpfr_poly_t poly2, mpfr_prec_t fb)
{
	long len1 = poly1->length;
	long len2 = poly2->length;
	long len_out = len1 + len2 - 1;
    mpfr_prec_t prec = poly1->prec;
	mpfr * p1, * p2;
	double slope1, slope2, inter1, inter2;
	long w;
    
	if (len1 == 0 || len2 == 0)
	{
		res->length = 0;
		return;
	}

	p1 = poly1->coeffs;
	p2 = poly2->coeffs;

	mpfr_poly_fit_length(res, len_out);
	
	_mpfr_poly_bound_newton(&inter1, &slope1, p1, len1, prec);
	_mpfr_poly_bound_newton(&inter2, &slope2, p2, len2, prec);
	w = (long)(((double)prec)/(slope1 - slope2));

	if (w >= 0.0)
		_mpfr_poly_mul_scale(res->coeffs, p1, len1, p2, len2, 
		                     slope1, inter1, inter2, w, prec, fb); 
    else
		_mpfr_poly_mul_scale(res->coeffs, p2, len2, p1, len1, 
		                     slope2, inter2, inter1, -w, prec, fb); 
	res->length = len_out;
}
