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
#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_poly.h"

void _mpfr_poly_mul_classical(mpfr * res, mpfr * in1, long len1, mpfr * in2, long len2, mp_bitcnt_t prec)
{
   mpfr * temp;
   long i;

   _mpfr_vec_scalar_mul_mpfr(res, in1, len1, in2);

   if (len2 == 1) return;
   
   _mpfr_vec_scalar_mul_mpfr(res + len1, in2 + 1, len2 - 1, in1 + len1 - 1);

   temp = _mpfr_vec_init(len1 - 1, prec);

   for (i = 1; i < len2; i++)
   {
      _mpfr_vec_scalar_mul_mpfr(temp, in1, len1 - 1, in2 + i);
	  _mpfr_vec_add(res + i, res + i, temp, len1 - 1);
   }

   _mpfr_vec_clear(temp, len1 - 1);
}

void mpfr_poly_mul_classical(mpfr_poly_t res, mpfr_poly_t poly1, mpfr_poly_t poly2)
{
	long len_out;
	mp_bitcnt_t prec = poly1->prec;

	if (poly1->length == 0 || poly2->length == 0)
	{
		res->length = 0;
		return;
	}

    len_out = poly1->length + poly2->length - 1;

	if (res == poly1 || res == poly2)
	{
		mpfr_poly_t temp;
		mpfr_poly_init2(temp, len_out, poly1->prec);
		
		if (poly1->length >= poly2->length)
			_mpfr_poly_mul_classical(temp->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, prec);
		else
			_mpfr_poly_mul_classical(temp->coeffs, poly2->coeffs, poly2->length, poly1->coeffs, poly1->length, prec);

        mpfr_poly_swap(temp, res);
		mpfr_poly_clear(temp);
	} else
	{
		mpfr_poly_fit_length(res, len_out);
		if (poly1->length >= poly2->length)
			_mpfr_poly_mul_classical(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, prec);
		else
			_mpfr_poly_mul_classical(res->coeffs, poly2->coeffs, poly2->length, poly1->coeffs, poly1->length, prec);
	}

	res->length = len_out;
}
