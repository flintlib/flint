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

void mpfr_poly_mul_FHT(mpfr_poly_t res, mpfr_poly_t poly1, mpfr_poly_t poly2)
{
	long len1 = poly1->length;
	long len2 = poly2->length;
	long len_out = len1 + len2 - 1;
    mp_bitcnt_t prec = poly1->prec;
	mpfr * t1, * t2;
    long i;

	long log_length = 0;
	long length;

	if (len1 == 0 || len2 == 0)
	{
		res->length = 0;
		return;
	}

	while ((1L<<log_length) < len_out) log_length++;
	length = (1L<<log_length);

	mpfr_poly_fit_length(res, length);
	t1 = res->coeffs;
    t2 = _mpfr_vec_init(length, prec);
	
	_mpfr_vec_copy(t1, poly1->coeffs, poly1->length);
	for (i = poly1->length; i < length; i++)
		mpfr_set_ui(t1 + i, 0, GMP_RNDN);

    _mpfr_vec_copy(t2, poly2->coeffs, poly2->length);
	for (i = poly2->length; i < length; i++)
		mpfr_set_ui(t2 + i, 0, GMP_RNDN);

	_mpfr_poly_convolution_FHT(t1, t2, log_length, prec);
	res->length = len_out;

    _mpfr_vec_clear(t2, length);
}
