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
/******************************************************************************

 Copyright (C) 2010 William Hart
 
******************************************************************************/

#ifndef FMPZ_VEC_H
#define FMPZ_VEC_H

#include <mpir.h>
#include "fmpz.h"

fmpz * _fmpz_vec_init(ulong length);

void _fmpz_vec_clear(fmpz * in1, ulong length);

void _fmpz_vec_randinit(void);

void _fmpz_vec_randclear(void);

void _fmpz_vec_randtest(fmpz * f, ulong length, mp_bitcnt_t bits_in);

void _fmpz_vec_print(fmpz * vec, ulong length);

void _fmpz_vec_zero(fmpz * vec1, ulong len1);

void _fmpz_vec_copy(fmpz * vec1, const fmpz * vec2, ulong len2);

int _fmpz_vec_equal(const fmpz * vec1, const fmpz * vec2, ulong length);

void _fmpz_vec_neg(fmpz * vec1, const fmpz * vec2, ulong len2);

void _fmpz_vec_scalar_mul_si(fmpz * vec1, 
							 const fmpz * vec2, ulong len2, long c);

void _fmpz_vec_scalar_mul_fmpz(fmpz * vec1, 
		            const fmpz * vec2, ulong len2, const fmpz_t x);

void _fmpz_vec_scalar_addmul_si(fmpz * vec1, 
						     const fmpz * vec2, ulong len2, long c);

void _fmpz_vec_scalar_addmul_fmpz(fmpz * poly1, 
	                const fmpz * poly2, ulong len2, const fmpz_t x);

void _fmpz_vec_add(fmpz * res, const fmpz * vec1, 
				                     const fmpz * vec2, ulong len2);

void _fmpz_vec_sub(fmpz * res, const fmpz * vec1, 
				                     const fmpz * vec2, ulong len2);

void _fmpz_vec_scalar_addmul_si_2exp(fmpz * vec1, 
	              const fmpz * vec2, ulong len2, long c, ulong exp);

void _fmpz_vec_scalar_mul_2exp(fmpz * vec1, 
                          const fmpz * vec2, ulong len2, ulong exp);

#endif






