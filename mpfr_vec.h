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

#ifndef MFPR_VEC_H
#define MPFR_VEC_H

#include <gmp.h>
#include <mpfr.h> 

#ifdef __cplusplus
 extern "C" {
#endif

mpfr * _mpfr_vec_init(slong length, mp_bitcnt_t prec);

void _mpfr_vec_clear(mpfr * vec, slong length);

void _mpfr_vec_zero(mpfr * vec, slong length);

void _mpfr_vec_set(mpfr * vec1, const mpfr * vec2, slong length);

void _mpfr_vec_add(mpfr * res, const mpfr * vec1, const mpfr * vec2, slong length);

void _mpfr_vec_scalar_mul_2exp(mpfr * res, const mpfr * vec, slong length, mp_bitcnt_t exp);

void _mpfr_vec_scalar_mul_mpfr(mpfr * res, const mpfr * vec, slong length, mpfr_t c);

void _mpfr_vec_scalar_product(mpfr_t res, const mpfr * vec1, const mpfr * vec2, slong length);

#ifdef __cplusplus
}
#endif

#endif






