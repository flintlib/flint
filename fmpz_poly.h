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

 Copyright (C) 2006, 2007, 2008, 2009, 2010 William Hart
 Copyright (C) 2009, Andy Novocin
 
******************************************************************************/

#ifndef FMPZ_POLY_H
#define FMPZ_POLY_H

#include <mpir.h>
#include "fmpz.h"

typedef struct
{
   fmpz * coeffs;
   ulong alloc;
   ulong length;
} fmpz_poly_struct;

typedef fmpz_poly_struct fmpz_poly_t[1];

void fmpz_poly_init(fmpz_poly_t poly);

void fmpz_poly_init2(fmpz_poly_t poly, const ulong alloc);

void fmpz_poly_realloc(fmpz_poly_t poly, const ulong alloc);

void fmpz_poly_fit_length(fmpz_poly_t poly, const ulong length);

void fmpz_poly_clear(fmpz_poly_t poly);

void _fmpz_poly_normalise(fmpz_poly_t poly);

static inline
void _fmpz_poly_set_length(fmpz_poly_t poly, const ulong length)
{
	if (poly->length > length) // demote coefficients beyond new length
   {
      ulong i;
      for (i = length; i < poly->length; i++)
			_fmpz_demote(poly->coeffs + i);	
   } 

	poly->length = length;
}

static inline
void fmpz_poly_truncate(fmpz_poly_t poly, const ulong length)
{
	if (poly->length > length) // only truncate if necessary
   {
      ulong i;
      for (i = length; i < poly->length; i++)
			_fmpz_demote(poly->coeffs + i);
		poly->length = length;
      _fmpz_poly_normalise(poly);
   }  
}

static inline
void fmpz_poly_zero(fmpz_poly_t poly)
{
   _fmpz_poly_set_length(poly, 0);
}

void fmpz_poly_randinit(void);

void fmpz_poly_randclear(void);

void fmpz_poly_randtest(fmpz_poly_t f, ulong length, mp_bitcnt_t bits_in);

void fmpz_poly_randtest_unsigned(fmpz_poly_t f, ulong length, 
                                                     mp_bitcnt_t bits_in);

void fmpz_poly_randtest_not_zero(fmpz_poly_t f, 
                                       ulong length, mp_bitcnt_t bits_in);

long fmpz_poly_get_coeff_si(const fmpz_poly_t poly, const ulong n);

void fmpz_poly_set_coeff_si(fmpz_poly_t poly, ulong n, const long x);

ulong fmpz_poly_get_coeff_ui(const fmpz_poly_t poly, const ulong n);

void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, ulong n, const ulong x);

void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, ulong n, const fmpz_t x);

void fmpz_poly_get_coeff_fmpz(fmpz_t x, 
                                    const fmpz_poly_t poly, const ulong n);

void fmpz_poly_print(fmpz_poly_t poly);

void fmpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2);

int fmpz_poly_equal(const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2);

void fmpz_poly_neg(fmpz_poly_t res, const fmpz_poly_t poly);

void _fmpz_poly_add(fmpz * res, const fmpz * poly1, ulong len1, 
					                      const fmpz * poly2, ulong len2);

void fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                   const fmpz_poly_t poly2);

void _fmpz_poly_sub(fmpz * res, const fmpz * poly1, ulong len1, 
					                        const fmpz * poly2, ulong len2);

void fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                   const fmpz_poly_t poly2);

void fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1, fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_mul_si(fmpz_poly_t poly1, fmpz_poly_t poly2, long x);

void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1, 
                                   const fmpz_poly_t poly2, const fmpz_t x);

void _fmpz_poly_scalar_addmul_fmpz(fmpz * poly1, 
						    const fmpz * poly2, ulong len2, const fmpz_t x);

void fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1, 
								   const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_shift_left(fmpz_poly_t res, 
                                     const fmpz_poly_t poly, const ulong n);

void fmpz_poly_shift_right(fmpz_poly_t res, 
                                     const fmpz_poly_t poly, const ulong n);

void _fmpz_poly_mul_classical(fmpz * res, const fmpz * poly1, ulong len1, 
							                const fmpz * poly2, ulong len2);

void fmpz_poly_mul_classical(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_classical(fmpz * res, const fmpz * poly1,
				   ulong len1, const fmpz * poly2, ulong len2, ulong trunc);

void fmpz_poly_mullow_classical(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong trunc);

void _fmpz_poly_mulhigh_classical(fmpz * res, const fmpz * poly1, 
				   ulong len1, const fmpz * poly2, ulong len2, ulong start);

void fmpz_poly_mulhigh_classical(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong start);

void _fmpz_poly_mulmid_classical(fmpz * res, const fmpz * poly1, 
							    ulong len1, const fmpz * poly2, ulong len2);

void fmpz_poly_mulmid_classical(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_mul_karatsuba(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mul_karatsuba(fmpz * res, const fmpz * poly1, 
							    ulong len1, const fmpz * poly2, ulong len2);

void _fmpz_poly_mullow_karatsuba_n(fmpz * res, const fmpz * poly1, 
								             const fmpz * poly2, ulong len);

void fmpz_poly_mullow_karatsuba_n(fmpz_poly_t res, 
            const fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong length);

void _fmpz_poly_mulhigh_karatsuba_n(fmpz * res, const fmpz * poly1, 
									         const fmpz * poly2, ulong len);

void fmpz_poly_mulhigh_karatsuba_n(fmpz_poly_t res, 
            const fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong length);

void _fmpz_poly_bit_pack(mp_limb_t * arr, const fmpz * poly,
						             ulong len, ulong bit_size, int negate);

void _fmpz_poly_bit_unpack(fmpz * poly, ulong length, 
						 const mp_limb_t * arr, ulong bit_size, int negate);

void _fmpz_poly_bit_unpack_unsigned(fmpz * poly, ulong length, 
						             const mp_limb_t * arr, ulong bit_size);

void _fmpz_poly_mul_KS(fmpz * res, const fmpz * poly1, ulong len1, 
					                        const fmpz * poly2, ulong len2);

void fmpz_poly_mul_KS(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mul(fmpz * res, const fmpz * poly1, 
				                ulong len1, const fmpz * poly2, ulong len2);

void fmpz_poly_mul(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_mullow_n(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong trunc);

void fmpz_poly_mulhigh_n(fmpz_poly_t res, 
                 const fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong n);

void _fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A, 
								  ulong A_len, const fmpz * B, ulong B_len);

static inline
void _fmpz_poly_div_basecase(fmpz * Q, const fmpz * A, ulong A_len,
							                    const fmpz * B, ulong B_len)
{
   _fmpz_poly_divrem_basecase(Q, NULL, A, A_len, B, B_len);
}

void _fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, 
					           const fmpz * A, const fmpz * B, ulong B_len);

void _fmpz_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
				  const fmpz * A, ulong A_len, const fmpz * B, ulong B_len);

void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, 
								  const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, ulong * d, 
				  const fmpz * A, ulong A_len, const fmpz * B, ulong B_len);

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
					   ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, ulong * d, 
				  const fmpz * A, ulong A_len, const fmpz * B, ulong B_len);

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
					   ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

static inline 
void _fmpz_poly_pseudo_divrem(fmpz * Q, fmpz * R, ulong * d, 
                      const fmpz * A, ulong A_len, const fmpz * B, ulong B_len)
{
    _fmpz_poly_pseudo_divrem_basecase(Q, R, d, A, A_len, B, B_len);
}

static inline 
void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R, 
                           ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
    fmpz_poly_pseudo_divrem_basecase(Q, R, d, A, B);
}

#endif

