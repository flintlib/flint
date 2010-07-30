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

=============================================================================*/
/******************************************************************************

 Copyright (C) 2010 Sebastian Pancratz
 Copyright (C) 2010 William Hart
 
******************************************************************************/

#ifndef FMPQ_POLY_H
#define FMPQ_POLY_H

#include <mpir.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*  Type definitions  ********************************************************/

typedef struct
{
    fmpz * coeffs;
    fmpz_t den;
    long alloc;
    long length;
} fmpq_poly_struct;

typedef fmpq_poly_struct fmpq_poly_t[1];

/*  Memory management  *******************************************************/

void fmpq_poly_init(fmpq_poly_t poly);

void fmpq_poly_init2(fmpq_poly_t poly, long alloc);

void fmpq_poly_realloc(fmpq_poly_t poly, long alloc);

void fmpq_poly_fit_length(fmpq_poly_t poly, long len);

void _fmpq_poly_set_length(fmpq_poly_t poly, long len);

void fmpq_poly_clear(fmpq_poly_t poly);

void _fmpq_poly_normalise(fmpq_poly_t poly);

void _fmpq_poly_canonicalise(fmpz * rpoly, fmpz_t den, long len);

void fmpq_poly_canonicalise(fmpq_poly_t poly);

/*  Accessing numerator and denominator  *************************************/

#define fmpq_poly_numref(poly)  ((poly)->coeffs)

#define fmpq_poly_denref(poly)  ((poly)->den)

/*  Polynomial parameters  ***************************************************/

static inline 
long fmpq_poly_degree(fmpq_poly_t poly)
{
    return poly->length - 1;
}

static inline 
long fmpq_poly_length(fmpq_poly_t poly)
{
    return poly->length;
}

/*  Randomisation  ***********************************************************/

void fmpq_poly_randinit(void);

void fmpq_poly_randclear(void);

void fmpq_poly_randtest(fmpq_poly_t f, long len, mp_bitcnt_t bits_in);

void fmpq_poly_randtest_unsigned(fmpq_poly_t f, long len, mp_bitcnt_t bits_in);

void fmpq_poly_randtest_not_zero(fmpq_poly_t f, long len, mp_bitcnt_t bits_in);

/*  Assignment and basic manipulation  ***************************************/

void fmpq_poly_set(fmpq_poly_t poly1, const fmpq_poly_t poly2);

void fmpq_poly_set_si(fmpq_poly_t poly, long x);

void fmpq_poly_set_ui(fmpq_poly_t poly, ulong x);

void fmpq_poly_set_fmpz(fmpq_poly_t poly, const fmpz_t x);

void fmpq_poly_set_mpz(fmpq_poly_t poly, const mpz_t x);

void fmpq_poly_set_mpq(fmpq_poly_t poly, const mpq_t x);

void fmpq_poly_zero(fmpq_poly_t poly);

void fmpq_poly_neg(fmpq_poly_t poly1, const fmpq_poly_t poly2);

void fmpq_poly_inv(fmpq_poly_t poly1, const fmpq_poly_t poly2);

void fmpq_poly_swap(fmpq_poly_t poly1, fmpq_poly_t poly2);

/*  Getting and setting coefficients  ****************************************/

void fmpq_poly_get_coeff_mpq(mpq_t x, const fmpq_poly_t poly, long n);

void fmpq_poly_set_coeff_si(fmpq_poly_t poly, long n, long x);

void fmpq_poly_set_coeff_ui(fmpq_poly_t poly, long n, ulong x);

void fmpq_poly_set_coeff_fmpz(fmpq_poly_t poly, long n, const fmpz_t x);

void fmpq_poly_set_coeff_mpz(fmpq_poly_t poly, long n, const mpz_t x);

void fmpq_poly_set_coeff_mpq(fmpq_poly_t poly, long n, const mpq_t x);

/*  Comparison  **************************************************************/

int fmpq_poly_equal(const fmpq_poly_t poly1, const fmpq_poly_t poly2);

int fmpq_poly_cmp(const fmpq_poly_t left, const fmpq_poly_t right);

static inline 
int fmpq_poly_is_zero(const fmpq_poly_t poly)
{
    return poly->length == 0;
}

static inline 
int fmpq_poly_is_one(const fmpq_poly_t poly)
{
    return (poly->length == 1) && (fmpz_equal(poly->coeffs, poly->den));
}

/*  Addition and subtraction  ************************************************/

void _fmpq_poly_add(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, long len1,
                    const fmpz * poly2, const fmpz_t den2, long len2);

void fmpq_poly_add(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2);

void _fmpq_poly_sub(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, long len1,
                    const fmpz * poly2, const fmpz_t den2, long len2);

void fmpq_poly_sub(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2);

/*  Scalar multiplication and division  **************************************/

void _fmpq_poly_scalar_mul_si(fmpz * rpoly, fmpz_t rden, 
                       const fmpz * poly, const fmpz_t den, long len, long c);

void _fmpq_poly_scalar_mul_ui(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, long len, ulong c);

void _fmpq_poly_scalar_mul_fmpz(fmpz * rpoly, fmpz_t rden, 
               const fmpz * poly, const fmpz_t den, long len, const fmpz_t c);

void _fmpq_poly_scalar_mul_mpq(fmpz * rpoly, fmpz_t rden, const fmpz * poly, 
                  const fmpz_t den, long len, const fmpz_t r, const fmpz_t s);

void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, long c);

void fmpq_poly_scalar_mul_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c);

void fmpq_poly_scalar_mul_fmpz(fmpq_poly_t rop, 
                               const fmpq_poly_t op, const fmpz_t c);

void fmpq_poly_scalar_mul_mpz(fmpq_poly_t rop, 
                              const fmpq_poly_t op, const mpz_t c);

void fmpq_poly_scalar_mul_mpq(fmpq_poly_t rop, 
                              const fmpq_poly_t op, const mpq_t c);

void _fmpq_poly_scalar_div_si(fmpz * rpoly, fmpz_t rden, 
                       const fmpz * poly, const fmpz_t den, long len, long c);

void _fmpq_poly_scalar_div_ui(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, long len, ulong c);

void _fmpq_poly_scalar_div_fmpz(fmpz * rpoly, fmpz_t rden, 
               const fmpz * poly, const fmpz_t den, long len, const fmpz_t c);

void _fmpq_poly_scalar_div_mpq(fmpz * rpoly, fmpz_t rden, const fmpz * poly, 
                  const fmpz_t den, long len, const fmpz_t r, const fmpz_t s);

void fmpq_poly_scalar_div_si(fmpq_poly_t rop, const fmpq_poly_t op, long c);

void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c);

void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, 
                               const fmpq_poly_t op, const fmpz_t c);

void fmpq_poly_scalar_div_mpz(fmpq_poly_t rop, 
                              const fmpq_poly_t op, const mpz_t c);

void fmpq_poly_scalar_div_mpq(fmpq_poly_t rop, 
                              const fmpq_poly_t op, const mpq_t c);

/*  Multiplication  **********************************************************/

void _fmpq_poly_mul(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, long len1, 
                    const fmpz * poly2, const fmpz_t den2, long len2);

void fmpq_poly_mul(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2);

/*  Powering  ****************************************************************/

void _fmpq_poly_pow(fmpz * rpoly, fmpz_t rden, const fmpz * poly, 
                                  const fmpz_t den, long len, long e);

void fmpq_poly_pow(fmpq_poly_t rpoly, const fmpq_poly_t poly, long e);

/*  Shifting  ****************************************************************/

void fmpq_poly_shift_left(fmpq_poly_t res, const fmpq_poly_t poly, long n);

void fmpq_poly_shift_right(fmpq_poly_t res, const fmpq_poly_t poly, long n);

/*  Euclidean division  ******************************************************/

void _fmpq_poly_divrem(fmpz * Q, fmpz_t q, fmpz * R, fmpz_t r,
                       const fmpz * A, const fmpz_t a, long lenA,
                       const fmpz * B, const fmpz_t b, long lenB);

void fmpq_poly_divrem(fmpq_poly_t Q, fmpq_poly_t R,
                      const fmpq_poly_t poly1, const fmpq_poly_t poly2);

/*  Derivative  **************************************************************/

void _fmpq_poly_derivative(fmpz * rpoly, fmpz_t rden, 
                           const fmpz * poly, const fmpz_t den, long len);

void fmpq_poly_derivative(fmpq_poly_t res, const fmpq_poly_t poly);

/*  Evaluation  **************************************************************/

void _fmpq_poly_evaluate_fmpz(mpq_t res, const fmpz * poly, 
                             const fmpz_t den, long len, const fmpz_t a);

void fmpq_poly_evaluate_fmpz(mpq_t res, const fmpq_poly_t poly, 
                             const fmpz_t a);

void _fmpq_poly_evaluate_mpq(mpq_t res, const fmpz * poly, 
                             const fmpz_t den, long len, const mpq_t a);

void fmpq_poly_evaluate_mpq(mpq_t res, const fmpq_poly_t poly, const mpq_t a);

/*  Gaussian content  ********************************************************/

void _fmpq_poly_content(mpq_t res, 
                        const fmpz * poly, const fmpz_t den, long len);

void fmpq_poly_content(mpq_t res, const fmpq_poly_t poly);

void _fmpq_poly_primitive_part(fmpz * rpoly, fmpz_t rden, 
                               const fmpz * poly, const fmpz_t den, long len);

void fmpq_poly_primitive_part(fmpq_poly_t res, const fmpq_poly_t poly);

int _fmpq_poly_is_monic(const fmpz * poly, const fmpz_t den, long len);

int fmpq_poly_is_monic(const fmpq_poly_t poly);

void _fmpq_poly_monic(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, long len);

void fmpq_poly_monic(fmpq_poly_t res, const fmpq_poly_t poly);

/*  Square-free  *************************************************************/

int _fmpq_poly_is_squarefree(const fmpz * poly, const fmpz_t den, long len);

int fmpq_poly_is_squarefree(const fmpq_poly_t poly);

/*  String conversion  *******************************************************/

char * fmpq_poly_to_string(const fmpq_poly_t poly);

char * fmpq_poly_to_string_pretty(const fmpq_poly_t poly, const char * var);

void fmpq_poly_print(const fmpq_poly_t poly);

void fmpq_poly_print_pretty(const fmpq_poly_t poly, const char * var);

#endif

