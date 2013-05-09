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

    Copyright (C) 2006, 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2009, 2011 Andy Novocin
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#ifndef FMPZ_POLY_H
#define FMPZ_POLY_H

#undef ulong /* interferes with system includes */
#include <stdio.h>
#define ulong unsigned long

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "nmod_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef struct
{
    fmpz * coeffs;
    len_t alloc;
    len_t length;
} fmpz_poly_struct;

typedef fmpz_poly_struct fmpz_poly_t[1];

typedef struct {
    fmpz c;
    fmpz_poly_struct *p;
    len_t *exp;
    len_t num;
    len_t alloc;
} fmpz_poly_factor_struct;

typedef fmpz_poly_factor_struct fmpz_poly_factor_t[1];

/*  Memory management ********************************************************/

void fmpz_poly_init(fmpz_poly_t poly);

void fmpz_poly_init2(fmpz_poly_t poly, len_t alloc);

void fmpz_poly_realloc(fmpz_poly_t poly, len_t alloc);

void fmpz_poly_fit_length(fmpz_poly_t poly, len_t len);

void fmpz_poly_clear(fmpz_poly_t poly);

void _fmpz_poly_normalise(fmpz_poly_t poly);

static __inline__
void _fmpz_poly_set_length(fmpz_poly_t poly, len_t newlen)
{
    if (poly->length > newlen)
    {
        len_t i;
        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = newlen;
}

/*  Polynomial parameters  ***************************************************/

static __inline__
len_t fmpz_poly_length(const fmpz_poly_t poly)
{
    return poly->length;
}

static __inline__
len_t fmpz_poly_degree(const fmpz_poly_t poly)
{
    return poly->length - 1;
}

/*  Assignment and basic manipulation  ***************************************/

void fmpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_set_ui(fmpz_poly_t poly, ulong c);

void fmpz_poly_set_si(fmpz_poly_t poly, len_t c);

void fmpz_poly_set_fmpz(fmpz_poly_t poly, const fmpz_t c);

void fmpz_poly_set_mpz(fmpz_poly_t poly, const mpz_t c);

int _fmpz_poly_set_str(fmpz * poly, const char * str);

int fmpz_poly_set_str(fmpz_poly_t poly, const char * str);

char * _fmpz_poly_get_str(const fmpz * poly, len_t len);

char * fmpz_poly_get_str(const fmpz_poly_t poly);

char * _fmpz_poly_get_str_pretty(const fmpz * poly, len_t len, const char * x);

char * fmpz_poly_get_str_pretty(const fmpz_poly_t poly, const char * x);

static __inline__
void fmpz_poly_zero(fmpz_poly_t poly)
{
   _fmpz_poly_set_length(poly, 0);
}

static __inline__
void fmpz_poly_one(fmpz_poly_t poly)
{
    fmpz_poly_set_ui(poly, 1UL);
}

void fmpz_poly_zero_coeffs(fmpz_poly_t poly, len_t i, len_t j);

void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2);

void _fmpz_poly_reverse(fmpz * res, const fmpz * poly, len_t len, len_t n);

void fmpz_poly_reverse(fmpz_poly_t res, const fmpz_poly_t poly, len_t n);

static __inline__
void fmpz_poly_truncate(fmpz_poly_t poly, len_t newlen)
{
    if (poly->length > newlen)
    {
        len_t i;
        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = newlen;
        _fmpz_poly_normalise(poly);
    }  
}

/*  Randomisation  ***********************************************************/

void fmpz_poly_randtest(fmpz_poly_t f, flint_rand_t state, 
                                                len_t len, mp_bitcnt_t bits);

void fmpz_poly_randtest_unsigned(fmpz_poly_t f, flint_rand_t state, 
                                                len_t len, mp_bitcnt_t bits);

void fmpz_poly_randtest_not_zero(fmpz_poly_t f, flint_rand_t state,
                                                len_t len, mp_bitcnt_t bits);

/*  Getting and setting coefficients  ****************************************/

len_t fmpz_poly_get_coeff_si(const fmpz_poly_t poly, len_t n);

void fmpz_poly_set_coeff_si(fmpz_poly_t poly, len_t n, len_t x);

ulong fmpz_poly_get_coeff_ui(const fmpz_poly_t poly, len_t n);

void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, len_t n, ulong x);

void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, len_t n, const fmpz_t x);

void fmpz_poly_get_coeff_fmpz(fmpz_t x, const fmpz_poly_t poly, len_t n);

#define fmpz_poly_get_coeff_ptr(poly, n) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

#define fmpz_poly_lead(poly) \
    ((poly)->length ? (poly)->coeffs + (poly)->length - 1 : NULL)

/*  Comparison  **************************************************************/

int fmpz_poly_equal(const fmpz_poly_t poly1, const fmpz_poly_t poly2);

#define fmpz_poly_is_zero(poly) \
    ((poly)->length == 0)

static __inline__ 
int _fmpz_poly_is_one(const fmpz *poly, len_t len)
{
    return (len > 0 && fmpz_is_one(poly) 
                    && _fmpz_vec_is_zero(poly + 1, len - 1));
}

static __inline__
int fmpz_poly_is_one(const fmpz_poly_t op)
{
    return (op->length) == 1 && (*(op->coeffs) == 1L);
}

static __inline__
int fmpz_poly_is_unit(const fmpz_poly_t op)
{
    return (op->length == 1) && (*(op->coeffs) == 1L || *(op->coeffs) == -1L);
}

static __inline__
int fmpz_poly_equal_fmpz(const fmpz_poly_t poly, const fmpz_t c)
{
	return ((poly->length == 0) && fmpz_is_zero(c)) ||
        ((poly->length == 1) && fmpz_equal(poly->coeffs, c));
}

/*  Addition and subtraction  ************************************************/

void _fmpz_poly_add(fmpz * res, const fmpz * poly1, len_t len1, 
                                             const fmpz * poly2, len_t len2);

void fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                   const fmpz_poly_t poly2);

void _fmpz_poly_sub(fmpz * res, const fmpz * poly1, len_t len1, 
                                             const fmpz * poly2, len_t len2);

void fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                   const fmpz_poly_t poly2);

void fmpz_poly_neg(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Scalar multiplication and division  **************************************/

void fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1, 
                             const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_mul_si(fmpz_poly_t poly1, 
                             const fmpz_poly_t poly2, len_t x);

void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1, 
                               const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1, 
                                   const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_submul_fmpz(fmpz_poly_t poly1, 
                                   const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_fdiv_ui(fmpz_poly_t poly1, 
                              const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_fdiv_si(fmpz_poly_t poly1, 
                              const fmpz_poly_t poly2, len_t x);

void fmpz_poly_scalar_fdiv_fmpz(fmpz_poly_t poly1, 
                                const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_tdiv_ui(fmpz_poly_t poly1, 
                              const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_tdiv_si(fmpz_poly_t poly1, 
                              const fmpz_poly_t poly2, len_t x);

void fmpz_poly_scalar_tdiv_fmpz(fmpz_poly_t poly1, 
                                const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_divexact_ui(fmpz_poly_t poly1, 
                                  const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_divexact_si(fmpz_poly_t poly1, 
                                  const fmpz_poly_t poly2, len_t x);

void fmpz_poly_scalar_divexact_fmpz(fmpz_poly_t poly1, 
                                    const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_fdiv_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                           ulong exp);

void fmpz_poly_scalar_tdiv_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                           ulong exp);

void fmpz_poly_scalar_mul_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                           ulong exp);

static __inline__ 
void fmpz_poly_scalar_mod_fmpz(fmpz_poly_t poly1, 
                               const fmpz_poly_t poly2, const fmpz_t x)
{
    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
    }
    else
    {
        fmpz_poly_fit_length(poly1, poly2->length);
        _fmpz_vec_scalar_mod_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, poly2->length);
        _fmpz_poly_normalise(poly1);
    }
}

static __inline__ 
void fmpz_poly_scalar_smod_fmpz(fmpz_poly_t poly1, 
                                const fmpz_poly_t poly2, const fmpz_t x)
{
    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
    }
    else
    {
        fmpz_poly_fit_length(poly1, poly2->length);
        _fmpz_vec_scalar_smod_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, poly2->length);
        _fmpz_poly_normalise(poly1);

    }
}

/*  Bit packing  *************************************************************/

void _fmpz_poly_bit_pack(mp_ptr arr, const fmpz * poly,
                                len_t len, mp_bitcnt_t bit_size, int negate);

int _fmpz_poly_bit_unpack(fmpz * poly, len_t len, 
                           mp_srcptr arr, mp_bitcnt_t bit_size, int negate);

void _fmpz_poly_bit_unpack_unsigned(fmpz * poly, len_t len, 
                                       mp_srcptr arr, mp_bitcnt_t bit_size);

void fmpz_poly_bit_pack(fmpz_t f, const fmpz_poly_t poly,
        mp_bitcnt_t bit_size);

void fmpz_poly_bit_unpack(fmpz_poly_t poly, const fmpz_t f,
        mp_bitcnt_t bit_size);

void fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly, const fmpz_t f,
        mp_bitcnt_t bit_size);


/*  Multiplication  **********************************************************/

void _fmpz_poly_mul_classical(fmpz * res, const fmpz * poly1, len_t len1, 
                                             const fmpz * poly2, len_t len2);

void fmpz_poly_mul_classical(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_classical(fmpz * res, const fmpz * poly1, len_t len1, 
                                     const fmpz * poly2, len_t len2, len_t n);

void fmpz_poly_mullow_classical(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                           const fmpz_poly_t poly2, len_t n);

void _fmpz_poly_mulhigh_classical(fmpz * res, const fmpz * poly1, 
                      len_t len1, const fmpz * poly2, len_t len2, len_t start);

void fmpz_poly_mulhigh_classical(fmpz_poly_t res, 
              const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t start);

void _fmpz_poly_mulmid_classical(fmpz * res, const fmpz * poly1, 
                                  len_t len1, const fmpz * poly2, len_t len2);

void fmpz_poly_mulmid_classical(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_mul_karatsuba(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mul_karatsuba(fmpz * res, const fmpz * poly1, 
                                  len_t len1, const fmpz * poly2, len_t len2);

void _fmpz_poly_mullow_karatsuba_n(fmpz * res, const fmpz * poly1, 
                                                const fmpz * poly2, len_t n);

void fmpz_poly_mullow_karatsuba_n(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t n);

void _fmpz_poly_mulhigh_karatsuba_n(fmpz * res, const fmpz * poly1, 
                                              const fmpz * poly2, len_t len);

void fmpz_poly_mulhigh_karatsuba_n(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t length);

void _fmpz_poly_mul_KS(fmpz * res, const fmpz * poly1, len_t len1, 
                                             const fmpz * poly2, len_t len2);

void fmpz_poly_mul_KS(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_KS(fmpz * res, const fmpz * poly1, len_t len1, 
                                     const fmpz * poly2, len_t len2, len_t n);

void fmpz_poly_mullow_KS(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                           const fmpz_poly_t poly2, len_t n);

void _fmpz_poly_mul_SS(fmpz * output, const fmpz * input1, len_t length1, 
                                         const fmpz * input2, len_t length2);

void fmpz_poly_mul_SS(fmpz_poly_t res,
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_SS(fmpz * output, const fmpz * input1, len_t length1, 
                                 const fmpz * input2, len_t length2, len_t n);

void fmpz_poly_mullow_SS(fmpz_poly_t res,
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t n);

void _fmpz_poly_mul(fmpz * res, const fmpz * poly1, 
                                  len_t len1, const fmpz * poly2, len_t len2);

void fmpz_poly_mul(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow(fmpz * res, const fmpz * poly1, len_t len1, 
                                     const fmpz * poly2, len_t len2, len_t n);

void fmpz_poly_mullow(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t n);

void fmpz_poly_mulhigh_n(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t n);

/* Squaring ******************************************************************/

void _fmpz_poly_sqr_KS(fmpz * rop, const fmpz * op, len_t len);

void fmpz_poly_sqr_KS(fmpz_poly_t rop, const fmpz_poly_t op);

void fmpz_poly_sqr_karatsuba(fmpz_poly_t rop, const fmpz_poly_t op);

void _fmpz_poly_sqr_karatsuba(fmpz * rop, const fmpz * op, len_t len);

void _fmpz_poly_sqr_classical(fmpz * rop, const fmpz * op, len_t len);

void fmpz_poly_sqr_classical(fmpz_poly_t rop, const fmpz_poly_t op);

void _fmpz_poly_sqr(fmpz * rop, const fmpz * op, len_t len);

void fmpz_poly_sqr(fmpz_poly_t rop, const fmpz_poly_t op);

void _fmpz_poly_sqrlow_KS(fmpz * res, const fmpz * poly, len_t len, len_t n);

void fmpz_poly_sqrlow_KS(fmpz_poly_t res, const fmpz_poly_t poly, len_t n);

void _fmpz_poly_sqrlow_karatsuba_n(fmpz * res, const fmpz * poly, len_t n);

void fmpz_poly_sqrlow_karatsuba_n(fmpz_poly_t res, const fmpz_poly_t poly, len_t n);

void _fmpz_poly_sqrlow_classical(fmpz * res, const fmpz * poly, len_t len, len_t n);

void fmpz_poly_sqrlow_classical(fmpz_poly_t res, const fmpz_poly_t poly, len_t n);

void _fmpz_poly_sqrlow(fmpz * res, const fmpz * poly, len_t len, len_t n);

void fmpz_poly_sqrlow(fmpz_poly_t res, const fmpz_poly_t poly, len_t n);

/*  Powering  ****************************************************************/

void _fmpz_poly_pow_multinomial(fmpz * res, const fmpz * poly, len_t len, ulong e);

void fmpz_poly_pow_multinomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_binomial(fmpz * res, const fmpz * poly, ulong e);

void fmpz_poly_pow_binomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_binexp(fmpz * res, const fmpz * poly, len_t len, ulong e);

void fmpz_poly_pow_binexp(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_addchains(fmpz * res, const fmpz * poly, len_t len, const int * a, int n);

void fmpz_poly_pow_addchains(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_small(fmpz * res, const fmpz * poly, len_t len, ulong e);

void _fmpz_poly_pow(fmpz * res, const fmpz * poly, len_t len, ulong e);

void fmpz_poly_pow(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_trunc(fmpz * res, const fmpz * poly, ulong e, len_t n);

void 
fmpz_poly_pow_trunc(fmpz_poly_t res, const fmpz_poly_t poly, ulong e, len_t n);

/*  Shifting  ****************************************************************/

void _fmpz_poly_shift_left(fmpz * res, const fmpz * poly, len_t len, len_t n);

void _fmpz_poly_shift_right(fmpz * res, const fmpz * poly, len_t len, len_t n);

void fmpz_poly_shift_left(fmpz_poly_t res, const fmpz_poly_t poly, len_t n);

void fmpz_poly_shift_right(fmpz_poly_t res, const fmpz_poly_t poly, len_t n);

/*  Norms  *******************************************************************/

void _fmpz_poly_2norm(fmpz_t res, const fmpz * poly, len_t len);

void fmpz_poly_2norm(fmpz_t res, const fmpz_poly_t poly);

mp_bitcnt_t _fmpz_poly_2norm_normalised_bits(const fmpz * poly, len_t len);

static __inline__ 
ulong fmpz_poly_max_limbs(const fmpz_poly_t poly)
{
    return _fmpz_vec_max_limbs(poly->coeffs, poly->length);
}

static __inline__ 
len_t fmpz_poly_max_bits(const fmpz_poly_t poly)
{
    return _fmpz_vec_max_bits(poly->coeffs, poly->length);
}

static __inline__ void
fmpz_poly_height(fmpz_t res, const fmpz_poly_t poly)
{
    _fmpz_vec_height(res, poly->coeffs, poly->length);
}

/*  Greatest common divisor  *************************************************/

void _fmpz_poly_gcd_subresultant(fmpz * res, const fmpz * poly1, len_t len1, 
                                              const fmpz * poly2, len_t len2);

void fmpz_poly_gcd_subresultant(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

int _fmpz_poly_gcd_heuristic(fmpz * res, const fmpz * poly1, len_t len1, 
                                              const fmpz * poly2, len_t len2);

int fmpz_poly_gcd_heuristic(fmpz_poly_t res,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_gcd_modular(fmpz * res, const fmpz * poly1, len_t len1, 
                                              const fmpz * poly2, len_t len2);

void fmpz_poly_gcd_modular(fmpz_poly_t res,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_gcd(fmpz * res, const fmpz * poly1, len_t len1, 
                                              const fmpz * poly2, len_t len2);

void fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

void _fmpz_poly_lcm(fmpz * res, const fmpz * poly1, len_t len1, 
                                              const fmpz * poly2, len_t len2);

void fmpz_poly_lcm(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

void _fmpz_poly_resultant(fmpz_t res, const fmpz * poly1, len_t len1, 
                                              const fmpz * poly2, len_t len2);

void fmpz_poly_resultant(fmpz_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

void _fmpz_poly_xgcd_modular(fmpz_t r, fmpz * s, fmpz * t, 
               const fmpz * poly1, len_t len1, const fmpz * poly2, len_t len2);

void fmpz_poly_xgcd_modular(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2);

static __inline__
void _fmpz_poly_xgcd(fmpz_t r, fmpz * s, fmpz * t, 
                const fmpz * poly1, len_t len1, const fmpz * poly2, len_t len2)
{
    _fmpz_poly_xgcd_modular(r, s, t, poly1, len1, poly2, len2);
}

static __inline__
void fmpz_poly_xgcd(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t,
                            const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    fmpz_poly_xgcd_modular(r, s, t, poly1, poly2);
}

/*  Gaussian content  ********************************************************/

void _fmpz_poly_content(fmpz_t res, const fmpz * poly, len_t len);

void fmpz_poly_content(fmpz_t res, const fmpz_poly_t poly);

void _fmpz_poly_primitive_part(fmpz * res, const fmpz * poly, len_t len);

void fmpz_poly_primitive_part(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Euclidean division  ******************************************************/

void _fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A, 
                                       len_t lenA, const fmpz * B, len_t lenB);

void fmpz_poly_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                                   const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, fmpz * W, 
                                 const fmpz * A, const fmpz * B, len_t lenB);

void _fmpz_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
                     const fmpz * A, len_t lenA, const fmpz * B, len_t lenB);

void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, 
                                   const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_divrem(fmpz * Q, fmpz * R, const fmpz * A, len_t lenA, 
                                           const fmpz * B, len_t lenB);

void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, 
                                                    const fmpz_poly_t B);

void _fmpz_poly_div_basecase(fmpz * Q, fmpz * R, const fmpz * A, len_t lenA,
                                                  const fmpz * B, len_t lenB);

void fmpz_poly_div_basecase(fmpz_poly_t Q, 
                                   const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_divremlow_divconquer_recursive(fmpz * Q, fmpz * QB, 
                                          const fmpz * A, const fmpz * B, len_t lenB);

void _fmpz_poly_div_divconquer_recursive(fmpz * Q, fmpz * temp, 
                                  const fmpz * A, const fmpz * B, len_t lenB);

void _fmpz_poly_div_divconquer(fmpz * Q, const fmpz * A, len_t lenA, 
                                                  const fmpz * B, len_t lenB);

void fmpz_poly_div_divconquer(fmpz_poly_t Q, 
                                   const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_div(fmpz * Q, const fmpz * A, len_t lenA, 
                              const fmpz * B, len_t lenB);

void fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A, 
                                  const fmpz_poly_t B);

void _fmpz_poly_rem_basecase(fmpz * Q, const fmpz * A, len_t lenA,
                                       const fmpz * B, len_t lenB);

void fmpz_poly_rem_basecase(fmpz_poly_t R, 
                                   const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_rem(fmpz * R, const fmpz * A, len_t lenA, 
                              const fmpz * B, len_t lenB);

void fmpz_poly_rem(fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void
fmpz_poly_div_root(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_t c);

void
_fmpz_poly_div_root(fmpz * Q, const fmpz * A, len_t len, const fmpz_t c);

/*  Power series division  ***************************************************/

void _fmpz_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, len_t n);

void fmpz_poly_inv_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, len_t n);

static __inline__ void 
_fmpz_poly_inv_series(fmpz * Qinv, const fmpz * Q, len_t n)
{
    _fmpz_poly_inv_series_newton(Qinv, Q, n);
}

static __inline__ void 
fmpz_poly_inv_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, len_t n)
{
    fmpz_poly_inv_series_newton(Qinv, Q, n);
}

void _fmpz_poly_div_series(fmpz * Q, const fmpz * A, const fmpz * B, len_t n);

void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A, 
                                         const fmpz_poly_t B, len_t n);

/*  Divisibility testing  ***************************************************/

int _fmpz_poly_divides(fmpz * q, const fmpz * a, 
                                         len_t len1, const fmpz * b, len_t len2);

int fmpz_poly_divides(fmpz_poly_t q, const fmpz_poly_t a, const fmpz_poly_t b);


/*  Pseudo division  *********************************************************/

void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, ulong * d, 
                    const fmpz * A, len_t A_len, const fmpz * B, len_t B_len);

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                        ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_divrem_divconquer(fmpz * Q, fmpz * R, ulong * d, 
                       const fmpz * A, len_t lenA, const fmpz * B, len_t lenB);

void fmpz_poly_pseudo_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, 
                        ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_divrem_cohen(fmpz * Q, fmpz * R, const fmpz * A, 
                                       len_t lenA, const fmpz * B, len_t lenB);

void fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R, 
                                   const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_rem_cohen(fmpz * R, const fmpz * A, len_t lenA, 
                                                  const fmpz * B, len_t lenB);

void fmpz_poly_pseudo_rem_cohen(fmpz_poly_t R, const fmpz_poly_t A, 
                                                        const fmpz_poly_t B);

static __inline__
void _fmpz_poly_pseudo_divrem(fmpz * Q, fmpz * R, ulong * d, 
                    const fmpz * A, len_t A_len, const fmpz * B, len_t B_len)
{
    _fmpz_poly_pseudo_divrem_divconquer(Q, R, d, A, A_len, B, B_len);
}

static __inline__
void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R, 
                       ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
    fmpz_poly_pseudo_divrem_divconquer(Q, R, d, A, B);
}

void _fmpz_poly_pseudo_div(fmpz * Q, ulong * d, const fmpz * A, len_t lenA, 
                                                const fmpz * B, len_t lenB);

void fmpz_poly_pseudo_div(fmpz_poly_t Q, ulong * d, const fmpz_poly_t A, 
                                                    const fmpz_poly_t B);

void _fmpz_poly_pseudo_rem(fmpz * R, ulong * d, const fmpz * A, len_t lenA, 
                                                const fmpz * B, len_t lenB);

void fmpz_poly_pseudo_rem(fmpz_poly_t R, ulong * d, const fmpz_poly_t A, 
                                                    const fmpz_poly_t B);

/*  Derivative  **************************************************************/

void _fmpz_poly_derivative(fmpz * rpoly, const fmpz * poly, len_t len);
 
void fmpz_poly_derivative(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Evaluation  **************************************************************/

void 
_fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz * poly, len_t len, 
                                                const fmpz_t a);

void fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz_poly_t poly, 
                                        const fmpz_t a);

void _fmpz_poly_evaluate_horner_fmpz(fmpz_t res, const fmpz * f, len_t len, 
                                                               const fmpz_t a);

void fmpz_poly_evaluate_horner_fmpz(fmpz_t res, const fmpz_poly_t f, 
                                                               const fmpz_t a);

void _fmpz_poly_evaluate_fmpz(fmpz_t res, const fmpz * f, len_t len, const fmpz_t a);

void fmpz_poly_evaluate_fmpz(fmpz_t res, const fmpz_poly_t f, const fmpz_t a);

void _fmpz_poly_evaluate_horner_mpq(fmpz_t rnum, fmpz_t rden, 
                                    const fmpz * f, len_t len, 
                                    const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_horner_mpq(mpq_t res, const fmpz_poly_t f, 
                                                                const mpq_t a);

void _fmpz_poly_evaluate_mpq(fmpz_t rnum, fmpz_t rden,
                             const fmpz * f, len_t len, 
                             const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_mpq(mpq_t res, const fmpz_poly_t f, const mpq_t a);

mp_limb_t _fmpz_poly_evaluate_mod(const fmpz * poly, len_t len, mp_limb_t a, 
                                  mp_limb_t n, mp_limb_t ninv);

mp_limb_t fmpz_poly_evaluate_mod(const fmpz_poly_t poly, mp_limb_t a, 
                                 mp_limb_t n);

void 
_fmpz_poly_evaluate_divconquer(fmpz * res, const fmpz * poly, len_t len, 
                               const fmpz_t x);

void 
fmpz_poly_evaluate_divconquer(fmpz_t res, 
                              const fmpz_poly_t poly, const fmpz_t x);

/*  Composition  *************************************************************/

void _fmpz_poly_compose_horner(fmpz * res, const fmpz * poly1, len_t len1, 
                                                const fmpz * poly2, len_t len2);

void fmpz_poly_compose_horner(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_compose_divconquer(fmpz * res, const fmpz * poly1, len_t len1, 
                                                const fmpz * poly2, len_t len2);

void fmpz_poly_compose_divconquer(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_compose(fmpz * res, const fmpz * poly1, len_t len1, 
                                                const fmpz * poly2, len_t len2);

void fmpz_poly_compose(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

/*  Taylor shift  ************************************************************/

void _fmpz_poly_taylor_shift_horner(fmpz * poly, const fmpz_t c, len_t n);

void fmpz_poly_taylor_shift_horner(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c);

void _fmpz_poly_taylor_shift_divconquer(fmpz * poly, const fmpz_t c, len_t n);

void fmpz_poly_taylor_shift_divconquer(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c);

void _fmpz_poly_taylor_shift(fmpz * poly, const fmpz_t c, len_t n);

void fmpz_poly_taylor_shift(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c);

/*  Power series composition and compositional inverse  **********************/

void
_fmpz_poly_compose_series_brent_kung(fmpz * res, const fmpz * poly1, len_t len1, 
                                      const fmpz * poly2, len_t len2, len_t n);

void
fmpz_poly_compose_series_brent_kung(fmpz_poly_t res, 
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t n);

void
_fmpz_poly_compose_series_horner(fmpz * res, const fmpz * poly1, len_t len1, 
                                      const fmpz * poly2, len_t len2, len_t n);

void
fmpz_poly_compose_series_horner(fmpz_poly_t res, 
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t n);

void
_fmpz_poly_compose_series(fmpz * res, const fmpz * poly1, len_t len1, 
                                      const fmpz * poly2, len_t len2, len_t n);

void
fmpz_poly_compose_series(fmpz_poly_t res, 
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, len_t n);

void
_fmpz_poly_revert_series_lagrange(fmpz * Qinv, const fmpz * Q, len_t n);

void
fmpz_poly_revert_series_lagrange(fmpz_poly_t Qinv, const fmpz_poly_t Q, len_t n);

void
_fmpz_poly_revert_series_lagrange_fast(fmpz * Qinv, const fmpz * Q, len_t n);

void
fmpz_poly_revert_series_lagrange_fast(fmpz_poly_t Qinv,
    const fmpz_poly_t Q, len_t n);

void
_fmpz_poly_revert_series_newton(fmpz * Qinv, const fmpz * Q, len_t n);

void
fmpz_poly_revert_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, len_t n);

void
_fmpz_poly_revert_series(fmpz * Qinv, const fmpz * Q, len_t n);

void
fmpz_poly_revert_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, len_t n);

/*  Square root  *************************************************************/

int _fmpz_poly_sqrt_classical(fmpz * res, const fmpz * poly, len_t len);

int fmpz_poly_sqrt_classical(fmpz_poly_t b, const fmpz_poly_t a);

int _fmpz_poly_sqrt(fmpz * res, const fmpz * poly, len_t len);

int fmpz_poly_sqrt(fmpz_poly_t b, const fmpz_poly_t a);


/*  Signature  ***************************************************************/

void _fmpz_poly_signature(len_t * r1, len_t * r2, fmpz * poly, len_t len);

void fmpz_poly_signature(len_t * r1, len_t * r2, fmpz_poly_t poly);

/*  Input and output  ********************************************************/

int fmpz_poly_fprint(FILE * file, const fmpz_poly_t poly);

int _fmpz_poly_fprint_pretty(FILE * file, 
                             const fmpz * poly, len_t len, const char * x);

int fmpz_poly_fprint_pretty(FILE * file, 
                                       const fmpz_poly_t poly, const char * x);

static __inline__
int fmpz_poly_print(const fmpz_poly_t poly)
{
    return fmpz_poly_fprint(stdout, poly);
}

static __inline__
int fmpz_poly_print_pretty(const fmpz_poly_t poly, const char * x)
{
    return fmpz_poly_fprint_pretty(stdout, poly, x);
}

int fmpz_poly_fread(FILE * file, fmpz_poly_t poly);

int fmpz_poly_fread_pretty(FILE *file, fmpz_poly_t poly, char **x);

static __inline__ 
int fmpz_poly_read(fmpz_poly_t poly)
{
    return fmpz_poly_fread(stdin, poly);
}

static __inline__ 
int fmpz_poly_read_pretty(fmpz_poly_t poly, char **x)
{
    return fmpz_poly_fread_pretty(stdin, poly, x);
}

static __inline__
void fmpz_poly_debug(const fmpz_poly_t poly)
{
    printf("(alloc = %ld, length = %ld, vec = ", poly->alloc, poly->length);
    if (poly->coeffs)
    {
        printf("{");
        _fmpz_vec_print(poly->coeffs, poly->alloc);
        printf("}");
    }
    else
    {
        printf("NULL");
    }
    printf(")");
    fflush(stdout);
}

/*  CRT  ********************************************************************/

void fmpz_poly_get_nmod_poly(nmod_poly_t res, const fmpz_poly_t poly);

void fmpz_poly_set_nmod_poly(fmpz_poly_t res, const nmod_poly_t poly);

void fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, const nmod_poly_t poly);

void
_fmpz_poly_CRT_ui_precomp(fmpz * res, const fmpz * poly1, len_t len1,
               const fmpz_t m1, mp_srcptr poly2, len_t len2, mp_limb_t m2,
                mp_limb_t m2inv, fmpz_t m1m2, mp_limb_t c, int sign);

void _fmpz_poly_CRT_ui(fmpz * res, const fmpz * poly1, len_t len1,
               const fmpz_t m1, mp_srcptr poly2, len_t len2, mp_limb_t m2,
                                                    mp_limb_t m2inv, int sign);

void fmpz_poly_CRT_ui(fmpz_poly_t res, const fmpz_poly_t poly1,
                                     const fmpz_t m1, const nmod_poly_t poly2,
                                        int sign);


/* Products *****************************************************************/

void _fmpz_poly_product_roots_fmpz_vec(fmpz * poly,
                                        const fmpz * xs, len_t n);

void fmpz_poly_product_roots_fmpz_vec(fmpz_poly_t poly,
                                        const fmpz * xs, len_t n);

/* Newton basis *************************************************************/

void _fmpz_poly_monomial_to_newton(fmpz * poly, const fmpz * roots, len_t n);

void _fmpz_poly_newton_to_monomial(fmpz * poly, const fmpz * roots, len_t n);


/* Multipoint evaluation and interpolation *********************************/

void
fmpz_poly_evaluate_fmpz_vec(fmpz * res, const fmpz_poly_t f,
                                const fmpz * a, len_t n);

void
fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, len_t n);

/* Hensel lifting ************************************************************/

void fmpz_poly_hensel_build_tree(len_t * link, fmpz_poly_t *v, fmpz_poly_t *w, 
                                 const nmod_poly_factor_t fac);

void fmpz_poly_hensel_lift(fmpz_poly_t Gout, fmpz_poly_t Hout, 
    fmpz_poly_t Aout, fmpz_poly_t Bout, 
    const fmpz_poly_t f, 
    const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1);

void _fmpz_poly_hensel_lift_without_inverse(fmpz *G, fmpz *H, 
    const fmpz *f, len_t lenF, 
    const fmpz *g, len_t lenG, const fmpz *h, len_t lenH, 
    const fmpz *a, len_t lenA, const fmpz *b, len_t lenB, 
    const fmpz_t p, const fmpz_t p1);

void fmpz_poly_hensel_lift_without_inverse(fmpz_poly_t Gout, fmpz_poly_t Hout, 
    const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1);

void _fmpz_poly_hensel_lift_only_inverse(fmpz *A, fmpz *B, 
    const fmpz *G, len_t lenG, const fmpz *H, len_t lenH, 
    const fmpz *a, len_t lenA, const fmpz *b, len_t lenB, 
    const fmpz_t p, const fmpz_t p1);

void fmpz_poly_hensel_lift_only_inverse(fmpz_poly_t Aout, fmpz_poly_t Bout, 
    const fmpz_poly_t G, const fmpz_poly_t H, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1);

void fmpz_poly_hensel_lift_tree_recursive(len_t *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, len_t j, len_t inv, 
    const fmpz_t p0, const fmpz_t p1);

void fmpz_poly_hensel_lift_tree(len_t *link, fmpz_poly_t *v, fmpz_poly_t *w, 
    fmpz_poly_t f, len_t r, const fmpz_t p, len_t e0, len_t e1, len_t inv);

len_t _fmpz_poly_hensel_start_lift(fmpz_poly_factor_t lifted_fac, len_t *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, const fmpz_poly_t f, 
    const nmod_poly_factor_t local_fac, len_t target_exp);

len_t _fmpz_poly_hensel_continue_lift(fmpz_poly_factor_t lifted_fac, 
    len_t *link, fmpz_poly_t *v, fmpz_poly_t *w, const fmpz_poly_t f, 
    len_t prev, len_t curr, len_t N, const fmpz_t p);

void fmpz_poly_hensel_lift_once(fmpz_poly_factor_t lifted_fac, 
                                const fmpz_poly_t f, 
                                const nmod_poly_factor_t local_fac, len_t N);

/* Some functions for backwards compatibility */

static __inline__ void fmpz_poly_scalar_mul_mpz(fmpz_poly_t poly1,
                               const fmpz_poly_t poly2, const mpz_t x)
{
    fmpz_t t;
    fmpz_init_set_readonly(t, x);
    fmpz_poly_scalar_mul_fmpz(poly1, poly2, t);
    fmpz_clear_readonly(t);
}

static __inline__ void fmpz_poly_scalar_divexact_mpz(fmpz_poly_t poly1,
                               const fmpz_poly_t poly2, const mpz_t x)
{
    fmpz_t t;
    fmpz_init_set_readonly(t, x);
    fmpz_poly_scalar_divexact_fmpz(poly1, poly2, t);
    fmpz_clear_readonly(t);
}

static __inline__ void fmpz_poly_scalar_fdiv_mpz(fmpz_poly_t poly1,
                               const fmpz_poly_t poly2, const mpz_t x)
{
    fmpz_t t;
    fmpz_init_set_readonly(t, x);
    fmpz_poly_scalar_fdiv_fmpz(poly1, poly2, t);
    fmpz_clear_readonly(t);
}

static __inline__ void fmpz_poly_set_coeff_mpz(fmpz_poly_t poly, len_t n,
    const mpz_t x)
{
    fmpz_t t;
    fmpz_init_set_readonly(t, x);
    fmpz_poly_set_coeff_fmpz(poly, n, t);
    fmpz_clear_readonly(t);
}

static __inline__ void fmpz_poly_get_coeff_mpz(mpz_t x, fmpz_poly_t poly, len_t n)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_poly_get_coeff_fmpz(t, poly, n);
    fmpz_get_mpz(x, t);
    fmpz_clear(t);
}

/* Roots */

void _fmpz_poly_bound_roots(fmpz_t bound, const fmpz * poly, len_t len);

void fmpz_poly_bound_roots(fmpz_t bound, const fmpz_poly_t poly);

#ifdef __cplusplus
}
#endif

#endif

