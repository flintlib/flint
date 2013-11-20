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

    Copyright (C) 2006, 2007, 2008, 2009, 2010, 2013 William Hart
    Copyright (C) 2009, 2011 Andy Novocin
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#ifndef FMPZ_POLY_H
#define FMPZ_POLY_H

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "nmod_poly.h"


#ifdef __cplusplus
 extern "C" {
#endif

#define FMPZ_POLY_INV_NEWTON_CUTOFF  32

/*  Type definitions *********************************************************/

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
} fmpz_poly_struct;

typedef fmpz_poly_struct fmpz_poly_t[1];

typedef struct
{
   fmpz ** powers;
   slong len;
} fmpz_poly_powers_precomp_struct;

typedef fmpz_poly_powers_precomp_struct fmpz_poly_powers_precomp_t[1];

typedef struct {
    fmpz c;
    fmpz_poly_struct *p;
    slong *exp;
    slong num;
    slong alloc;
} fmpz_poly_factor_struct;

typedef fmpz_poly_factor_struct fmpz_poly_factor_t[1];

/*  Memory management ********************************************************/

void fmpz_poly_init(fmpz_poly_t poly);

void fmpz_poly_init2(fmpz_poly_t poly, slong alloc);

void fmpz_poly_realloc(fmpz_poly_t poly, slong alloc);

void fmpz_poly_fit_length(fmpz_poly_t poly, slong len);

void fmpz_poly_clear(fmpz_poly_t poly);

void _fmpz_poly_normalise(fmpz_poly_t poly);

static __inline__
void _fmpz_poly_set_length(fmpz_poly_t poly, slong newlen)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = newlen;
}

/*  Polynomial parameters  ***************************************************/

static __inline__
slong fmpz_poly_length(const fmpz_poly_t poly)
{
    return poly->length;
}

static __inline__
slong fmpz_poly_degree(const fmpz_poly_t poly)
{
    return poly->length - 1;
}

/*  Assignment and basic manipulation  ***************************************/

void fmpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_set_ui(fmpz_poly_t poly, ulong c);

void fmpz_poly_set_si(fmpz_poly_t poly, slong c);

void fmpz_poly_set_fmpz(fmpz_poly_t poly, const fmpz_t c);

void fmpz_poly_set_mpz(fmpz_poly_t poly, const mpz_t c);

int _fmpz_poly_set_str(fmpz * poly, const char * str);

int fmpz_poly_set_str(fmpz_poly_t poly, const char * str);

char * _fmpz_poly_get_str(const fmpz * poly, slong len);

char * fmpz_poly_get_str(const fmpz_poly_t poly);

char * _fmpz_poly_get_str_pretty(const fmpz * poly, slong len, const char * x);

char * fmpz_poly_get_str_pretty(const fmpz_poly_t poly, const char * x);

static __inline__
void fmpz_poly_zero(fmpz_poly_t poly)
{
   _fmpz_poly_set_length(poly, 0);
}

static __inline__
void fmpz_poly_one(fmpz_poly_t poly)
{
    fmpz_poly_set_ui(poly, UWORD(1));
}

void fmpz_poly_zero_coeffs(fmpz_poly_t poly, slong i, slong j);

void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2);

void _fmpz_poly_reverse(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_poly_reverse(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

static __inline__
void fmpz_poly_truncate(fmpz_poly_t poly, slong newlen)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = newlen;
        _fmpz_poly_normalise(poly);
    }  
}

/*  Randomisation  ***********************************************************/

void fmpz_poly_randtest(fmpz_poly_t f, flint_rand_t state, 
                                                slong len, mp_bitcnt_t bits);

void fmpz_poly_randtest_unsigned(fmpz_poly_t f, flint_rand_t state, 
                                                slong len, mp_bitcnt_t bits);

void fmpz_poly_randtest_not_zero(fmpz_poly_t f, flint_rand_t state,
                                                slong len, mp_bitcnt_t bits);

/*  Getting and setting coefficients  ****************************************/

slong fmpz_poly_get_coeff_si(const fmpz_poly_t poly, slong n);

void fmpz_poly_set_coeff_si(fmpz_poly_t poly, slong n, slong x);

ulong fmpz_poly_get_coeff_ui(const fmpz_poly_t poly, slong n);

void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, slong n, ulong x);

void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, slong n, const fmpz_t x);

void fmpz_poly_get_coeff_fmpz(fmpz_t x, const fmpz_poly_t poly, slong n);

#define fmpz_poly_get_coeff_ptr(poly, n) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

#define fmpz_poly_lead(poly) \
    ((poly)->length ? (poly)->coeffs + (poly)->length - 1 : NULL)

/*  Comparison  **************************************************************/

int fmpz_poly_equal(const fmpz_poly_t poly1, const fmpz_poly_t poly2);

#define fmpz_poly_is_zero(poly) \
    ((poly)->length == 0)

static __inline__ 
int _fmpz_poly_is_one(const fmpz *poly, slong len)
{
    return (len > 0 && fmpz_is_one(poly) 
                    && _fmpz_vec_is_zero(poly + 1, len - 1));
}

static __inline__
int fmpz_poly_is_one(const fmpz_poly_t op)
{
    return (op->length) == 1 && (*(op->coeffs) == WORD(1));
}

static __inline__
int fmpz_poly_is_unit(const fmpz_poly_t op)
{
    return (op->length == 1) && (*(op->coeffs) == WORD(1) || *(op->coeffs) == WORD(-1));
}

static __inline__
int fmpz_poly_equal_fmpz(const fmpz_poly_t poly, const fmpz_t c)
{
	return ((poly->length == 0) && fmpz_is_zero(c)) ||
        ((poly->length == 1) && fmpz_equal(poly->coeffs, c));
}

/*  Addition and subtraction  ************************************************/

void _fmpz_poly_add(fmpz * res, const fmpz * poly1, slong len1, 
                                             const fmpz * poly2, slong len2);

void fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                   const fmpz_poly_t poly2);

void _fmpz_poly_sub(fmpz * res, const fmpz * poly1, slong len1, 
                                             const fmpz * poly2, slong len2);

void fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                   const fmpz_poly_t poly2);

void fmpz_poly_neg(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Scalar multiplication and division  **************************************/

void fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1, 
                             const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_mul_si(fmpz_poly_t poly1, 
                             const fmpz_poly_t poly2, slong x);

void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1, 
                               const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1, 
                                   const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_submul_fmpz(fmpz_poly_t poly1, 
                                   const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_fdiv_ui(fmpz_poly_t poly1, 
                              const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_fdiv_si(fmpz_poly_t poly1, 
                              const fmpz_poly_t poly2, slong x);

void fmpz_poly_scalar_fdiv_fmpz(fmpz_poly_t poly1, 
                                const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_tdiv_ui(fmpz_poly_t poly1, 
                              const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_tdiv_si(fmpz_poly_t poly1, 
                              const fmpz_poly_t poly2, slong x);

void fmpz_poly_scalar_tdiv_fmpz(fmpz_poly_t poly1, 
                                const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_divexact_ui(fmpz_poly_t poly1, 
                                  const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_divexact_si(fmpz_poly_t poly1, 
                                  const fmpz_poly_t poly2, slong x);

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
                                slong len, mp_bitcnt_t bit_size, int negate);

int _fmpz_poly_bit_unpack(fmpz * poly, slong len, 
                           mp_srcptr arr, mp_bitcnt_t bit_size, int negate);

void _fmpz_poly_bit_unpack_unsigned(fmpz * poly, slong len, 
                                       mp_srcptr arr, mp_bitcnt_t bit_size);

void fmpz_poly_bit_pack(fmpz_t f, const fmpz_poly_t poly,
        mp_bitcnt_t bit_size);

void fmpz_poly_bit_unpack(fmpz_poly_t poly, const fmpz_t f,
        mp_bitcnt_t bit_size);

void fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly, const fmpz_t f,
        mp_bitcnt_t bit_size);


/*  Multiplication  **********************************************************/

void _fmpz_poly_mul_classical(fmpz * res, const fmpz * poly1, slong len1, 
                                             const fmpz * poly2, slong len2);

void fmpz_poly_mul_classical(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_classical(fmpz * res, const fmpz * poly1, slong len1, 
                                     const fmpz * poly2, slong len2, slong n);

void fmpz_poly_mullow_classical(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                           const fmpz_poly_t poly2, slong n);

void _fmpz_poly_mulhigh_classical(fmpz * res, const fmpz * poly1, 
                      slong len1, const fmpz * poly2, slong len2, slong start);

void fmpz_poly_mulhigh_classical(fmpz_poly_t res, 
              const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong start);

void _fmpz_poly_mulmid_classical(fmpz * res, const fmpz * poly1, 
                                  slong len1, const fmpz * poly2, slong len2);

void fmpz_poly_mulmid_classical(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_mul_karatsuba(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mul_karatsuba(fmpz * res, const fmpz * poly1, 
                                  slong len1, const fmpz * poly2, slong len2);

void _fmpz_poly_mullow_karatsuba_n(fmpz * res, const fmpz * poly1, 
                                                const fmpz * poly2, slong n);

void fmpz_poly_mullow_karatsuba_n(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void _fmpz_poly_mulhigh_karatsuba_n(fmpz * res, const fmpz * poly1, 
                                              const fmpz * poly2, slong len);

void fmpz_poly_mulhigh_karatsuba_n(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong length);

void _fmpz_poly_mul_KS(fmpz * res, const fmpz * poly1, slong len1, 
                                             const fmpz * poly2, slong len2);

void fmpz_poly_mul_KS(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_KS(fmpz * res, const fmpz * poly1, slong len1, 
                                     const fmpz * poly2, slong len2, slong n);

void fmpz_poly_mullow_KS(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                           const fmpz_poly_t poly2, slong n);

void _fmpz_poly_mul_SS(fmpz * output, const fmpz * input1, slong length1, 
                                         const fmpz * input2, slong length2);

void fmpz_poly_mul_SS(fmpz_poly_t res,
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_SS(fmpz * output, const fmpz * input1, slong length1, 
                                 const fmpz * input2, slong length2, slong n);

void fmpz_poly_mullow_SS(fmpz_poly_t res,
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void _fmpz_poly_mul(fmpz * res, const fmpz * poly1, 
                                  slong len1, const fmpz * poly2, slong len2);

void fmpz_poly_mul(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow(fmpz * res, const fmpz * poly1, slong len1, 
                                     const fmpz * poly2, slong len2, slong n);

void fmpz_poly_mullow(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void fmpz_poly_mulhigh_n(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

/* Squaring ******************************************************************/

void _fmpz_poly_sqr_KS(fmpz * rop, const fmpz * op, slong len);

void fmpz_poly_sqr_KS(fmpz_poly_t rop, const fmpz_poly_t op);

void fmpz_poly_sqr_karatsuba(fmpz_poly_t rop, const fmpz_poly_t op);

void _fmpz_poly_sqr_karatsuba(fmpz * rop, const fmpz * op, slong len);

void _fmpz_poly_sqr_classical(fmpz * rop, const fmpz * op, slong len);

void fmpz_poly_sqr_classical(fmpz_poly_t rop, const fmpz_poly_t op);

void _fmpz_poly_sqr(fmpz * rop, const fmpz * op, slong len);

void fmpz_poly_sqr(fmpz_poly_t rop, const fmpz_poly_t op);

void _fmpz_poly_sqrlow_KS(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_poly_sqrlow_KS(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

void _fmpz_poly_sqrlow_karatsuba_n(fmpz * res, const fmpz * poly, slong n);

void fmpz_poly_sqrlow_karatsuba_n(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

void _fmpz_poly_sqrlow_classical(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_poly_sqrlow_classical(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

void _fmpz_poly_sqrlow(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_poly_sqrlow(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

/*  Powering  ****************************************************************/

void _fmpz_poly_pow_multinomial(fmpz * res, const fmpz * poly, slong len, ulong e);

void fmpz_poly_pow_multinomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_binomial(fmpz * res, const fmpz * poly, ulong e);

void fmpz_poly_pow_binomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_binexp(fmpz * res, const fmpz * poly, slong len, ulong e);

void fmpz_poly_pow_binexp(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_addchains(fmpz * res, const fmpz * poly, slong len, const int * a, int n);

void fmpz_poly_pow_addchains(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_small(fmpz * res, const fmpz * poly, slong len, ulong e);

void _fmpz_poly_pow(fmpz * res, const fmpz * poly, slong len, ulong e);

void fmpz_poly_pow(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_trunc(fmpz * res, const fmpz * poly, ulong e, slong n);

void 
fmpz_poly_pow_trunc(fmpz_poly_t res, const fmpz_poly_t poly, ulong e, slong n);

/*  Shifting  ****************************************************************/

void _fmpz_poly_shift_left(fmpz * res, const fmpz * poly, slong len, slong n);

void _fmpz_poly_shift_right(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_poly_shift_left(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

void fmpz_poly_shift_right(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

/*  Norms  *******************************************************************/

void _fmpz_poly_2norm(fmpz_t res, const fmpz * poly, slong len);

void fmpz_poly_2norm(fmpz_t res, const fmpz_poly_t poly);

mp_bitcnt_t _fmpz_poly_2norm_normalised_bits(const fmpz * poly, slong len);

static __inline__ 
ulong fmpz_poly_max_limbs(const fmpz_poly_t poly)
{
    return _fmpz_vec_max_limbs(poly->coeffs, poly->length);
}

static __inline__ 
slong fmpz_poly_max_bits(const fmpz_poly_t poly)
{
    return _fmpz_vec_max_bits(poly->coeffs, poly->length);
}

static __inline__ void
fmpz_poly_height(fmpz_t res, const fmpz_poly_t poly)
{
    _fmpz_vec_height(res, poly->coeffs, poly->length);
}

/*  Greatest common divisor  *************************************************/

void _fmpz_poly_gcd_subresultant(fmpz * res, const fmpz * poly1, slong len1, 
                                              const fmpz * poly2, slong len2);

void fmpz_poly_gcd_subresultant(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

int _fmpz_poly_gcd_heuristic(fmpz * res, const fmpz * poly1, slong len1, 
                                              const fmpz * poly2, slong len2);

int fmpz_poly_gcd_heuristic(fmpz_poly_t res,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_gcd_modular(fmpz * res, const fmpz * poly1, slong len1, 
                                              const fmpz * poly2, slong len2);

void fmpz_poly_gcd_modular(fmpz_poly_t res,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_gcd(fmpz * res, const fmpz * poly1, slong len1, 
                                              const fmpz * poly2, slong len2);

void fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

void _fmpz_poly_lcm(fmpz * res, const fmpz * poly1, slong len1, 
                                              const fmpz * poly2, slong len2);

void fmpz_poly_lcm(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

void _fmpz_poly_resultant(fmpz_t res, const fmpz * poly1, slong len1, 
                                              const fmpz * poly2, slong len2);

void fmpz_poly_resultant(fmpz_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

void _fmpz_poly_xgcd_modular(fmpz_t r, fmpz * s, fmpz * t, 
               const fmpz * poly1, slong len1, const fmpz * poly2, slong len2);

void fmpz_poly_xgcd_modular(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2);

static __inline__
void _fmpz_poly_xgcd(fmpz_t r, fmpz * s, fmpz * t, 
                const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
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

void _fmpz_poly_content(fmpz_t res, const fmpz * poly, slong len);

void fmpz_poly_content(fmpz_t res, const fmpz_poly_t poly);

void _fmpz_poly_primitive_part(fmpz * res, const fmpz * poly, slong len);

void fmpz_poly_primitive_part(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Square-free  *************************************************************/

int _fmpz_poly_is_squarefree(const fmpz * poly, slong len);

int fmpz_poly_is_squarefree(const fmpz_poly_t poly);

/*  Euclidean division  ******************************************************/

void _fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A, 
                                       slong lenA, const fmpz * B, slong lenB);

void fmpz_poly_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                                   const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, fmpz * W, 
                                 const fmpz * A, const fmpz * B, slong lenB);

void _fmpz_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
                     const fmpz * A, slong lenA, const fmpz * B, slong lenB);

void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, 
                                   const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_divrem(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, 
                                           const fmpz * B, slong lenB);

void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, 
                                                          const fmpz_poly_t B);

void _fmpz_poly_div_basecase(fmpz * Q, fmpz * R, const fmpz * A, slong lenA,
                                                   const fmpz * B, slong lenB);

void fmpz_poly_div_basecase(fmpz_poly_t Q, 
                                     const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_divremlow_divconquer_recursive(fmpz * Q, fmpz * QB, 
                                   const fmpz * A, const fmpz * B, slong lenB);

void _fmpz_poly_div_divconquer_recursive(fmpz * Q, fmpz * temp, 
                                   const fmpz * A, const fmpz * B, slong lenB);

void _fmpz_poly_div_divconquer(fmpz * Q, const fmpz * A, slong lenA, 
                                                   const fmpz * B, slong lenB);

void fmpz_poly_div_divconquer(fmpz_poly_t Q, 
                                     const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_div(fmpz * Q, const fmpz * A, slong lenA, 
                                                   const fmpz * B, slong lenB);

void fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_preinvert(fmpz * B_inv, const fmpz * B, slong n);

void fmpz_poly_preinvert(fmpz_poly_t B_inv, const fmpz_poly_t B);

void _fmpz_poly_div_preinv(fmpz * Q, const fmpz * A, slong len1, 
                               const fmpz * B, const fmpz * B_inv, slong len2);

void fmpz_poly_div_preinv(fmpz_poly_t Q, const fmpz_poly_t A, 
                                 const fmpz_poly_t B, const fmpz_poly_t B_inv);

void _fmpz_poly_divrem_preinv(fmpz * Q, fmpz * A, slong len1, 
                               const fmpz * B, const fmpz * B_inv, slong len2);

void fmpz_poly_divrem_preinv(fmpz_poly_t Q, fmpz_poly_t R, 
            const fmpz_poly_t A, const fmpz_poly_t B, const fmpz_poly_t B_inv);

fmpz ** _fmpz_poly_powers_precompute(const fmpz * B, slong len);

void fmpz_poly_powers_precompute(fmpz_poly_powers_precomp_t pinv, 
                                                             fmpz_poly_t poly);

void _fmpz_poly_powers_clear(fmpz ** powers, slong len);

void fmpz_poly_powers_clear(fmpz_poly_powers_precomp_t pinv);

void _fmpz_poly_rem_powers_precomp(fmpz * A, slong m, 
                                const fmpz * B, slong n, fmpz ** const powers);

void fmpz_poly_rem_powers_precomp(fmpz_poly_t R, 
                             const fmpz_poly_t A, const fmpz_poly_t B, 
                                       const fmpz_poly_powers_precomp_t B_inv);

void _fmpz_poly_rem_basecase(fmpz * Q, const fmpz * A, slong lenA,
                                                   const fmpz * B, slong lenB);

void fmpz_poly_rem_basecase(fmpz_poly_t R, 
                                     const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_rem(fmpz * R, const fmpz * A, slong lenA, 
                              const fmpz * B, slong lenB);

void fmpz_poly_rem(fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void
fmpz_poly_div_root(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_t c);

void
_fmpz_poly_div_root(fmpz * Q, const fmpz * A, slong len, const fmpz_t c);

/*  Power series division  ***************************************************/

void _fmpz_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, slong n);

void fmpz_poly_inv_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n);

static __inline__ void 
_fmpz_poly_inv_series(fmpz * Qinv, const fmpz * Q, slong n)
{
    _fmpz_poly_inv_series_newton(Qinv, Q, n);
}

static __inline__ void 
fmpz_poly_inv_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
{
    fmpz_poly_inv_series_newton(Qinv, Q, n);
}

void _fmpz_poly_div_series(fmpz * Q, const fmpz * A, const fmpz * B, slong n);

void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A, 
                                         const fmpz_poly_t B, slong n);

/*  Divisibility testing  ***************************************************/

int _fmpz_poly_divides(fmpz * q, const fmpz * a, 
                                         slong len1, const fmpz * b, slong len2);

int fmpz_poly_divides(fmpz_poly_t q, const fmpz_poly_t a, const fmpz_poly_t b);


/*  Pseudo division  *********************************************************/

void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, 
                   ulong * d, const fmpz * A, slong A_len, 
                        const fmpz * B, slong B_len, const fmpz_preinvn_t inv);

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                          ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_divrem_divconquer(fmpz * Q, fmpz * R, 
                    ulong * d, const fmpz * A, slong lenA, 
                         const fmpz * B, slong lenB, const fmpz_preinvn_t inv);

void fmpz_poly_pseudo_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, 
                          ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_divrem_cohen(fmpz * Q, fmpz * R, const fmpz * A, 
                                       slong lenA, const fmpz * B, slong lenB);

void fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R, 
                                     const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_rem_cohen(fmpz * R, const fmpz * A, slong lenA, 
                                                   const fmpz * B, slong lenB);

void fmpz_poly_pseudo_rem_cohen(fmpz_poly_t R, const fmpz_poly_t A, 
                                                          const fmpz_poly_t B);

static __inline__
void _fmpz_poly_pseudo_divrem(fmpz * Q, fmpz * R, 
                    ulong * d, const fmpz * A, slong A_len, 
                         const fmpz * B, slong B_len, const fmpz_preinvn_t inv)
{
    _fmpz_poly_pseudo_divrem_divconquer(Q, R, d, A, A_len, B, B_len, inv);
}

static __inline__
void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R, 
                           ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
    fmpz_poly_pseudo_divrem_divconquer(Q, R, d, A, B);
}

void _fmpz_poly_pseudo_div(fmpz * Q, ulong * d, const fmpz * A, slong lenA, 
                         const fmpz * B, slong lenB, const fmpz_preinvn_t inv);

void fmpz_poly_pseudo_div(fmpz_poly_t Q, ulong * d, const fmpz_poly_t A, 
                                                          const fmpz_poly_t B);

void _fmpz_poly_pseudo_rem(fmpz * R, ulong * d, const fmpz * A, slong lenA, 
                         const fmpz * B, slong lenB, const fmpz_preinvn_t inv);

void fmpz_poly_pseudo_rem(fmpz_poly_t R, ulong * d, const fmpz_poly_t A, 
                                                          const fmpz_poly_t B);

/*  Derivative  **************************************************************/

void _fmpz_poly_derivative(fmpz * rpoly, const fmpz * poly, slong len);
 
void fmpz_poly_derivative(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Evaluation  **************************************************************/

void 
_fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz * poly, slong len, 
                                                const fmpz_t a);

void fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz_poly_t poly, 
                                        const fmpz_t a);

void _fmpz_poly_evaluate_horner_fmpz(fmpz_t res, const fmpz * f, slong len, 
                                                               const fmpz_t a);

void fmpz_poly_evaluate_horner_fmpz(fmpz_t res, const fmpz_poly_t f, 
                                                               const fmpz_t a);

void _fmpz_poly_evaluate_fmpz(fmpz_t res, const fmpz * f, slong len, const fmpz_t a);

void fmpz_poly_evaluate_fmpz(fmpz_t res, const fmpz_poly_t f, const fmpz_t a);

void _fmpz_poly_evaluate_horner_mpq(fmpz_t rnum, fmpz_t rden, 
                                    const fmpz * f, slong len, 
                                    const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_horner_mpq(mpq_t res, const fmpz_poly_t f, 
                                                                const mpq_t a);

void _fmpz_poly_evaluate_mpq(fmpz_t rnum, fmpz_t rden,
                             const fmpz * f, slong len, 
                             const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_mpq(mpq_t res, const fmpz_poly_t f, const mpq_t a);

mp_limb_t _fmpz_poly_evaluate_mod(const fmpz * poly, slong len, mp_limb_t a, 
                                  mp_limb_t n, mp_limb_t ninv);

mp_limb_t fmpz_poly_evaluate_mod(const fmpz_poly_t poly, mp_limb_t a, 
                                 mp_limb_t n);

void 
_fmpz_poly_evaluate_divconquer(fmpz * res, const fmpz * poly, slong len, 
                               const fmpz_t x);

void 
fmpz_poly_evaluate_divconquer(fmpz_t res, 
                              const fmpz_poly_t poly, const fmpz_t x);

/*  Composition  *************************************************************/

void _fmpz_poly_compose_horner(fmpz * res, const fmpz * poly1, slong len1, 
                                                const fmpz * poly2, slong len2);

void fmpz_poly_compose_horner(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_compose_divconquer(fmpz * res, const fmpz * poly1, slong len1, 
                                                const fmpz * poly2, slong len2);

void fmpz_poly_compose_divconquer(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_compose(fmpz * res, const fmpz * poly1, slong len1, 
                                                const fmpz * poly2, slong len2);

void fmpz_poly_compose(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

/*  Taylor shift  ************************************************************/

void _fmpz_poly_taylor_shift_horner(fmpz * poly, const fmpz_t c, slong n);

void fmpz_poly_taylor_shift_horner(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c);

void _fmpz_poly_taylor_shift_divconquer(fmpz * poly, const fmpz_t c, slong n);

void fmpz_poly_taylor_shift_divconquer(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c);

void _fmpz_poly_taylor_shift(fmpz * poly, const fmpz_t c, slong n);

void fmpz_poly_taylor_shift(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c);

/*  Power series composition and compositional inverse  **********************/

void
_fmpz_poly_compose_series_brent_kung(fmpz * res, const fmpz * poly1, slong len1, 
                                      const fmpz * poly2, slong len2, slong n);

void
fmpz_poly_compose_series_brent_kung(fmpz_poly_t res, 
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void
_fmpz_poly_compose_series_horner(fmpz * res, const fmpz * poly1, slong len1, 
                                      const fmpz * poly2, slong len2, slong n);

void
fmpz_poly_compose_series_horner(fmpz_poly_t res, 
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void
_fmpz_poly_compose_series(fmpz * res, const fmpz * poly1, slong len1, 
                                      const fmpz * poly2, slong len2, slong n);

void
fmpz_poly_compose_series(fmpz_poly_t res, 
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void
_fmpz_poly_revert_series_lagrange(fmpz * Qinv, const fmpz * Q, slong n);

void
fmpz_poly_revert_series_lagrange(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n);

void
_fmpz_poly_revert_series_lagrange_fast(fmpz * Qinv, const fmpz * Q, slong n);

void
fmpz_poly_revert_series_lagrange_fast(fmpz_poly_t Qinv,
    const fmpz_poly_t Q, slong n);

void
_fmpz_poly_revert_series_newton(fmpz * Qinv, const fmpz * Q, slong n);

void
fmpz_poly_revert_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n);

void
_fmpz_poly_revert_series(fmpz * Qinv, const fmpz * Q, slong n);

void
fmpz_poly_revert_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n);

/*  Square root  *************************************************************/

int _fmpz_poly_sqrt_classical(fmpz * res, const fmpz * poly, slong len);

int fmpz_poly_sqrt_classical(fmpz_poly_t b, const fmpz_poly_t a);

int _fmpz_poly_sqrt(fmpz * res, const fmpz * poly, slong len);

int fmpz_poly_sqrt(fmpz_poly_t b, const fmpz_poly_t a);


/*  Signature  ***************************************************************/

void _fmpz_poly_signature(slong * r1, slong * r2, const fmpz * poly, slong len);

void fmpz_poly_signature(slong * r1, slong * r2, const fmpz_poly_t poly);

/*  Input and output  ********************************************************/

int fmpz_poly_fprint(FILE * file, const fmpz_poly_t poly);

int _fmpz_poly_fprint_pretty(FILE * file, 
                             const fmpz * poly, slong len, const char * x);

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
    flint_printf("(alloc = %wd, length = %wd, vec = ", poly->alloc, poly->length);
    if (poly->coeffs)
    {
        flint_printf("{");
        _fmpz_vec_print(poly->coeffs, poly->alloc);
        flint_printf("}");
    }
    else
    {
        flint_printf("NULL");
    }
    flint_printf(")");
    fflush(stdout);
}

/*  CRT  ********************************************************************/

void fmpz_poly_get_nmod_poly(nmod_poly_t res, const fmpz_poly_t poly);

void fmpz_poly_set_nmod_poly(fmpz_poly_t res, const nmod_poly_t poly);

void fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, const nmod_poly_t poly);

void
_fmpz_poly_CRT_ui_precomp(fmpz * res, const fmpz * poly1, slong len1,
               const fmpz_t m1, mp_srcptr poly2, slong len2, mp_limb_t m2,
                mp_limb_t m2inv, fmpz_t m1m2, mp_limb_t c, int sign);

void _fmpz_poly_CRT_ui(fmpz * res, const fmpz * poly1, slong len1,
               const fmpz_t m1, mp_srcptr poly2, slong len2, mp_limb_t m2,
                                                    mp_limb_t m2inv, int sign);

void fmpz_poly_CRT_ui(fmpz_poly_t res, const fmpz_poly_t poly1,
                                     const fmpz_t m1, const nmod_poly_t poly2,
                                        int sign);


/* Products *****************************************************************/

void _fmpz_poly_product_roots_fmpz_vec(fmpz * poly,
                                        const fmpz * xs, slong n);

void fmpz_poly_product_roots_fmpz_vec(fmpz_poly_t poly,
                                        const fmpz * xs, slong n);

/* Newton basis *************************************************************/

void _fmpz_poly_monomial_to_newton(fmpz * poly, const fmpz * roots, slong n);

void _fmpz_poly_newton_to_monomial(fmpz * poly, const fmpz * roots, slong n);


/* Multipoint evaluation and interpolation *********************************/

void
fmpz_poly_evaluate_fmpz_vec(fmpz * res, const fmpz_poly_t f,
                                const fmpz * a, slong n);

void
fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n);

/* Hensel lifting ************************************************************/

void fmpz_poly_hensel_build_tree(slong * link, fmpz_poly_t *v, fmpz_poly_t *w, 
                                 const nmod_poly_factor_t fac);

void fmpz_poly_hensel_lift(fmpz_poly_t Gout, fmpz_poly_t Hout, 
    fmpz_poly_t Aout, fmpz_poly_t Bout, 
    const fmpz_poly_t f, 
    const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1);

void _fmpz_poly_hensel_lift_without_inverse(fmpz *G, fmpz *H, 
    const fmpz *f, slong lenF, 
    const fmpz *g, slong lenG, const fmpz *h, slong lenH, 
    const fmpz *a, slong lenA, const fmpz *b, slong lenB, 
    const fmpz_t p, const fmpz_t p1);

void fmpz_poly_hensel_lift_without_inverse(fmpz_poly_t Gout, fmpz_poly_t Hout, 
    const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_poly_t h, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1);

void _fmpz_poly_hensel_lift_only_inverse(fmpz *A, fmpz *B, 
    const fmpz *G, slong lenG, const fmpz *H, slong lenH, 
    const fmpz *a, slong lenA, const fmpz *b, slong lenB, 
    const fmpz_t p, const fmpz_t p1);

void fmpz_poly_hensel_lift_only_inverse(fmpz_poly_t Aout, fmpz_poly_t Bout, 
    const fmpz_poly_t G, const fmpz_poly_t H, 
    const fmpz_poly_t a, const fmpz_poly_t b, 
    const fmpz_t p, const fmpz_t p1);

void fmpz_poly_hensel_lift_tree_recursive(slong *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, slong j, slong inv, 
    const fmpz_t p0, const fmpz_t p1);

void fmpz_poly_hensel_lift_tree(slong *link, fmpz_poly_t *v, fmpz_poly_t *w, 
    fmpz_poly_t f, slong r, const fmpz_t p, slong e0, slong e1, slong inv);

slong _fmpz_poly_hensel_start_lift(fmpz_poly_factor_t lifted_fac, slong *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, const fmpz_poly_t f, 
    const nmod_poly_factor_t local_fac, slong target_exp);

slong _fmpz_poly_hensel_continue_lift(fmpz_poly_factor_t lifted_fac, 
    slong *link, fmpz_poly_t *v, fmpz_poly_t *w, const fmpz_poly_t f, 
    slong prev, slong curr, slong N, const fmpz_t p);

void fmpz_poly_hensel_lift_once(fmpz_poly_factor_t lifted_fac, 
                                const fmpz_poly_t f, 
                                const nmod_poly_factor_t local_fac, slong N);

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

static __inline__ void fmpz_poly_set_coeff_mpz(fmpz_poly_t poly, slong n,
    const mpz_t x)
{
    fmpz_t t;
    fmpz_init_set_readonly(t, x);
    fmpz_poly_set_coeff_fmpz(poly, n, t);
    fmpz_clear_readonly(t);
}

static __inline__ void fmpz_poly_get_coeff_mpz(mpz_t x, const fmpz_poly_t poly, slong n)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_poly_get_coeff_fmpz(t, poly, n);
    fmpz_get_mpz(x, t);
    fmpz_clear(t);
}

/* Roots */

void _fmpz_poly_bound_roots(fmpz_t bound, const fmpz * poly, slong len);

void fmpz_poly_bound_roots(fmpz_t bound, const fmpz_poly_t poly);

#ifdef __cplusplus
}
#endif

#include "fmpz_poly_factor.h"

#endif

