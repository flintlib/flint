/*
    Copyright (C) 2006, 2007, 2008, 2009, 2010, 2013 William Hart
    Copyright (C) 2009, 2011 Andy Novocin
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_POLY_H
#define FMPZ_POLY_H

#ifdef FMPZ_POLY_INLINES_C
#define FMPZ_POLY_INLINE
#else
#define FMPZ_POLY_INLINE static inline
#endif

#include "fmpz_types.h"
#include "nmod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FMPZ_POLY_INV_NEWTON_CUTOFF 32
#define FMPZ_POLY_SQRT_DIVCONQUER_CUTOFF 16
#define FMPZ_POLY_SQRTREM_DIVCONQUER_CUTOFF 16

/*  Type definitions *********************************************************/

typedef struct
{
   fmpz ** powers;
   slong len;
} fmpz_poly_powers_precomp_struct;

typedef fmpz_poly_powers_precomp_struct fmpz_poly_powers_precomp_t[1];

typedef struct
{
   ulong ** jj; /* used by fft_convolution_precache */
   slong n;
   slong len2;
   slong loglen;
   slong bits2;
   slong limbs;
   fmpz_poly_t poly2;
} fmpz_poly_mul_precache_struct;

typedef fmpz_poly_mul_precache_struct fmpz_poly_mul_precache_t[1];

/*  Memory management ********************************************************/

void fmpz_poly_init(fmpz_poly_t poly);

void fmpz_poly_init2(fmpz_poly_t poly, slong alloc);

void fmpz_poly_realloc(fmpz_poly_t poly, slong alloc);

void fmpz_poly_fit_length(fmpz_poly_t poly, slong len);

void fmpz_poly_clear(fmpz_poly_t poly);

void _fmpz_poly_normalise(fmpz_poly_t poly);

void _fmpz_poly_set_length(fmpz_poly_t poly, slong newlen);

FMPZ_POLY_INLINE
void fmpz_poly_attach_truncate(fmpz_poly_t trunc,
                                               const fmpz_poly_t poly, slong n)
{
   trunc->coeffs = poly->coeffs;
   trunc->alloc = poly->alloc;
   trunc->length = FLINT_MIN(poly->length, n);
}

FMPZ_POLY_INLINE
void fmpz_poly_attach_shift(fmpz_poly_t trunc, const fmpz_poly_t poly, slong n)
{
   trunc->coeffs = poly->coeffs + n;
   trunc->alloc = poly->alloc - n;
   trunc->length = FLINT_MAX(poly->length - n, 0);
}

/*  Polynomial parameters  ***************************************************/

FMPZ_POLY_INLINE
slong fmpz_poly_length(const fmpz_poly_t poly)
{
    return poly->length;
}

FMPZ_POLY_INLINE
slong fmpz_poly_degree(const fmpz_poly_t poly)
{
    return poly->length - 1;
}

/*  Assignment and basic manipulation  ***************************************/

void fmpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_set_ui(fmpz_poly_t poly, ulong c);

void fmpz_poly_set_si(fmpz_poly_t poly, slong c);

void fmpz_poly_set_fmpz(fmpz_poly_t poly, const fmpz_t c);

int _fmpz_poly_set_str(fmpz * poly, const char * str);

int fmpz_poly_set_str(fmpz_poly_t poly, const char * str);

char * _fmpz_poly_get_str(const fmpz * poly, slong len);

char * fmpz_poly_get_str(const fmpz_poly_t poly);

char * _fmpz_poly_get_str_pretty(const fmpz * poly, slong len, const char * x);

char * fmpz_poly_get_str_pretty(const fmpz_poly_t poly, const char * x);

FMPZ_POLY_INLINE
void fmpz_poly_zero(fmpz_poly_t poly)
{
   _fmpz_poly_set_length(poly, 0);
}

FMPZ_POLY_INLINE
void fmpz_poly_one(fmpz_poly_t poly)
{
    fmpz_poly_set_ui(poly, UWORD(1));
}

void fmpz_poly_zero_coeffs(fmpz_poly_t poly, slong i, slong j);

void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2);

void _fmpz_poly_reverse(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_poly_reverse(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

ulong _fmpz_poly_deflation(const fmpz* a, slong len);

FMPZ_POLY_INLINE
ulong fmpz_poly_deflation(const fmpz_poly_t input)
{
    return _fmpz_poly_deflation(input->coeffs, input->length);
}

void fmpz_poly_deflate(fmpz_poly_t result, const fmpz_poly_t input,
                                                              ulong deflation);

void fmpz_poly_inflate(fmpz_poly_t result, const fmpz_poly_t input,
                                                              ulong inflation);

void fmpz_poly_truncate(fmpz_poly_t poly, slong newlen);

void fmpz_poly_set_trunc(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

/*  Randomisation  ***********************************************************/

void fmpz_poly_randtest(fmpz_poly_t f, flint_rand_t state,
                                                slong len, flint_bitcnt_t bits);

void fmpz_poly_randtest_unsigned(fmpz_poly_t f, flint_rand_t state,
                                                slong len, flint_bitcnt_t bits);

void fmpz_poly_randtest_not_zero(fmpz_poly_t f, flint_rand_t state,
                                                slong len, flint_bitcnt_t bits);

void fmpz_poly_randtest_no_real_root(fmpz_poly_t p, flint_rand_t state,
                                                slong len, flint_bitcnt_t bits);

void fmpz_poly_randtest_irreducible1(fmpz_poly_t pol, flint_rand_t state, slong len, flint_bitcnt_t bits);
void fmpz_poly_randtest_irreducible2(fmpz_poly_t pol, flint_rand_t state, slong len, flint_bitcnt_t bits);
void fmpz_poly_randtest_irreducible(fmpz_poly_t pol, flint_rand_t state, slong len, flint_bitcnt_t bits);

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

int fmpz_poly_equal_trunc(const fmpz_poly_t poly1,
                                             const fmpz_poly_t poly2, slong n);

#define fmpz_poly_is_zero(poly) \
    ((poly)->length == 0)

int _fmpz_poly_is_one(const fmpz *poly, slong len);

FMPZ_POLY_INLINE
int fmpz_poly_is_one(const fmpz_poly_t op)
{
    return (op->length) == 1 && (*(op->coeffs) == WORD(1));
}

FMPZ_POLY_INLINE
int fmpz_poly_is_unit(const fmpz_poly_t op)
{
    return (op->length == 1) && (*(op->coeffs) == WORD(1) || *(op->coeffs) == WORD(-1));
}

FMPZ_POLY_INLINE
int fmpz_poly_is_gen(const fmpz_poly_t op)
{
    return (op->length) == 2 && (*(op->coeffs + 1) == WORD(1)) && (*(op->coeffs + 0) == WORD(0));
}

int fmpz_poly_equal_fmpz(const fmpz_poly_t poly, const fmpz_t c);

/*  Addition and subtraction  ************************************************/

void _fmpz_poly_add(fmpz * res, const fmpz * poly1, slong len1,
                                               const fmpz * poly2, slong len2);

void fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t poly1,
                                                      const fmpz_poly_t poly2);

void fmpz_poly_add_series(fmpz_poly_t res, const fmpz_poly_t poly1,
                                             const fmpz_poly_t poly2, slong n);

void _fmpz_poly_sub(fmpz * res, const fmpz * poly1, slong len1,
                                               const fmpz * poly2, slong len2);

void fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1,
                                                      const fmpz_poly_t poly2);

void fmpz_poly_sub_series(fmpz_poly_t res, const fmpz_poly_t poly1,
                                             const fmpz_poly_t poly2, slong n);

void fmpz_poly_neg(fmpz_poly_t res, const fmpz_poly_t poly);

void fmpz_poly_add_si(fmpz_poly_t res, const fmpz_poly_t poly, slong c);
void fmpz_poly_sub_si(fmpz_poly_t res, const fmpz_poly_t poly, slong c);
void fmpz_poly_si_sub(fmpz_poly_t res, slong c, const fmpz_poly_t poly);
void fmpz_poly_add_fmpz(fmpz_poly_t res, const fmpz_poly_t poly, fmpz_t c);
void fmpz_poly_sub_fmpz(fmpz_poly_t res, const fmpz_poly_t poly, fmpz_t c);
void fmpz_poly_fmpz_sub(fmpz_poly_t res, fmpz_t c, const fmpz_poly_t poly);

/*  Scalar absolute value multiplication and division  ***********************/

void fmpz_poly_scalar_abs(fmpz_poly_t res, const fmpz_poly_t poly);

void fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1,
                             const fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_mul_si(fmpz_poly_t poly1,
                             const fmpz_poly_t poly2, slong x);

void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1,
                               const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1,
                                   const fmpz_poly_t poly2, const fmpz_t x);

void fmpz_poly_scalar_addmul_si(fmpz_poly_t poly1,
                                   const fmpz_poly_t poly2, slong x);

void fmpz_poly_scalar_addmul_ui(fmpz_poly_t poly1,
                                   const fmpz_poly_t poly2, ulong x);

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

void fmpz_poly_scalar_mod_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x);
void fmpz_poly_scalar_smod_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x);

slong _fmpz_poly_remove_content_2exp(fmpz * pol, slong len);

void _fmpz_poly_scale_2exp(fmpz * pol, slong len, slong k);

/*  Bit packing  *************************************************************/

void _fmpz_poly_bit_pack(nn_ptr arr, const fmpz * poly,
                                slong len, flint_bitcnt_t bit_size, int negate);

int _fmpz_poly_bit_unpack(fmpz * poly, slong len,
                           nn_srcptr arr, flint_bitcnt_t bit_size, int negate);

void _fmpz_poly_bit_unpack_unsigned(fmpz * poly, slong len,
                                       nn_srcptr arr, flint_bitcnt_t bit_size);

void fmpz_poly_bit_pack(fmpz_t f, const fmpz_poly_t poly,
        flint_bitcnt_t bit_size);

void fmpz_poly_bit_unpack(fmpz_poly_t poly, const fmpz_t f,
        flint_bitcnt_t bit_size);

void fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly, const fmpz_t f,
        flint_bitcnt_t bit_size);


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

void _fmpz_poly_mullow_karatsuba(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, slong len2, slong n);

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

void _fmpz_poly_mulhigh(fmpz * res,
                          const fmpz * poly1, slong len1,
                                   const fmpz * poly2, slong len2, slong start);

/* FFT precached multiplication **********************************************/

void fmpz_poly_mul_SS_precache_init(fmpz_poly_mul_precache_t pre,
                             slong len1, slong bits1, const fmpz_poly_t poly2);

void fmpz_poly_mul_precache_clear(fmpz_poly_mul_precache_t pre);

void _fmpz_poly_mullow_SS_precache(fmpz * output,
   const fmpz * input1, slong len1, fmpz_poly_mul_precache_t pre, slong trunc);

void fmpz_poly_mullow_SS_precache(fmpz_poly_t res,
               const fmpz_poly_t poly1, fmpz_poly_mul_precache_t pre, slong n);

FMPZ_POLY_INLINE void fmpz_poly_mul_SS_precache(fmpz_poly_t res,
                         const fmpz_poly_t poly1, fmpz_poly_mul_precache_t pre)
{
    fmpz_poly_mullow_SS_precache(res, poly1, pre,
		                  FLINT_MAX(poly1->length + pre->len2 - 1, 0));
}

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

void _fmpz_poly_sqrlow_karatsuba(fmpz * res, const fmpz * poly, slong len, slong n);

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

flint_bitcnt_t _fmpz_poly_2norm_normalised_bits(const fmpz * poly, slong len);

ulong fmpz_poly_max_limbs(const fmpz_poly_t poly);
slong fmpz_poly_max_bits(const fmpz_poly_t poly);
void fmpz_poly_height(fmpz_t res, const fmpz_poly_t poly);
slong _fmpz_poly_hamming_weight(const fmpz * a, slong len);

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

void _fmpz_poly_resultant_euclidean(fmpz_t res, const fmpz * poly1, slong len1,
                                               const fmpz * poly2, slong len2);

void fmpz_poly_resultant_euclidean(fmpz_t res, const fmpz_poly_t poly1,
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_resultant_modular(fmpz_t res, const fmpz * poly1, slong len1,
                                               const fmpz * poly2, slong len2);

void fmpz_poly_resultant_modular(fmpz_t res, const fmpz_poly_t poly1,
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_resultant(fmpz_t res, const fmpz * poly1, slong len1,
                                               const fmpz * poly2, slong len2);

void fmpz_poly_resultant(fmpz_t res, const fmpz_poly_t poly1,
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_resultant_modular_div(fmpz_t res,
       const fmpz * poly1, slong len1,
       const fmpz * poly2, slong len2, const fmpz_t divisor, slong nbits);

void fmpz_poly_resultant_modular_div(fmpz_t res,
               const fmpz_poly_t poly1, const fmpz_poly_t poly2,
               const fmpz_t divisor, slong nbits);

void _fmpz_poly_xgcd_modular(fmpz_t r, fmpz * s, fmpz * t,
               const fmpz * poly1, slong len1, const fmpz * poly2, slong len2);

void fmpz_poly_xgcd_modular(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t,
                             const fmpz_poly_t poly1, const fmpz_poly_t poly2);

FMPZ_POLY_INLINE
void _fmpz_poly_xgcd(fmpz_t r, fmpz * s, fmpz * t,
                const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
{
    _fmpz_poly_xgcd_modular(r, s, t, poly1, len1, poly2, len2);
}

FMPZ_POLY_INLINE
void fmpz_poly_xgcd(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t,
                            const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    fmpz_poly_xgcd_modular(r, s, t, poly1, poly2);
}

/*  Discriminant  ********************************************************/

void _fmpz_poly_discriminant(fmpz_t res, const fmpz * poly, slong len);

void fmpz_poly_discriminant(fmpz_t res, const fmpz_poly_t poly);

/*  Gaussian content  ********************************************************/

void _fmpz_poly_content(fmpz_t res, const fmpz * poly, slong len);

void fmpz_poly_content(fmpz_t res, const fmpz_poly_t poly);

void _fmpz_poly_primitive_part(fmpz * res, const fmpz * poly, slong len);

void fmpz_poly_primitive_part(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Square-free  *************************************************************/

int _fmpz_poly_is_squarefree(const fmpz * poly, slong len);

int fmpz_poly_is_squarefree(const fmpz_poly_t poly);

/*  Euclidean division  ******************************************************/

int _fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A,
                            slong lenA, const fmpz * B, slong lenB, int exact);

void fmpz_poly_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R,
                                   const fmpz_poly_t A, const fmpz_poly_t B);

int _fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ,
              fmpz * W, const fmpz * A, const fmpz * B, slong lenB, int exact);

int _fmpz_poly_divrem_divconquer(fmpz * Q, fmpz * R,
            const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact);

void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R,
                                     const fmpz_poly_t A, const fmpz_poly_t B);

int _fmpz_poly_divrem(fmpz * Q, fmpz * R, const fmpz * A, slong lenA,
                                        const fmpz * B, slong lenB, int exact);

void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R,
                                     const fmpz_poly_t A, const fmpz_poly_t B);

int _fmpz_poly_div_basecase(fmpz * Q, fmpz * R, const fmpz * A,
                            slong lenA, const fmpz * B, slong lenB, int exact);

void fmpz_poly_div_basecase(fmpz_poly_t Q,
                                     const fmpz_poly_t A, const fmpz_poly_t B);

int _fmpz_poly_divremlow_divconquer_recursive(fmpz * Q, fmpz * QB,
                        const fmpz * A, const fmpz * B, slong lenB, int exact);

int _fmpz_poly_div_divconquer_recursive(fmpz * Q, fmpz * temp,
                        const fmpz * A, const fmpz * B, slong lenB, int exact);

int _fmpz_poly_div_divconquer(fmpz * Q, const fmpz * A, slong lenA,
                                        const fmpz * B, slong lenB, int exact);

void fmpz_poly_div_divconquer(fmpz_poly_t Q,
                                     const fmpz_poly_t A, const fmpz_poly_t B);

int _fmpz_poly_div(fmpz * Q, const fmpz * A, slong lenA,
                                        const fmpz * B, slong lenB, int exact);

void fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A,
                                                          const fmpz_poly_t B);

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

void _fmpz_poly_divexact(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB);
void fmpz_poly_divexact(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

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

void fmpz_poly_div_root(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_t c);

void _fmpz_poly_div_root(fmpz * Q, const fmpz * A, slong len, const fmpz_t c);

/*  Power series division  ***************************************************/

void _fmpz_poly_inv_series_basecase(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n);

void fmpz_poly_inv_series_basecase(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n);

void _fmpz_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n);

void fmpz_poly_inv_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n);

void _fmpz_poly_inv_series(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n);

void fmpz_poly_inv_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n);

void _fmpz_poly_div_series_basecase(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, slong n);

void _fmpz_poly_div_series_divconquer(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, slong n);

void _fmpz_poly_div_series(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, slong n);

void fmpz_poly_div_series_basecase(fmpz_poly_t Q,
		    const fmpz_poly_t A, const fmpz_poly_t B, slong n);

void fmpz_poly_div_series_divconquer(fmpz_poly_t Q,
		    const fmpz_poly_t A, const fmpz_poly_t B, slong n);

void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A,
                                         const fmpz_poly_t B, slong n);

/*  Divisibility testing  ***************************************************/

int _fmpz_poly_divides(fmpz * q, const fmpz * a,
                                       slong len1, const fmpz * b, slong len2);

int fmpz_poly_divides(fmpz_poly_t q,
		                     const fmpz_poly_t a, const fmpz_poly_t b);

slong fmpz_poly_remove(fmpz_poly_t res, const fmpz_poly_t poly1,
		                                      const fmpz_poly_t poly2);

/*  Pseudo division  *********************************************************/

#ifdef FMPZ_H
void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R,
                   ulong * d, const fmpz * A, slong A_len,
                        const fmpz * B, slong B_len, const fmpz_preinvn_t inv);

void _fmpz_poly_pseudo_divrem_divconquer(fmpz * Q, fmpz * R,
                    ulong * d, const fmpz * A, slong lenA,
                         const fmpz * B, slong lenB, const fmpz_preinvn_t inv);

FMPZ_POLY_INLINE
void _fmpz_poly_pseudo_divrem(fmpz * Q, fmpz * R,
                    ulong * d, const fmpz * A, slong A_len,
                         const fmpz * B, slong B_len, const fmpz_preinvn_t inv)
{
    _fmpz_poly_pseudo_divrem_divconquer(Q, R, d, A, A_len, B, B_len, inv);
}

void _fmpz_poly_pseudo_div(fmpz * Q, ulong * d, const fmpz * A,
             slong lenA, const fmpz * B, slong lenB, const fmpz_preinvn_t inv);

void _fmpz_poly_pseudo_rem(fmpz * R, ulong * d, const fmpz * A,
             slong lenA, const fmpz * B, slong lenB, const fmpz_preinvn_t inv);
#endif

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R,
                          ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

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

FMPZ_POLY_INLINE
void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R,
                           ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
    fmpz_poly_pseudo_divrem_divconquer(Q, R, d, A, B);
}

void fmpz_poly_pseudo_div(fmpz_poly_t Q, ulong * d,
                                     const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_pseudo_rem(fmpz_poly_t R, ulong * d,
                                     const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_divlow_smodp(fmpz * res,
            const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_t p, slong n);

void fmpz_poly_divhigh_smodp(fmpz * res,
            const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_t p, slong n);

/*  Derivative  **************************************************************/

void _fmpz_poly_derivative(fmpz * rpoly, const fmpz * poly, slong len);

void fmpz_poly_derivative(fmpz_poly_t res, const fmpz_poly_t poly);

void _fmpz_poly_nth_derivative(fmpz * rpoly, const fmpz * poly, ulong n, slong len);

void fmpz_poly_nth_derivative(fmpz_poly_t res, const fmpz_poly_t poly, ulong n);

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

void _fmpz_poly_evaluate_horner_fmpq(fmpz_t rnum, fmpz_t rden,
                                    const fmpz * f, slong len,
                                    const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_horner_fmpq(fmpq_t res, const fmpz_poly_t f,
                                                               const fmpq_t a);

void _fmpz_poly_evaluate_divconquer_fmpq(fmpz_t rnum, fmpz_t rden,
                                    const fmpz * f, slong len,
                                    const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_divconquer_fmpq(fmpq_t res,
                                          const fmpz_poly_t f, const fmpq_t a);

void _fmpz_poly_evaluate_fmpq(fmpz_t rnum, fmpz_t rden,
                                    const fmpz * f, slong len,
                                    const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_fmpq(fmpq_t res,
                                          const fmpz_poly_t f, const fmpq_t a);

ulong _fmpz_poly_evaluate_mod(const fmpz * poly, slong len,
                                     ulong a, ulong n, ulong ninv);

ulong fmpz_poly_evaluate_mod(const fmpz_poly_t poly, ulong a,
                                 ulong n);

double _fmpz_poly_evaluate_horner_d(const fmpz * poly, slong n,
                                                                     double d);

double fmpz_poly_evaluate_horner_d(const fmpz_poly_t poly, double d);

double _fmpz_poly_evaluate_horner_d_2exp(slong * exp,
                                         const fmpz * poly, slong n, double d);

double fmpz_poly_evaluate_horner_d_2exp(slong * exp,
                                             const fmpz_poly_t poly, double d);

double _fmpz_poly_evaluate_horner_d_2exp2(slong * exp, const fmpz * poly,
                                      slong n, double d, slong dexp);

double fmpz_poly_evaluate_horner_d_2exp2(slong * exp,
		     const fmpz_poly_t poly, double d, slong dexp);

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

void _fmpz_poly_taylor_shift_multi_mod(fmpz * poly, const fmpz_t c, slong n);

FMPZ_POLY_INLINE
void fmpz_poly_taylor_shift_multi_mod(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c)
{
    if (f != g)
        fmpz_poly_set(g, f);
    _fmpz_poly_taylor_shift_multi_mod(g->coeffs, c, g->length);
}

void _fmpz_poly_taylor_shift(fmpz * poly, const fmpz_t c, slong n);

void fmpz_poly_taylor_shift(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c);

/*  Power series composition and compositional inverse  **********************/

void _fmpz_poly_compose_series_brent_kung(fmpz * res, const fmpz * poly1, slong len1,
                                      const fmpz * poly2, slong len2, slong n);

void fmpz_poly_compose_series_brent_kung(fmpz_poly_t res,
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void _fmpz_poly_compose_series_horner(fmpz * res, const fmpz * poly1, slong len1,
                                      const fmpz * poly2, slong len2, slong n);

void fmpz_poly_compose_series_horner(fmpz_poly_t res,
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void _fmpz_poly_compose_series(fmpz * res, const fmpz * poly1, slong len1,
                                      const fmpz * poly2, slong len2, slong n);

void fmpz_poly_compose_series(fmpz_poly_t res,
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n);

void _fmpz_poly_revert_series(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n);
void fmpz_poly_revert_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n);

/*  Square root  *************************************************************/

int _fmpz_poly_sqrtrem_classical(fmpz * res, fmpz * r,
                                                 const fmpz * poly, slong len);

int fmpz_poly_sqrtrem_classical(fmpz_poly_t b,
                                           fmpz_poly_t r, const fmpz_poly_t a);

int _fmpz_poly_sqrtrem_divconquer(fmpz * res, fmpz * r,
                                    const fmpz * poly, slong len, fmpz * temp);

int fmpz_poly_sqrtrem_divconquer(fmpz_poly_t b,
                                           fmpz_poly_t r, const fmpz_poly_t a);

int _fmpz_poly_sqrt_classical(fmpz * res, const fmpz * poly,
                                                         slong len, int exact);

int fmpz_poly_sqrt_classical(fmpz_poly_t b, const fmpz_poly_t a);

int _fmpz_poly_sqrt_divconquer(fmpz * res, const fmpz * poly,
                                                         slong len, int exact);

int fmpz_poly_sqrt_divconquer(fmpz_poly_t b, const fmpz_poly_t a);

int _fmpz_poly_sqrt_KS(fmpz *rop, const fmpz *op, slong len);

int fmpz_poly_sqrt_KS(fmpz_poly_t b, const fmpz_poly_t a);

int _fmpz_poly_sqrt(fmpz * res, const fmpz * poly, slong len);

int fmpz_poly_sqrt(fmpz_poly_t b, const fmpz_poly_t a);

int _fmpz_poly_sqrt_series(fmpz * res,
                                        const fmpz * poly, slong len, slong n);

int fmpz_poly_sqrt_series(fmpz_poly_t b,
                                                 const fmpz_poly_t a, slong n);

/* Power sums ****************************************************************/

void _fmpz_poly_power_sums_naive(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_poly_power_sums_naive(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

void fmpz_poly_power_sums(fmpz_poly_t res, const fmpz_poly_t poly, slong n);

void _fmpz_poly_power_sums_to_poly(fmpz * res, const fmpz * poly, slong len);

void fmpz_poly_power_sums_to_poly(fmpz_poly_t res, const fmpz_poly_t Q);

/*  Signature  ***************************************************************/

void _fmpz_poly_signature(slong * r1, slong * r2, const fmpz * poly, slong len);

void fmpz_poly_signature(slong * r1, slong * r2, const fmpz_poly_t poly);

/*  Input and output  ********************************************************/

#ifdef FLINT_HAVE_FILE
int _fmpz_poly_fprint(FILE * file, const fmpz * poly, slong len);
int fmpz_poly_fprint(FILE * file, const fmpz_poly_t poly);

int _fmpz_poly_fprint_pretty(FILE * file, const fmpz * poly, slong len, const char * x);
int fmpz_poly_fprint_pretty(FILE * file, const fmpz_poly_t poly, const char * x);

int fmpz_poly_fread(FILE * file, fmpz_poly_t poly);

int fmpz_poly_fread_pretty(FILE *file, fmpz_poly_t poly, char **x);
#endif

int _fmpz_poly_print_pretty(const fmpz * poly, slong len, const char * x);
int fmpz_poly_print_pretty(const fmpz_poly_t poly, const char * x);
int _fmpz_poly_print(const fmpz * poly, slong n);
int fmpz_poly_print(const fmpz_poly_t poly);

int fmpz_poly_read(fmpz_poly_t poly);
int fmpz_poly_read_pretty(fmpz_poly_t poly, char **x);

void fmpz_poly_debug(const fmpz_poly_t poly);

/*  CRT  ********************************************************************/

void fmpz_poly_get_nmod_poly(nmod_poly_t res, const fmpz_poly_t poly);

void fmpz_poly_set_nmod_poly(fmpz_poly_t res, const nmod_poly_t poly);

void fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, const nmod_poly_t poly);

void
_fmpz_poly_CRT_ui_precomp(fmpz * res, const fmpz * poly1, slong len1,
               const fmpz_t m1, nn_srcptr poly2, slong len2, ulong m2,
                ulong m2inv, fmpz_t m1m2, ulong c, int sign);

void _fmpz_poly_CRT_ui(fmpz * res, const fmpz * poly1, slong len1,
               const fmpz_t m1, nn_srcptr poly2, slong len2, ulong m2,
                                                    ulong m2inv, int sign);

void fmpz_poly_CRT_ui(fmpz_poly_t res, const fmpz_poly_t poly1,
                                     const fmpz_t m1, const nmod_poly_t poly2,
                                        int sign);


/* Products *****************************************************************/

void _fmpz_poly_product_roots_fmpz_vec(fmpz * poly,
                                        const fmpz * xs, slong n);

void fmpz_poly_product_roots_fmpz_vec(fmpz_poly_t poly,
                                        const fmpz * xs, slong n);

void _fmpz_poly_product_roots_fmpq_vec(fmpz * poly,
                                        const fmpq * xs, slong n);

void fmpz_poly_product_roots_fmpq_vec(fmpz_poly_t poly,
                                        const fmpq * xs, slong n);


/* Newton basis *************************************************************/

void _fmpz_poly_monomial_to_newton(fmpz * poly, const fmpz * roots, slong n);

void _fmpz_poly_newton_to_monomial(fmpz * poly, const fmpz * roots, slong n);


/* Multipoint evaluation and interpolation *********************************/

void fmpz_poly_evaluate_fmpz_vec(fmpz * res, const fmpz_poly_t f,
                                const fmpz * a, slong n);

void fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n);

int _fmpz_poly_interpolate_newton(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n);
int fmpz_poly_interpolate_newton(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n);
int _fmpz_poly_interpolate_multi_mod(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n);
int fmpz_poly_interpolate_multi_mod(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n);
int _fmpz_poly_interpolate(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n);
int fmpz_poly_interpolate(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n);

void _fmpz_poly_interpolate_exact_newton(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n);
void fmpz_poly_interpolate_exact_newton(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n);
void _fmpz_poly_interpolate_exact(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n);
void fmpz_poly_interpolate_exact(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n);

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

/* Roots */

void _fmpz_poly_bound_roots(fmpz_t bound, const fmpz * poly, slong len);

void fmpz_poly_bound_roots(fmpz_t bound, const fmpz_poly_t poly);

void _fmpz_poly_num_real_roots_sturm(slong * n_neg, slong * n_pos, const fmpz * pol, slong len);

slong fmpz_poly_num_real_roots_sturm(const fmpz_poly_t poly);

slong _fmpz_poly_num_real_roots(const fmpz * pol, slong len);

slong fmpz_poly_num_real_roots(const fmpz_poly_t poly);

/* CLD bounds */

void fmpz_poly_CLD_bound(fmpz_t res, const fmpz_poly_t f, slong n);

/* Special polynomials */

void _fmpz_poly_cyclotomic(fmpz * a, ulong n, nn_ptr factors,
                                        slong num_factors, ulong phi);
void fmpz_poly_cyclotomic(fmpz_poly_t poly, ulong n);

ulong _fmpz_poly_is_cyclotomic(const fmpz * poly, slong len);

ulong fmpz_poly_is_cyclotomic(const fmpz_poly_t poly);

void _fmpz_poly_cos_minpoly(fmpz * f, ulong n);

void fmpz_poly_cos_minpoly(fmpz_poly_t f, ulong n);

void _fmpz_poly_swinnerton_dyer(fmpz * T, ulong n);

void fmpz_poly_swinnerton_dyer(fmpz_poly_t poly, ulong n);

void _fmpz_poly_chebyshev_t(fmpz * coeffs, ulong n);

void fmpz_poly_chebyshev_t(fmpz_poly_t poly, ulong n);

void _fmpz_poly_chebyshev_u(fmpz * coeffs, ulong n);

void fmpz_poly_chebyshev_u(fmpz_poly_t poly, ulong n);

void _fmpz_poly_legendre_pt(fmpz * coeffs, ulong n);

void fmpz_poly_legendre_pt(fmpz_poly_t poly, ulong n);

void _fmpz_poly_hermite_h(fmpz * coeffs, ulong n);

void fmpz_poly_hermite_h(fmpz_poly_t poly, ulong n);

void _fmpz_poly_hermite_he(fmpz * coeffs, ulong n);

void fmpz_poly_hermite_he(fmpz_poly_t poly, ulong n);

void _fmpz_poly_fibonacci(fmpz * coeffs, ulong n);

void fmpz_poly_fibonacci(fmpz_poly_t poly, ulong n);

void _fmpz_poly_eta_qexp(fmpz * f, slong e, slong n);

void fmpz_poly_eta_qexp(fmpz_poly_t f, slong e, slong n);

void _fmpz_poly_theta_qexp(fmpz * f, slong e, slong n);

void fmpz_poly_theta_qexp(fmpz_poly_t f, slong e, slong n);

void fmpz_poly_eulerian_polynomial(fmpz_poly_t poly, ulong n);

/* Declare some old functions dead */

#define fmpz_poly_scalar_mul_mpz _Pragma("GCC error \"'fmpz_poly_scalar_mul_mpz' is deprecated. Use 'fmpz_poly_scalar_mul_fmpz' instead.\"")
#define fmpz_poly_scalar_divexact_mpz _Pragma("GCC error \"'fmpz_poly_scalar_divexact_mpz' is deprecated. Use 'fmpz_poly_scalar_divexact_fmpz' instead.\"")
#define fmpz_poly_scalar_fdiv_mpz _Pragma("GCC error \"'fmpz_poly_scalar_fdiv_mpz' is deprecated. Use 'fmpz_poly_scalar_fdiv_fmpz' instead.\"")
#define fmpz_poly_set_coeff_mpz _Pragma("GCC error \"'fmpz_poly_set_coeff_mpz' is deprecated. Use 'fmpz_poly_set_coeff_fmpz' instead.\"")
#define fmpz_poly_get_coeff_mpz _Pragma("GCC error \"'fmpz_poly_get_coeff_mpz' is deprecated. Use 'fmpz_poly_get_coeff_fmpz' instead.\"")
#define fmpz_poly_set_mpz _Pragma("GCC error \"'fmpz_poly_set_mpz' is deprecated. Use 'fmpz_poly_set_fmpz' instead.\"")
#define fmpz_poly_evaluate_mpq _Pragma("GCC error \"'fmpz_poly_evaluate_mpq' is deprecated. Use 'fmpz_poly_evaluate_fmpq' instead.\"")

#ifdef __cplusplus
}
#endif

#endif
