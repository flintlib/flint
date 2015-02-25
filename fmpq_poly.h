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
    Copyright (C) 2011 Fredrik Johansson
 
******************************************************************************/

#ifndef FMPQ_POLY_H
#define FMPQ_POLY_H

#ifdef FMPQ_POLY_INLINES_C
#define FMPQ_POLY_INLINE FLINT_DLL
#else
#define FMPQ_POLY_INLINE static __inline__
#endif

#include <gmp.h>
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Type definitions  ********************************************************/

typedef struct
{
    fmpz * coeffs;
    fmpz_t den;
    slong alloc;
    slong length;
} fmpq_poly_struct;

typedef fmpq_poly_struct fmpq_poly_t[1];

typedef struct
{
   fmpq_poly_struct * powers;
   slong len;
} fmpq_poly_powers_precomp_struct;

typedef fmpq_poly_powers_precomp_struct fmpq_poly_powers_precomp_t[1];

#define WEAK_CANONICALISE_CUTOFF 25600

/*  Memory management  *******************************************************/

FLINT_DLL void fmpq_poly_init(fmpq_poly_t poly);

FLINT_DLL void fmpq_poly_init2(fmpq_poly_t poly, slong alloc);

FLINT_DLL void fmpq_poly_realloc(fmpq_poly_t poly, slong alloc);

FLINT_DLL void fmpq_poly_fit_length(fmpq_poly_t poly, slong len);

FLINT_DLL void _fmpq_poly_set_length(fmpq_poly_t poly, slong len);

FLINT_DLL void fmpq_poly_clear(fmpq_poly_t poly);

FLINT_DLL void _fmpq_poly_normalise(fmpq_poly_t poly);

/*  Accessing numerator and denominator  *************************************/

#define fmpq_poly_numref(poly)  ((poly)->coeffs)

#define fmpq_poly_denref(poly)  ((poly)->den)

FMPQ_POLY_INLINE void
fmpq_poly_get_numerator(fmpz_poly_t res, const fmpq_poly_t poly)
{
    fmpz_poly_fit_length(res, poly->length);
    _fmpz_vec_set(res->coeffs, poly->coeffs, poly->length);
    _fmpz_poly_set_length(res, poly->length);
}

FMPQ_POLY_INLINE void
fmpq_poly_get_denominator(fmpz_t den, const fmpq_poly_t poly)
{
   fmpz_set(den, fmpq_poly_denref(poly));
}

/*  Canonicalisation  *************************************/

FLINT_DLL void _fmpq_poly_canonicalise(fmpz * rpoly, fmpz_t den, slong len);

FLINT_DLL void fmpq_poly_canonicalise(fmpq_poly_t poly);

FLINT_DLL int _fmpq_poly_is_canonical(const fmpz * poly, const fmpz_t den, slong len);

FLINT_DLL int fmpq_poly_is_canonical(const fmpq_poly_t poly);

/*  Polynomial parameters  ***************************************************/

FMPQ_POLY_INLINE
slong fmpq_poly_degree(const fmpq_poly_t poly)
{
    return poly->length - 1;
}

FMPQ_POLY_INLINE
slong fmpq_poly_length(const fmpq_poly_t poly)
{
    return poly->length;
}

/*  Randomisation  ***********************************************************/

FLINT_DLL void fmpq_poly_randtest(fmpq_poly_t f, flint_rand_t state, 
                                                slong len, mp_bitcnt_t bits_in);

FLINT_DLL void fmpq_poly_randtest_unsigned(fmpq_poly_t f, flint_rand_t state, 
                                                slong len, mp_bitcnt_t bits_in);

FLINT_DLL void fmpq_poly_randtest_not_zero(fmpq_poly_t f, flint_rand_t state,
                                                slong len, mp_bitcnt_t bits_in);

/*  Assignment and basic manipulation  ***************************************/

FLINT_DLL void fmpq_poly_set(fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void fmpq_poly_set_si(fmpq_poly_t poly, slong x);

FLINT_DLL void fmpq_poly_set_ui(fmpq_poly_t poly, ulong x);

FLINT_DLL void fmpq_poly_set_fmpz(fmpq_poly_t poly, const fmpz_t x);

FLINT_DLL void fmpq_poly_set_fmpq(fmpq_poly_t poly, const fmpq_t x);

FLINT_DLL void fmpq_poly_set_mpz(fmpq_poly_t poly, const mpz_t x);

FLINT_DLL void fmpq_poly_set_mpq(fmpq_poly_t poly, const mpq_t x);

FLINT_DLL void fmpq_poly_set_fmpz_poly(fmpq_poly_t rop, const fmpz_poly_t op);

FLINT_DLL void _fmpq_poly_set_array_mpq(fmpz * poly, 
                                         fmpz_t den, const mpq_t * a, slong n);

FLINT_DLL void fmpq_poly_set_array_mpq(fmpq_poly_t poly, 
                                                     const mpq_t * a, slong n);

FLINT_DLL int _fmpq_poly_set_str(fmpz * poly, fmpz_t den, const char * str);

FLINT_DLL int fmpq_poly_set_str(fmpq_poly_t poly, const char * str);

FLINT_DLL char * fmpq_poly_get_str(const fmpq_poly_t poly);

FLINT_DLL char * fmpq_poly_get_str_pretty(const fmpq_poly_t poly, 
                                                             const char * var);

FLINT_DLL char * _fmpq_poly_get_str_pretty(const fmpz *poly, 
                                 const fmpz_t den, slong len, const char *var);
								 
FLINT_DLL void fmpq_poly_zero(fmpq_poly_t poly);

FMPQ_POLY_INLINE void fmpq_poly_one(fmpq_poly_t poly)
{
    fmpq_poly_fit_length(poly, 1);
    _fmpq_poly_set_length(poly, 1);
    fmpz_one(poly->coeffs);
    fmpz_one(poly->den);
}

FLINT_DLL void fmpq_poly_neg(fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void fmpq_poly_inv(fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void fmpq_poly_swap(fmpq_poly_t poly1, fmpq_poly_t poly2);

FMPQ_POLY_INLINE 
void fmpq_poly_truncate(fmpq_poly_t poly, slong n)
{
    if (poly->length > n)
    {
        slong i;
        for (i = n; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = n;
        fmpq_poly_canonicalise(poly);
    }
}

FLINT_DLL void fmpq_poly_set_trunc(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void fmpq_poly_get_slice(fmpq_poly_t rop, 
                         const fmpq_poly_t op, slong i, slong j);

FLINT_DLL void fmpq_poly_reverse(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

/*  Getting and setting coefficients  ****************************************/

FLINT_DLL void fmpq_poly_get_coeff_fmpq(fmpq_t x, const fmpq_poly_t poly, slong n);

FLINT_DLL void fmpq_poly_get_coeff_mpq(mpq_t x, const fmpq_poly_t poly, slong n);

FLINT_DLL void fmpq_poly_set_coeff_si(fmpq_poly_t poly, slong n, slong x);

FLINT_DLL void fmpq_poly_set_coeff_ui(fmpq_poly_t poly, slong n, ulong x);

FLINT_DLL void fmpq_poly_set_coeff_fmpz(fmpq_poly_t poly, slong n, const fmpz_t x);

FLINT_DLL void fmpq_poly_set_coeff_fmpq(fmpq_poly_t poly, slong n, const fmpq_t x);

FLINT_DLL void fmpq_poly_set_coeff_mpz(fmpq_poly_t poly, slong n, const mpz_t x);

FLINT_DLL void fmpq_poly_set_coeff_mpq(fmpq_poly_t poly, slong n, const mpq_t x);

/*  Comparison  **************************************************************/

FLINT_DLL int fmpq_poly_equal(const fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL int fmpq_poly_cmp(const fmpq_poly_t left, const fmpq_poly_t right);

FLINT_DLL int _fmpq_poly_equal_trunc(const fmpz * poly1, const fmpz_t den1, slong len1, 
                           const fmpz * poly2, const fmpz_t den2, slong len2, slong n);

FLINT_DLL int fmpq_poly_equal_trunc(const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n);

FLINT_DLL int fmpq_poly_cmp(const fmpq_poly_t left, const fmpq_poly_t right);

FMPQ_POLY_INLINE
int fmpq_poly_is_zero(const fmpq_poly_t poly)
{
    return poly->length == WORD(0);
}

FMPQ_POLY_INLINE
int fmpq_poly_is_one(const fmpq_poly_t poly)
{
    return (poly->length == WORD(1)) && (fmpz_equal(poly->coeffs, poly->den));
}

FMPQ_POLY_INLINE
int fmpq_poly_is_x(const fmpq_poly_t op)
{
    return (op->length) == 2 && (*(op->coeffs + 1) == WORD(1)) && (*(op->coeffs + 0) == WORD(0)) && (*(op->den) == WORD(1));
}

/*  Inlines, see inlines.c  ************************************************/

FLINT_DLL void fmpq_poly_add_si(fmpq_poly_t res, const fmpq_poly_t poly, slong c);

FLINT_DLL void fmpq_poly_sub_si(fmpq_poly_t res, const fmpq_poly_t poly, slong c);

FLINT_DLL void fmpq_poly_si_sub(fmpq_poly_t res, slong c, const fmpq_poly_t poly);

FLINT_DLL void fmpq_poly_add_fmpz(fmpq_poly_t res, const fmpq_poly_t poly, fmpz_t c);

FLINT_DLL void fmpq_poly_sub_fmpz(fmpq_poly_t res, const fmpq_poly_t poly, fmpz_t c);

FLINT_DLL void fmpq_poly_fmpz_sub(fmpq_poly_t res, fmpz_t c, const fmpq_poly_t poly);

FLINT_DLL void fmpq_poly_add_fmpq(fmpq_poly_t res, const fmpq_poly_t poly, fmpq_t c);

FLINT_DLL void fmpq_poly_sub_fmpq(fmpq_poly_t res, const fmpq_poly_t poly, fmpq_t c);

FLINT_DLL void fmpq_poly_fmpq_sub(fmpq_poly_t res, fmpq_t c, const fmpq_poly_t poly);

/*  Addition and subtraction  ************************************************/

FLINT_DLL void _fmpq_poly_add(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1,
                    const fmpz * poly2, const fmpz_t den2, slong len2);

FLINT_DLL void fmpq_poly_add(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void _fmpq_poly_add_can(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1,
                    const fmpz * poly2, const fmpz_t den2, slong len2, int can);

FLINT_DLL void fmpq_poly_add_can(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2, int can);

FLINT_DLL void _fmpq_poly_sub(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1,
                    const fmpz * poly2, const fmpz_t den2, slong len2);

FLINT_DLL void fmpq_poly_add_series(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n);

FLINT_DLL void _fmpq_poly_add_series_can(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1,
                    const fmpz * poly2, const fmpz_t den2, slong len2, slong n, int can);

FLINT_DLL void fmpq_poly_add_series_can(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n, int can);

FLINT_DLL void _fmpq_poly_sub(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1,
                    const fmpz * poly2, const fmpz_t den2, slong len2);

FLINT_DLL void fmpq_poly_sub(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void _fmpq_poly_sub_can(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1,
                    const fmpz * poly2, const fmpz_t den2, slong len2, int can);

FLINT_DLL void fmpq_poly_sub_can(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2, int can);

FLINT_DLL void _fmpq_poly_sub_series(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1,
                    const fmpz * poly2, const fmpz_t den2, slong len2, slong n);

FLINT_DLL void fmpq_poly_sub_series(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n);

FLINT_DLL void _fmpq_poly_sub_series_can(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1,
                    const fmpz * poly2, const fmpz_t den2, slong len2, slong n, int can);

FLINT_DLL void fmpq_poly_sub_series_can(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n, int can);

/*  Scalar multiplication and division  **************************************/

FLINT_DLL void _fmpq_poly_scalar_mul_si(fmpz * rpoly, fmpz_t rden, 
                       const fmpz * poly, const fmpz_t den, slong len, slong c);

FLINT_DLL void _fmpq_poly_scalar_mul_ui(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, slong len, ulong c);

FLINT_DLL void _fmpq_poly_scalar_mul_fmpz(fmpz * rpoly, fmpz_t rden, 
               const fmpz * poly, const fmpz_t den, slong len, const fmpz_t c);

FLINT_DLL void _fmpq_poly_scalar_mul_fmpq(fmpz * rpoly, fmpz_t rden, const fmpz * poly, 
                  const fmpz_t den, slong len, const fmpz_t r, const fmpz_t s);

FLINT_DLL void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c);

FLINT_DLL void fmpq_poly_scalar_mul_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c);

FLINT_DLL void fmpq_poly_scalar_mul_fmpz(fmpq_poly_t rop, 
                               const fmpq_poly_t op, const fmpz_t c);

FLINT_DLL void fmpq_poly_scalar_mul_fmpq(fmpq_poly_t rop, 
                               const fmpq_poly_t op, const fmpq_t c);

FLINT_DLL void fmpq_poly_scalar_mul_mpz(fmpq_poly_t rop, 
                              const fmpq_poly_t op, const mpz_t c);

FLINT_DLL void fmpq_poly_scalar_mul_mpq(fmpq_poly_t rop, 
                              const fmpq_poly_t op, const mpq_t c);

FLINT_DLL void _fmpq_poly_scalar_div_si(fmpz * rpoly, fmpz_t rden, 
                       const fmpz * poly, const fmpz_t den, slong len, slong c);

FLINT_DLL void _fmpq_poly_scalar_div_ui(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, slong len, ulong c);

FLINT_DLL void _fmpq_poly_scalar_div_fmpz(fmpz * rpoly, fmpz_t rden, 
               const fmpz * poly, const fmpz_t den, slong len, const fmpz_t c);

FLINT_DLL void _fmpq_poly_scalar_div_fmpq(fmpz * rpoly, fmpz_t rden, const fmpz * poly, 
                  const fmpz_t den, slong len, const fmpz_t r, const fmpz_t s);

FLINT_DLL void fmpq_poly_scalar_div_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c);

FLINT_DLL void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c);

FLINT_DLL void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, 
                               const fmpq_poly_t op, const fmpz_t c);

FLINT_DLL void fmpq_poly_scalar_div_fmpq(fmpq_poly_t rop, 
                               const fmpq_poly_t op, const fmpq_t c);

FLINT_DLL void fmpq_poly_scalar_div_mpz(fmpq_poly_t rop, 
                              const fmpq_poly_t op, const mpz_t c);

FLINT_DLL void fmpq_poly_scalar_div_mpq(fmpq_poly_t rop, 
                              const fmpq_poly_t op, const mpq_t c);

/*  Multiplication  **********************************************************/

FLINT_DLL void _fmpq_poly_mul(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1, 
                    const fmpz * poly2, const fmpz_t den2, slong len2);

FLINT_DLL void fmpq_poly_mul(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void _fmpq_poly_mullow(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, slong len1, 
                    const fmpz * poly2, const fmpz_t den2, slong len2, slong n);

FLINT_DLL void fmpq_poly_mullow(fmpq_poly_t res, 
                   const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n);

FMPQ_POLY_INLINE 
void fmpq_poly_addmul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2)
{
    fmpq_poly_t t;

    fmpq_poly_init(t);
    fmpq_poly_mul(t, op1, op2);
    fmpq_poly_add(rop, rop, t);
    fmpq_poly_clear(t);
}

FMPQ_POLY_INLINE 
void fmpq_poly_submul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2)
{
    fmpq_poly_t t;

    fmpq_poly_init(t);
    fmpq_poly_mul(t, op1, op2);
    fmpq_poly_sub(rop, rop, t);
    fmpq_poly_clear(t);
}

/*  Powering  ****************************************************************/

FLINT_DLL void _fmpq_poly_pow(fmpz * rpoly, fmpz_t rden, const fmpz * poly, 
                                  const fmpz_t den, slong len, ulong e);

FLINT_DLL void fmpq_poly_pow(fmpq_poly_t rpoly, const fmpq_poly_t poly, ulong e);

/*  Shifting  ****************************************************************/

FLINT_DLL void fmpq_poly_shift_left(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void fmpq_poly_shift_right(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

/*  Euclidean division  ******************************************************/

FLINT_DLL void _fmpq_poly_divrem(fmpz * Q, fmpz_t q, fmpz * R, fmpz_t r,
                       const fmpz * A, const fmpz_t a, slong lenA,
  const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv);

FLINT_DLL void fmpq_poly_divrem(fmpq_poly_t Q, fmpq_poly_t R,
                      const fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void _fmpq_poly_div(fmpz * Q, fmpz_t q,
                       const fmpz * A, const fmpz_t a, slong lenA,
  const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv);

FLINT_DLL void fmpq_poly_div(fmpq_poly_t Q,
                      const fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void _fmpq_poly_rem(fmpz * R, fmpz_t r,
                       const fmpz * A, const fmpz_t a, slong lenA,
  const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv);

FLINT_DLL void fmpq_poly_rem(fmpq_poly_t R,
                      const fmpq_poly_t poly1, const fmpq_poly_t poly2);

/*  Precomputed inverse  *****************************************************/

FLINT_DLL fmpq_poly_struct * _fmpq_poly_powers_precompute(const fmpz * B, 
                                                 const fmpz_t denB, slong len);

FLINT_DLL void fmpq_poly_powers_precompute(fmpq_poly_powers_precomp_t pinv, 
                                                             fmpq_poly_t poly);

FLINT_DLL void _fmpq_poly_powers_clear(fmpq_poly_struct * powers, slong len);

FLINT_DLL void fmpq_poly_powers_clear(fmpq_poly_powers_precomp_t pinv);

FLINT_DLL void _fmpq_poly_rem_powers_precomp(fmpz * A, fmpz_t denA, slong m, 
                              const fmpz * B, const fmpz_t denB, slong n, 
                                              fmpq_poly_struct * const powers);

FLINT_DLL void fmpq_poly_rem_powers_precomp(fmpq_poly_t R, const fmpq_poly_t A, 
                  const fmpq_poly_t B, const fmpq_poly_powers_precomp_t B_inv);

/*  Power series division  ***************************************************/

FLINT_DLL void _fmpq_poly_inv_series_newton(fmpz * Qinv, fmpz_t Qinvden, 
                           const fmpz * Q, const fmpz_t Qden, slong Qlen, slong n);

FLINT_DLL void fmpq_poly_inv_series_newton(fmpq_poly_t Qinv, const fmpq_poly_t Q, slong n);

FMPQ_POLY_INLINE void 
_fmpq_poly_inv_series(fmpz * Qinv, fmpz_t Qinvden, 
                      const fmpz * Q, const fmpz_t Qden, slong Qlen, slong n)
{
    _fmpq_poly_inv_series_newton(Qinv, Qinvden, Q, Qden, Qlen, n);
}

FMPQ_POLY_INLINE void 
fmpq_poly_inv_series(fmpq_poly_t Qinv, const fmpq_poly_t Q, slong n)
{
    fmpq_poly_inv_series_newton(Qinv, Q, n);
}

FLINT_DLL void _fmpq_poly_div_series(fmpz * Q, fmpz_t denQ, 
        const fmpz * A, const fmpz_t denA, slong lenA,
        const fmpz * B, const fmpz_t denB, slong lenB, slong n);

FLINT_DLL void fmpq_poly_div_series(fmpq_poly_t Q, const fmpq_poly_t A, 
                                         const fmpq_poly_t B, slong n);

/*  Greatest common divisor **************************************************/

FLINT_DLL void _fmpq_poly_gcd(fmpz *G, fmpz_t denG, 
                    const fmpz *A, slong lenA, const fmpz *B, slong lenB);

FLINT_DLL void fmpq_poly_gcd(fmpq_poly_t G, const fmpq_poly_t A, const fmpq_poly_t B);

FLINT_DLL void _fmpq_poly_xgcd(fmpz *G, fmpz_t denG, 
                     fmpz *S, fmpz_t denS, fmpz *T, fmpz_t denT, 
                     const fmpz *A, const fmpz_t denA, slong lenA, 
                     const fmpz *B, const fmpz_t denB, slong lenB);

FLINT_DLL void fmpq_poly_xgcd(fmpq_poly_t G, fmpq_poly_t S, fmpq_poly_t T, 
                    const fmpq_poly_t A, const fmpq_poly_t B);

FLINT_DLL void _fmpq_poly_lcm(fmpz *G, fmpz_t denG, 
                    const fmpz *A, slong lenA, const fmpz *B, slong lenB);

FLINT_DLL void fmpq_poly_lcm(fmpq_poly_t L, const fmpq_poly_t A, const fmpq_poly_t B);

FLINT_DLL void _fmpq_poly_resultant(fmpz_t rnum, fmpz_t rden, 
                          const fmpz *poly1, const fmpz_t den1, slong len1, 
                          const fmpz *poly2, const fmpz_t den2, slong len2);

FLINT_DLL void fmpq_poly_resultant(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g);

/*  Derivative and integral  *************************************************/

FLINT_DLL void _fmpq_poly_derivative(fmpz * rpoly, fmpz_t rden, 
                           const fmpz * poly, const fmpz_t den, slong len);

FLINT_DLL void fmpq_poly_derivative(fmpq_poly_t res, const fmpq_poly_t poly);

FLINT_DLL void _fmpq_poly_integral(fmpz * rpoly, fmpz_t rden, 
                           const fmpz * poly, const fmpz_t den, slong len);

FLINT_DLL void fmpq_poly_integral(fmpq_poly_t res, const fmpq_poly_t poly);

/*  Square roots  ************************************************************/

FLINT_DLL void  _fmpq_poly_invsqrt_series(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, slong len, slong n);

FLINT_DLL void fmpq_poly_invsqrt_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void  _fmpq_poly_sqrt_series(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, slong len, slong n);

FLINT_DLL void fmpq_poly_sqrt_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);


/*  Transcendental functions  ************************************************/

FLINT_DLL void _fmpq_poly_log_series(fmpz * g, fmpz_t gden, 
                       const fmpz * f, const fmpz_t fden, slong flen, slong n);

FLINT_DLL void fmpq_poly_log_series(fmpq_poly_t res, const fmpq_poly_t f, slong n);

FLINT_DLL void _fmpq_poly_exp_series(fmpz * g, fmpz_t gden,
                        const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_exp_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_atan_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_atan_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_atanh_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_atanh_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_asin_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_asin_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_asinh_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_asinh_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_tan_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_tan_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_sin_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_sin_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_cos_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_cos_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_sinh_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_sinh_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_cosh_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_cosh_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_tanh_series(fmpz * g, fmpz_t gden,
                            const fmpz * h, const fmpz_t hden, slong hlen, slong n);

FLINT_DLL void fmpq_poly_tanh_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n);

/*  Evaluation  **************************************************************/

FLINT_DLL void _fmpq_poly_evaluate_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz * poly, 
                              const fmpz_t den, slong len, const fmpz_t a);

FLINT_DLL void fmpq_poly_evaluate_fmpz(fmpq_t res, const fmpq_poly_t poly, 
                             const fmpz_t a);

FLINT_DLL void _fmpq_poly_evaluate_fmpq(fmpz_t rnum, fmpz_t rden, 
                        const fmpz * poly, const fmpz_t den, slong len, 
                        const fmpz_t anum, const fmpz_t aden);

FLINT_DLL void fmpq_poly_evaluate_fmpq(fmpq_t res, const fmpq_poly_t poly, const fmpq_t a);

FLINT_DLL void fmpq_poly_evaluate_mpz(mpq_t res, const fmpq_poly_t poly, const mpz_t a);

FLINT_DLL void fmpq_poly_evaluate_mpq(mpq_t res, const fmpq_poly_t poly, const mpq_t a);

/*  Interpolation ************************************************************/

FLINT_DLL void _fmpq_poly_interpolate_fmpz_vec(fmpz * poly, fmpz_t den,
                                    const fmpz * xs, const fmpz * ys, slong n);

FLINT_DLL void fmpq_poly_interpolate_fmpz_vec(fmpq_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n);

/*  Composition  *************************************************************/

FLINT_DLL void _fmpq_poly_compose(fmpz * res, fmpz_t den, 
                             const fmpz * poly1, const fmpz_t den1, slong len1, 
                             const fmpz * poly2, const fmpz_t den2, slong len2);

FLINT_DLL void fmpq_poly_compose(fmpq_poly_t res, 
                             const fmpq_poly_t poly1, const fmpq_poly_t poly2);

FLINT_DLL void _fmpq_poly_rescale(fmpz * res, fmpz_t denr, const fmpz * poly, 
             const fmpz_t den, slong len, const fmpz_t xnum, const fmpz_t xden);

FLINT_DLL void fmpq_poly_rescale(fmpq_poly_t res, 
                       const fmpq_poly_t poly, const fmpq_t x);

/*  Power series composition  ************************************************/

FLINT_DLL void _fmpq_poly_compose_series_horner(fmpz * res, fmpz_t den, const fmpz * poly1,
        const fmpz_t den1, slong len1, const fmpz * poly2,
        const fmpz_t den2, slong len2, slong n);

FLINT_DLL void fmpq_poly_compose_series_horner(fmpq_poly_t res, 
                    const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n);

FLINT_DLL void _fmpq_poly_compose_series_brent_kung(fmpz * res, fmpz_t den,
        const fmpz * poly1, const fmpz_t den1, slong len1,
        const fmpz * poly2, const fmpz_t den2, slong len2, slong n);

FLINT_DLL void fmpq_poly_compose_series_brent_kung(fmpq_poly_t res, 
                    const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n);

FLINT_DLL void _fmpq_poly_compose_series(fmpz * res, fmpz_t den,
        const fmpz * poly1, const fmpz_t den1, slong len1,
        const fmpz * poly2, const fmpz_t den2, slong len2, slong n);

FLINT_DLL void fmpq_poly_compose_series(fmpq_poly_t res, 
                    const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n);

/*  Power series reversion  ************************************************/

FLINT_DLL void _fmpq_poly_revert_series_lagrange(fmpz * res, fmpz_t den,
        const fmpz * poly1, const fmpz_t den1, slong len1, slong n);

FLINT_DLL void fmpq_poly_revert_series_lagrange(fmpq_poly_t res, 
                    const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_revert_series_lagrange_fast(fmpz * res, fmpz_t den,
        const fmpz * poly1, const fmpz_t den1, slong len1, slong n);

FLINT_DLL void fmpq_poly_revert_series_lagrange_fast(fmpq_poly_t res, 
                    const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_revert_series_newton(fmpz * res, fmpz_t den,
        const fmpz * poly1, const fmpz_t den1, slong len1, slong n);

FLINT_DLL void fmpq_poly_revert_series_newton(fmpq_poly_t res, 
                    const fmpq_poly_t poly, slong n);

FLINT_DLL void _fmpq_poly_revert_series(fmpz * res, fmpz_t den,
        const fmpz * poly1, const fmpz_t den1, slong len1, slong n);

FLINT_DLL void fmpq_poly_revert_series(fmpq_poly_t res, 
                    const fmpq_poly_t poly, slong n);

/*  Gaussian content  ********************************************************/

FLINT_DLL void _fmpq_poly_content(fmpq_t res, 
                        const fmpz * poly, const fmpz_t den, slong len);

FLINT_DLL void fmpq_poly_content(fmpq_t res, const fmpq_poly_t poly);

FLINT_DLL void _fmpq_poly_primitive_part(fmpz * rpoly, fmpz_t rden, 
                               const fmpz * poly, const fmpz_t den, slong len);

FLINT_DLL void fmpq_poly_primitive_part(fmpq_poly_t res, const fmpq_poly_t poly);

FLINT_DLL int _fmpq_poly_is_monic(const fmpz * poly, const fmpz_t den, slong len);

FLINT_DLL int fmpq_poly_is_monic(const fmpq_poly_t poly);

FLINT_DLL void _fmpq_poly_make_monic(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, slong len);

FLINT_DLL void fmpq_poly_make_monic(fmpq_poly_t res, const fmpq_poly_t poly);

/*  Square-free  *************************************************************/

FLINT_DLL int fmpq_poly_is_squarefree(const fmpq_poly_t poly);

/*  Input and output *********************************************************/

FLINT_DLL int fmpq_poly_debug(const fmpq_poly_t poly);

FLINT_DLL int _fmpq_poly_fprint(FILE * file, 
                      const fmpz * poly, const fmpz_t den, slong len);

int fmpq_poly_fprint(FILE * file, const fmpq_poly_t poly);

FLINT_DLL int _fmpq_poly_fprint_pretty(FILE * file, 
                             const fmpz *poly, const fmpz_t den, slong len, 
                             const char * x);

int fmpq_poly_fprint_pretty(FILE * file, 
                            const fmpq_poly_t poly, const char * var);

FMPQ_POLY_INLINE
int _fmpq_poly_print(const fmpz * poly, const fmpz_t den, slong len)
{
    return _fmpq_poly_fprint(stdout, poly, den, len);
}

FMPQ_POLY_INLINE
int fmpq_poly_print(const fmpq_poly_t poly)
{
    return fmpq_poly_fprint(stdout, poly);
}

FMPQ_POLY_INLINE 
int _fmpq_poly_print_pretty(const fmpz *poly, const fmpz_t den, slong len, 
                            const char * x)
{
    return _fmpq_poly_fprint_pretty(stdout, poly, den, len, x);
}

FMPQ_POLY_INLINE 
int fmpq_poly_print_pretty(const fmpq_poly_t poly, const char * var)
{
    return fmpq_poly_fprint_pretty(stdout, poly, var);
}

FLINT_DLL int fmpq_poly_fread(FILE * file, fmpq_poly_t poly);

FMPQ_POLY_INLINE
int fmpq_poly_read(fmpq_poly_t poly)
{
    return fmpq_poly_fread(stdin, poly);
}

#ifdef __cplusplus
}
#endif

#endif

