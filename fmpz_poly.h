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
    Copyright (C) 2009 Andy Novocin
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#ifndef FMPZ_POLY_H
#define FMPZ_POLY_H

#include <mpir.h>
#include "fmpz.h"

/*  Type definitions *********************************************************/

typedef struct
{
    fmpz * coeffs;
    long alloc;
    long length;
} fmpz_poly_struct;

typedef fmpz_poly_struct fmpz_poly_t[1];

/*  Memory management ********************************************************/

void fmpz_poly_init(fmpz_poly_t poly);

void fmpz_poly_init2(fmpz_poly_t poly, long alloc);

void fmpz_poly_realloc(fmpz_poly_t poly, long alloc);

void fmpz_poly_fit_length(fmpz_poly_t poly, long len);

void fmpz_poly_clear(fmpz_poly_t poly);

void _fmpz_poly_normalise(fmpz_poly_t poly);

static __inline__
void _fmpz_poly_set_length(fmpz_poly_t poly, long newlen)
{
    if (poly->length > newlen)
    {
        long i;
        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = newlen;
}

/*  Polynomial parameters  ***************************************************/

static __inline__
long fmpz_poly_length(const fmpz_poly_t poly)
{
    return poly->length;
}

static __inline__
long fmpz_poly_degree(const fmpz_poly_t poly)
{
    return poly->length - 1;
}

/*  Assignment and basic manipulation  ***************************************/

void fmpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_set_ui(fmpz_poly_t poly, ulong c);

void fmpz_poly_set_si(fmpz_poly_t poly, long c);

void fmpz_poly_set_fmpz(fmpz_poly_t poly, const fmpz_t c);

void fmpz_poly_set_mpz(fmpz_poly_t poly, const mpz_t c);

int _fmpz_poly_set_str(fmpz * poly, const char * str);

int fmpz_poly_set_str(fmpz_poly_t poly, const char * str);

char * _fmpz_poly_get_str(const fmpz * poly, long len);

char * fmpz_poly_get_str(const fmpz_poly_t poly);

char * _fmpz_poly_get_str_pretty(const fmpz * poly, long len, const char * x);

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

void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2);

void _fmpz_poly_reverse(fmpz * res, const fmpz * poly, long len, long n);

void fmpz_poly_reverse(fmpz_poly_t res, const fmpz_poly_t poly, long n);

static __inline__
void fmpz_poly_truncate(fmpz_poly_t poly, long newlen)
{
    if (poly->length > newlen)
    {
        long i;
        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = newlen;
        _fmpz_poly_normalise(poly);
    }  
}

/*  Randomisation  ***********************************************************/

void fmpz_poly_randinit(fmpz_randstate_t state);

void fmpz_poly_randclear(fmpz_randstate_t state);

void fmpz_poly_randtest(fmpz_poly_t f, fmpz_randstate_t state, 
                                                long len, mp_bitcnt_t bits);

void fmpz_poly_randtest_unsigned(fmpz_poly_t f, fmpz_randstate_t state, 
                                                long len, mp_bitcnt_t bits);

void fmpz_poly_randtest_not_zero(fmpz_poly_t f, fmpz_randstate_t state,
                                                long len, mp_bitcnt_t bits);

/*  Getting and setting coefficients  ****************************************/

long fmpz_poly_get_coeff_si(const fmpz_poly_t poly, long n);

void fmpz_poly_set_coeff_si(fmpz_poly_t poly, long n, long x);

ulong fmpz_poly_get_coeff_ui(const fmpz_poly_t poly, long n);

void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, long n, ulong x);

void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, long n, const fmpz_t x);

void fmpz_poly_get_coeff_fmpz(fmpz_t x, const fmpz_poly_t poly, long n);

/*  Comparison  **************************************************************/

int fmpz_poly_equal(const fmpz_poly_t poly1, const fmpz_poly_t poly2);

/*  Addition and subtraction  ************************************************/

void _fmpz_poly_add(fmpz * res, const fmpz * poly1, long len1, 
                                             const fmpz * poly2, long len2);

void fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                   const fmpz_poly_t poly2);

void _fmpz_poly_sub(fmpz * res, const fmpz * poly1, long len1, 
                                             const fmpz * poly2, long len2);

void fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                   const fmpz_poly_t poly2);

void fmpz_poly_neg(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Scalar multiplication  ***************************************************/

void fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1, fmpz_poly_t poly2, ulong x);

void fmpz_poly_scalar_mul_si(fmpz_poly_t poly1, fmpz_poly_t poly2, long x);

void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1, 
                                   const fmpz_poly_t poly2, const fmpz_t x);

void _fmpz_poly_scalar_addmul_fmpz(fmpz * poly1, 
                             const fmpz * poly2, long len2, const fmpz_t x);

void fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1, 
                                   const fmpz_poly_t poly2, const fmpz_t x);

/*  Bit packing  *************************************************************/

void _fmpz_poly_bit_pack(mp_ptr arr, const fmpz * poly,
                                long len, mp_bitcnt_t bit_size, int negate);

void _fmpz_poly_bit_unpack(fmpz * poly, long len, 
                           mp_srcptr arr, mp_bitcnt_t bit_size, int negate);

void _fmpz_poly_bit_unpack_unsigned(fmpz * poly, long len, 
                                       mp_srcptr arr, mp_bitcnt_t bit_size);

/*  Multiplication  **********************************************************/

void _fmpz_poly_mul_classical(fmpz * res, const fmpz * poly1, long len1, 
                                             const fmpz * poly2, long len2);

void fmpz_poly_mul_classical(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_classical(fmpz * res, const fmpz * poly1,
                      long len1, const fmpz * poly2, long len2, long trunc);

void fmpz_poly_mullow_classical(fmpz_poly_t res, 
              const fmpz_poly_t poly1, const fmpz_poly_t poly2, long trunc);

void _fmpz_poly_mulhigh_classical(fmpz * res, const fmpz * poly1, 
                      long len1, const fmpz * poly2, long len2, long start);

void fmpz_poly_mulhigh_classical(fmpz_poly_t res, 
              const fmpz_poly_t poly1, const fmpz_poly_t poly2, long start);

void _fmpz_poly_mulmid_classical(fmpz * res, const fmpz * poly1, 
                                  long len1, const fmpz * poly2, long len2);

void fmpz_poly_mulmid_classical(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_mul_karatsuba(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mul_karatsuba(fmpz * res, const fmpz * poly1, 
                                  long len1, const fmpz * poly2, long len2);

void _fmpz_poly_mullow_karatsuba_n(fmpz * res, const fmpz * poly1, 
                                              const fmpz * poly2, long len);

void fmpz_poly_mullow_karatsuba_n(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2, long length);

void _fmpz_poly_mulhigh_karatsuba_n(fmpz * res, const fmpz * poly1, 
                                              const fmpz * poly2, long len);

void fmpz_poly_mulhigh_karatsuba_n(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2, long length);

void _fmpz_poly_mul_KS(fmpz * res, const fmpz * poly1, long len1, 
                                             const fmpz * poly2, long len2);

void fmpz_poly_mul_KS(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_KS(fmpz * res, const fmpz * poly1, long len1, 
                                 const fmpz * poly2, long len2, long trunc);

void fmpz_poly_mullow_KS(fmpz_poly_t res, 
              const fmpz_poly_t poly1, const fmpz_poly_t poly2, long trunc);

void _fmpz_poly_mul(fmpz * res, const fmpz * poly1, 
                                  long len1, const fmpz * poly2, long len2);

void fmpz_poly_mul(fmpz_poly_t res, 
                          const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void _fmpz_poly_mullow_n(fmpz * res, const fmpz * poly1, 
                                            const fmpz * poly2, long trunc);

void fmpz_poly_mullow_n(fmpz_poly_t res, 
              const fmpz_poly_t poly1, const fmpz_poly_t poly2, long trunc);

void fmpz_poly_mulhigh_n(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2, long n);

/*  Powering  ****************************************************************/

void _fmpz_poly_pow_multinomial(fmpz * res, const fmpz * poly, long len, ulong e);

void fmpz_poly_pow_multinomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_binomial(fmpz * res, const fmpz * poly, ulong e);

void fmpz_poly_pow_binomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_binexp(fmpz * res, const fmpz * poly, long len, ulong e);

void fmpz_poly_pow_binexp(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_addchains(fmpz * res, const fmpz * poly, long len, const int * a, int n);

void fmpz_poly_pow_addchains(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

void _fmpz_poly_pow_small(fmpz * res, const fmpz * poly, long len, ulong e);

void _fmpz_poly_pow(fmpz * res, const fmpz * poly, long len, ulong e);

void fmpz_poly_pow(fmpz_poly_t res, const fmpz_poly_t poly, ulong e);

/*  Shifting  ****************************************************************/

void _fmpz_poly_shift_left(fmpz * res, const fmpz * poly, long len, long n);

void _fmpz_poly_shift_right(fmpz * res, const fmpz * poly, long len, long n);

void fmpz_poly_shift_left(fmpz_poly_t res, const fmpz_poly_t poly, long n);

void fmpz_poly_shift_right(fmpz_poly_t res, const fmpz_poly_t poly, long n);

/*  Norms  *******************************************************************/

void fmpz_poly_2norm(fmpz_t res, const fmpz_poly_t poly);

/*  Greatest common divisor  *************************************************/

void _fmpz_poly_gcd_subresultant(fmpz * res, const fmpz * poly1, long len1, 
                                              const fmpz * poly2, long len2);

void fmpz_poly_gcd_subresultant(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

void _fmpz_poly_gcd(fmpz * res, const fmpz * poly1, long len1, 
                                              const fmpz * poly2, long len2);

void fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2);

/*  Gaussian content  ********************************************************/

void _fmpz_poly_content(fmpz_t res, const fmpz * poly, long len);

void fmpz_poly_content(fmpz_t res, const fmpz_poly_t poly);

void _fmpz_poly_primitive_part(fmpz * res, const fmpz * poly, long len);

void fmpz_poly_primitive_part(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Euclidean division  ******************************************************/

void _fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A, 
                                      long lenA, const fmpz * B, long lenB);

void fmpz_poly_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                                  const fmpz_poly_t A, const fmpz_poly_t B);

static __inline__
void _fmpz_poly_div_basecase(fmpz * Q, const fmpz * A, long A_len,
                                                const fmpz * B, long B_len)
{
   _fmpz_poly_divrem_basecase(Q, NULL, A, A_len, B, B_len);
}

void _fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, 
                                const fmpz * A, const fmpz * B, long B_len);

void _fmpz_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
                    const fmpz * A, long A_len, const fmpz * B, long B_len);

void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, 
                                  const fmpz_poly_t A, const fmpz_poly_t B);

/*  Pseudo division  *********************************************************/

void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, ulong * d, 
                    const fmpz * A, long A_len, const fmpz * B, long B_len);

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                       ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, ulong * d, 
                    const fmpz * A, long A_len, const fmpz * B, long B_len);

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                       ulong * d, const fmpz_poly_t A, const fmpz_poly_t B);

static __inline__
void _fmpz_poly_pseudo_divrem(fmpz * Q, fmpz * R, ulong * d, 
                    const fmpz * A, long A_len, const fmpz * B, long B_len)
{
    _fmpz_poly_pseudo_divrem_basecase(Q, R, d, A, A_len, B, B_len);
}

static __inline__
void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R, 
                       ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
    fmpz_poly_pseudo_divrem_basecase(Q, R, d, A, B);
}

void _fmpz_poly_pseudo_divrem_cohen(fmpz * Q, fmpz * R, const fmpz * A, 
                                       long lenA, const fmpz * B, long lenB);

void fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R, 
                                               fmpz_poly_t A, fmpz_poly_t B);

void _fmpz_poly_pseudo_rem_cohen(fmpz * R, const fmpz * A, long lenA, 
                                                  const fmpz * B, long lenB);

void fmpz_poly_pseudo_rem_cohen(fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B);

/*  Derivative  **************************************************************/

void _fmpz_poly_derivative(fmpz * rpoly, const fmpz * poly, long len);
 
void fmpz_poly_derivative(fmpz_poly_t res, const fmpz_poly_t poly);

/*  Evaluation  **************************************************************/

void _fmpz_poly_evaluate_horner(fmpz_t res, const fmpz * f, long len, 
                                                               const fmpz_t a);

void fmpz_poly_evaluate_horner(fmpz_t res, const fmpz_poly_t f, 
                                                               const fmpz_t a);

void _fmpz_poly_evaluate(fmpz_t res, const fmpz * f, long len, const fmpz_t a);

void fmpz_poly_evaluate(fmpz_t res, const fmpz_poly_t f, const fmpz_t a);

void _fmpz_poly_evaluate_horner_mpq(fmpz_t rnum, fmpz_t rden, 
                                    const fmpz * f, long len, 
                                    const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_horner_mpq(mpq_t res, const fmpz_poly_t f, 
                                                                const mpq_t a);

void _fmpz_poly_evaluate_mpq(fmpz_t rnum, fmpz_t rden,
                             const fmpz * f, long len, 
                             const fmpz_t anum, const fmpz_t aden);

void fmpz_poly_evaluate_mpq(mpq_t res, const fmpz_poly_t f, const mpq_t a);

/*  Composition  *************************************************************/

void _fmpz_poly_compose_horner(fmpz * res, const fmpz * poly1, long len1, 
                                                const fmpz * poly2, long len2);

void fmpz_poly_compose_horner(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_compose_divconquer(fmpz * res, const fmpz * poly1, long len1, 
                                                const fmpz * poly2, long len2);

void fmpz_poly_compose_divconquer(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

void _fmpz_poly_compose(fmpz * res, const fmpz * poly1, long len1, 
                                                const fmpz * poly2, long len2);

void fmpz_poly_compose(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                      const fmpz_poly_t poly2);

/*  Signature  ***************************************************************/

void _fmpz_poly_signature(ulong * r1, ulong * r2, fmpz * poly, long len);

void fmpz_poly_signature(ulong * r1, ulong * r2, fmpz_poly_t poly);

/*  Printing  ****************************************************************/

void fmpz_poly_print(const fmpz_poly_t poly);

#endif

