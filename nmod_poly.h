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

    Copyright (C) 2007, David Howden.
    Copyright (C) 2010, 2011 William Hart

******************************************************************************/

#ifndef NMOD_POLY_H
#define NMOD_POLY_H

#include <stdio.h>
#include <mpir.h>
#include "nmod_vec.h"
#include "ulong_extras.h"

#define NMOD_DIVREM_DIVCONQUER_CUTOFF  300
#define NMOD_DIV_DIVCONQUER_CUTOFF     300 /* Must be <= NMOD_DIV_DIVCONQUER_CUTOFF */

static __inline__
long NMOD_DIVREM_BC_ITCH(long lenA, long lenB, nmod_t mod)
{
    mp_bitcnt_t bits;

    bits =
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);
    
    if (bits <= FLINT_BITS)
        return lenA;
    else if (bits <= 2 * FLINT_BITS)
        return 2*(lenA + lenB - 1);
    else
        return 3*(lenA + lenB - 1);
}

static __inline__
long NMOD_DIV_BC_ITCH(long lenA, long lenB, nmod_t mod)
{
    mp_bitcnt_t bits;

    bits =
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(lenA - lenB + 1);
    
    if (bits <= FLINT_BITS)
        return lenA - lenB + 1;
    else if (bits <= 2 * FLINT_BITS)
        return 2*lenA;
    else
        return 3*lenA;
}

static __inline__
long NMOD_DIVREM_DC_ITCH(long lenB, nmod_t mod)
{
    long i = 0;
    
    while (lenB > NMOD_DIVREM_DIVCONQUER_CUTOFF + i)
    {
        lenB = (lenB + 1)/2;
        i++;
    }
    if (lenB > NMOD_DIVREM_DIVCONQUER_CUTOFF)
        lenB = NMOD_DIVREM_DIVCONQUER_CUTOFF;

    return NMOD_DIVREM_BC_ITCH(2*lenB - 1, lenB, mod) + 2*lenB - 1;
}

typedef struct
{
    mp_ptr coeffs;
    long alloc;
    long length;
    nmod_t mod;
} nmod_poly_struct;

typedef nmod_poly_struct nmod_poly_t[1];

/* Memory management  ********************************************************/

void nmod_poly_init(nmod_poly_t poly, mp_limb_t n);

void nmod_poly_init_preinv(nmod_poly_t poly, mp_limb_t n, mp_limb_t ninv);

void nmod_poly_init2(nmod_poly_t poly, mp_limb_t n, long alloc);

void nmod_poly_init2_preinv(nmod_poly_t poly, 
                                      mp_limb_t n, mp_limb_t ninv, long alloc);

void nmod_poly_realloc(nmod_poly_t poly, long alloc);

void nmod_poly_clear(nmod_poly_t poly);

void nmod_poly_fit_length(nmod_poly_t poly, long alloc);

static __inline__
void _nmod_poly_normalise(nmod_poly_t poly)
{
    while (poly->length && (poly->coeffs[poly->length - 1] == 0L))
        poly->length--;
}

/* Polynomial parameters  ****************************************************/

static __inline__
long nmod_poly_length(const nmod_poly_t poly)
{
    return poly->length;
}

static __inline__
long nmod_poly_degree(const nmod_poly_t poly)
{
    return poly->length - 1;
}

static __inline__
mp_limb_t nmod_poly_modulus(const nmod_poly_t poly)
{
    return poly->mod.n;
}

static __inline__
mp_bitcnt_t nmod_poly_max_bits(const nmod_poly_t poly)
{
    return _nmod_vec_max_bits(poly->coeffs, poly->length);
}

/* Assignment and basic manipulation  ****************************************/

static __inline__
void nmod_poly_set(nmod_poly_t a, const nmod_poly_t b)
{
    if (a != b)
    {
        nmod_poly_fit_length(a, b->length);
        mpn_copyi(a->coeffs, b->coeffs, b->length);
        a->length = b->length;
   }
}

static __inline__
void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)
{
    long t;
    mp_ptr tp;

    t = poly1->alloc;
    poly1->alloc = poly2->alloc;
    poly2->alloc = t;

    t = poly1->length;
    poly1->length = poly2->length;
    poly2->length = t;

    tp = poly1->coeffs;
    poly1->coeffs = poly2->coeffs;
    poly2->coeffs = tp;
}

static __inline__
void nmod_poly_zero(nmod_poly_t res)
{
    res->length = 0;
}

static __inline__
void nmod_poly_truncate(nmod_poly_t poly, long len)
{
    if (poly->length > len)
    {
        poly->length = len;
        _nmod_poly_normalise(poly);
    }
}

void _nmod_poly_reverse(mp_ptr output, mp_srcptr input, long len, long m);

void nmod_poly_reverse(nmod_poly_t output, const nmod_poly_t input, long m);

/* Randomisation  ************************************************************/

void nmod_poly_randtest(nmod_poly_t poly, flint_rand_t state, long len);

/* Getting and setting coefficients  *****************************************/

static __inline__
ulong nmod_poly_get_coeff_ui(const nmod_poly_t poly, ulong j)
{
    return (j >= poly->length) ? 0 : poly->coeffs[j];
}

void nmod_poly_set_coeff_ui(nmod_poly_t poly, ulong j, ulong c);

/* Input and output  *********************************************************/

char * nmod_poly_get_str(const nmod_poly_t poly);

int nmod_poly_set_str(const char * s, nmod_poly_t poly);

static __inline__
int nmod_poly_print(const nmod_poly_t a)
{
    int r;
    long i;

    r = printf("%ld %lu", a->length, a->mod.n);

    if (a->length == 0)
        return r;
    else
        if (r > 0)
            r = printf(" ");

    for (i = 0; (r > 0) && (i < a->length); i++)
        r = printf(" %lu", a->coeffs[i]);

    return r;
}

int nmod_poly_fread(FILE * f, nmod_poly_t poly);

static __inline__
int nmod_poly_fprint(FILE * f, const nmod_poly_t poly)
{
    char *s;
    int r;

    s = nmod_poly_get_str(poly);
    r = fputs(s, f);
    free(s);

    return (r < 0) ? r : 1;
}

static __inline__
int nmod_poly_read(nmod_poly_t poly)
{
    return nmod_poly_fread(stdin, poly);
}

/* Comparison  ***************************************************************/

static __inline__
int nmod_poly_equal(const nmod_poly_t a, const nmod_poly_t b)
{
    if (a->length != b->length)
        return 0;

    if (a != b)
        if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
            return 0;

   return 1;
}

static __inline__
int nmod_poly_is_zero(const nmod_poly_t poly)
{
    return (poly->length == 0);
}

/* Shifting  *****************************************************************/

void _nmod_poly_shift_left(mp_ptr res, mp_srcptr poly, long len, long k);

void nmod_poly_shift_left(nmod_poly_t res, const nmod_poly_t poly, long k);

void _nmod_poly_shift_right(mp_ptr res, mp_srcptr poly, long len, long k);

void nmod_poly_shift_right(nmod_poly_t res, const nmod_poly_t poly, long k);

/* Addition and subtraction  *************************************************/

void _nmod_poly_add(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1, 
                                                      const nmod_poly_t poly2);

void _nmod_poly_sub(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, 
                                                      const nmod_poly_t poly2);

void nmod_poly_neg(nmod_poly_t res, const nmod_poly_t poly1);

/* Scalar multiplication and division  ***************************************/

void nmod_poly_scalar_mul_nmod(nmod_poly_t res, 
                                         const nmod_poly_t poly1, mp_limb_t c);

void _nmod_poly_make_monic(mp_ptr output, 
                                   mp_srcptr input, long len, nmod_t mod);

void nmod_poly_make_monic(nmod_poly_t output, const nmod_poly_t input);

/* Bit packing and unpacking  ************************************************/

void _nmod_poly_bit_pack(mp_ptr res, mp_srcptr poly, 
                                              long len, mp_bitcnt_t bits);

void _nmod_poly_bit_unpack(mp_ptr res, long len, 
                                  mp_srcptr mpn, mp_bitcnt_t bits, nmod_t mod);

/* Multiplication  ***********************************************************/

void _nmod_poly_mul_classical(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_mul_classical(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow_classical(mp_ptr res, mp_srcptr poly1, long len1, 
                           mp_srcptr poly2, long len2, long trunc, nmod_t mod);

void nmod_poly_mullow_classical(nmod_poly_t res, 
                 const nmod_poly_t poly1, const nmod_poly_t poly2, long trunc);

void _nmod_poly_mulhigh_classical(mp_ptr res, mp_srcptr poly1, long len1, 
                           mp_srcptr poly2, long len2, long start, nmod_t mod);

void nmod_poly_mulhigh_classical(nmod_poly_t res, 
                 const nmod_poly_t poly1, const nmod_poly_t poly2, long start);

void _nmod_poly_mul_KS(mp_ptr out, mp_srcptr in1, long len1, 
                       mp_srcptr in2, long len2, mp_bitcnt_t bits, nmod_t mod);

void nmod_poly_mul_KS(nmod_poly_t res, 
           const nmod_poly_t poly1, const nmod_poly_t poly2, mp_bitcnt_t bits);

void _nmod_poly_mullow_KS(mp_ptr out, mp_srcptr in1, long len1,
               mp_srcptr in2, long len2, mp_bitcnt_t bits, long n, nmod_t mod);

void nmod_poly_mullow_KS(nmod_poly_t res, const nmod_poly_t poly1, 
                            const nmod_poly_t poly2, mp_bitcnt_t bits, long n);

void _nmod_poly_mul(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_mul(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_mullow(mp_ptr res, mp_srcptr poly1, long len1, 
                           mp_srcptr poly2, long len2, long trunc, nmod_t mod);

void nmod_poly_mullow(nmod_poly_t res, const nmod_poly_t poly1, 
                                          const nmod_poly_t poly2, long trunc);

void _nmod_poly_mulhigh(mp_ptr res, mp_srcptr poly1, long len1, 
                               mp_srcptr poly2, long len2, long n, nmod_t mod);

void nmod_poly_mulhigh(nmod_poly_t res, const nmod_poly_t poly1, 
                                              const nmod_poly_t poly2, long n);

/* Powering  *****************************************************************/

void _nmod_poly_pow_binexp(mp_ptr res, 
                              mp_srcptr poly, long len, ulong e, nmod_t mod);

void nmod_poly_pow_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e);

void _nmod_poly_pow(mp_ptr res, mp_srcptr poly, long len, ulong e, nmod_t mod);

void nmod_poly_pow(nmod_poly_t res, const nmod_poly_t poly, ulong e);

void _nmod_poly_pow_trunc_binexp(mp_ptr res, mp_srcptr poly, 
                                              ulong e, long trunc, nmod_t mod);

void nmod_poly_pow_trunc_binexp(nmod_poly_t res, 
                                  const nmod_poly_t poly, ulong e, long trunc);

void _nmod_poly_pow_trunc(mp_ptr res, mp_srcptr poly, 
                                              ulong e, long trunc, nmod_t mod);

void nmod_poly_pow_trunc(nmod_poly_t res, 
                                  const nmod_poly_t poly, ulong e, long trunc);

/* Division  *****************************************************************/

void _nmod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, mp_ptr W,
                 mp_srcptr A, long A_len, mp_srcptr B, long B_len, nmod_t mod);

void nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, 
                                     const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_divrem_divconquer_recursive(mp_ptr Q, mp_ptr BQ, mp_ptr W,  
                    mp_ptr V, mp_srcptr A, mp_srcptr B, long lenB, nmod_t mod);

void _nmod_poly_divrem_divconquer(mp_ptr Q, mp_ptr R, 
                   mp_srcptr A, long lenA, mp_srcptr B, long lenB, nmod_t mod);

void nmod_poly_divrem_divconquer(nmod_poly_t Q, nmod_poly_t R,
                                     const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_divrem(mp_ptr Q, mp_ptr R, mp_srcptr A, long lenA, 
                                           mp_srcptr B, long lenB, nmod_t mod);

void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R,
                                     const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_div_basecase(mp_ptr Q, mp_ptr W, mp_srcptr A, long A_len, 
                                          mp_srcptr B, long B_len, nmod_t mod);

void nmod_poly_div_basecase(nmod_poly_t Q, const nmod_poly_t A,
                                                          const nmod_poly_t B);

void _nmod_poly_div_divconquer_recursive(mp_ptr Q, mp_ptr W, mp_ptr V,
                              mp_srcptr A, mp_srcptr B, long lenB, nmod_t mod);

void _nmod_poly_div_divconquer(mp_ptr Q, mp_srcptr A, long lenA, 
                                           mp_srcptr B, long lenB, nmod_t mod);

void nmod_poly_div_divconquer(nmod_poly_t Q,
                                     const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_div(mp_ptr Q, mp_srcptr A, long lenA, 
                                           mp_srcptr B, long lenB, nmod_t mod);

void nmod_poly_div(nmod_poly_t Q, const nmod_poly_t A, const nmod_poly_t B);

void _nmod_poly_inv_series_basecase(mp_ptr Qinv, 
                                              mp_srcptr Q, long n, nmod_t mod);

void nmod_poly_inv_series_basecase(nmod_poly_t Qinv, 
                                                  const nmod_poly_t Q, long n);

void _nmod_poly_inv_series_newton(mp_ptr Qinv, 
                                              mp_srcptr Q, long n, nmod_t mod);

void nmod_poly_inv_series_newton(nmod_poly_t Qinv, 
                                                  const nmod_poly_t Q, long n);

static __inline__
void _nmod_poly_inv_series(mp_ptr Qinv, mp_srcptr Q, long n, nmod_t mod)
{
    _nmod_poly_inv_series_newton(Qinv, Q, n, mod);
}

static __inline__
void nmod_poly_inv_series(nmod_poly_t Qinv, const nmod_poly_t Q, long n)
{
    nmod_poly_inv_series_newton(Qinv, Q, n);
}

void _nmod_poly_div_series(mp_ptr Q, mp_srcptr A, mp_srcptr B, 
                                                          long n, nmod_t mod);

void nmod_poly_div_series(nmod_poly_t Q, const nmod_poly_t A, 
                                                 const nmod_poly_t B, long n);

void _nmod_poly_div_newton(mp_ptr Q, mp_srcptr A, long Alen, 
                                          mp_srcptr B, long Blen, nmod_t mod);

void nmod_poly_div_newton(nmod_poly_t Q, const nmod_poly_t A,
                                                         const nmod_poly_t B);

void _nmod_poly_divrem_newton(mp_ptr Q, mp_ptr R, 
                  mp_srcptr A, long Alen, mp_srcptr B, long Blen, nmod_t mod);

void nmod_poly_divrem_newton(nmod_poly_t Q, nmod_poly_t R, 
                                    const nmod_poly_t A, const nmod_poly_t B);

/* Derivative  ***************************************************************/

void _nmod_poly_derivative(mp_ptr x_prime, mp_srcptr x, long len, nmod_t mod);

void nmod_poly_derivative(nmod_poly_t x_prime, const nmod_poly_t x);

void _nmod_poly_integral(mp_ptr x_int, mp_srcptr x, long len, nmod_t mod);

void nmod_poly_integral(nmod_poly_t x_int, const nmod_poly_t x);

/* Evaluation  ***************************************************************/

mp_limb_t _nmod_poly_evaluate_nmod(mp_srcptr poly, 
                                           long len, mp_limb_t c, nmod_t mod);

mp_limb_t nmod_poly_evaluate_nmod(const nmod_poly_t poly, mp_limb_t c);

/* Composition  **************************************************************/

void _nmod_poly_compose_horner(mp_ptr res, mp_srcptr poly1, 
                            long len1, mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_compose_horner(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_compose_divconquer(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);
void nmod_poly_compose_divconquer(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

void _nmod_poly_compose(mp_ptr res, mp_srcptr poly1, long len1, 
                                       mp_srcptr poly2, long len2, nmod_t mod);

void nmod_poly_compose(nmod_poly_t res, 
                             const nmod_poly_t poly1, const nmod_poly_t poly2);

long _nmod_poly_gcd_euclidean(mp_ptr G, 
                   mp_srcptr A, long lenA, mp_srcptr B, long lenB, nmod_t mod);

void nmod_poly_gcd_euclidean(nmod_poly_t G, 
                                     const nmod_poly_t A, const nmod_poly_t B);

/* Square roots **************************************************************/

void _nmod_poly_invsqrt_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);

void nmod_poly_invsqrt_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_sqrt_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);

void nmod_poly_sqrt_series(nmod_poly_t g, const nmod_poly_t h, long n);

/* Transcendental functions **************************************************/

void _nmod_poly_atan_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_atan_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_tan_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_tan_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_asin_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_asin_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_sin_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_sin_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_cos_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_cos_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_asinh_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_asinh_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_atanh_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_atanh_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_sinh_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_sinh_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_cosh_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_cosh_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_tanh_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_tanh_series(nmod_poly_t g, const nmod_poly_t h, long n);

void _nmod_poly_log_series_monomial_ui(mp_ptr res, mp_limb_t coeff,
                ulong power, long n, nmod_t mod);
void nmod_poly_log_series_monomial_ui(nmod_poly_t res, mp_limb_t coeff,
                ulong power, long n);

void _nmod_poly_log_series(mp_ptr res, mp_srcptr f, long n, nmod_t mod);
void nmod_poly_log_series(nmod_poly_t res, const nmod_poly_t f, long n);

void _nmod_poly_exp_series_monomial_ui(mp_ptr res, mp_limb_t coeff,
                ulong power, long n, nmod_t mod);
void nmod_poly_exp_series_monomial_ui(nmod_poly_t res, mp_limb_t coeff,
                ulong power, long n);

void
__nmod_poly_exp_series_prealloc(mp_ptr f, mp_ptr g, mp_srcptr h,
    mp_srcptr hprime, mp_ptr T, mp_ptr U, long n, nmod_t mod, int extend);

void
_nmod_poly_exp_series_basecase(mp_ptr f, mp_srcptr h,
                                    long hlen, long n, nmod_t mod);
void nmod_poly_exp_series_basecase(nmod_poly_t f, const nmod_poly_t h, long n);

void _nmod_poly_exp_series(mp_ptr f, mp_srcptr h, long n, nmod_t mod);
void nmod_poly_exp_series(nmod_poly_t f, const nmod_poly_t h, long n);

#endif
