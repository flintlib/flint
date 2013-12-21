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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#ifndef FMPZ_MOD_POLY_H
#define FMPZ_MOD_POLY_H

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
    fmpz p;
} fmpz_mod_poly_struct;

typedef fmpz_mod_poly_struct fmpz_mod_poly_t[1];

/*  Initialisation and memory management *************************************/

void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_t p);

void fmpz_mod_poly_init2(fmpz_mod_poly_t poly, const fmpz_t p, slong alloc);

void fmpz_mod_poly_clear(fmpz_mod_poly_t poly);

void fmpz_mod_poly_realloc(fmpz_mod_poly_t poly, slong alloc);

void fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len);

/*  Normalisation and truncation *********************************************/

void _fmpz_mod_poly_normalise(fmpz_mod_poly_t poly);

static __inline__ 
void _fmpz_mod_poly_set_length(fmpz_mod_poly_t poly, slong len)
{
    if (poly->length > len)
    {
        slong i;

        for (i = len; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = len;
}

static __inline__ 
void fmpz_mod_poly_truncate(fmpz_mod_poly_t poly, slong len)
{
    if (poly->length > len)
    {
        slong i;

        for (i = len; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = len;
        _fmpz_mod_poly_normalise(poly);
    }  
}

/*  Randomisation ************************************************************/

void
fmpz_mod_poly_randtest(fmpz_mod_poly_t f, flint_rand_t state, slong len);

void
fmpz_mod_poly_randtest_irreducible(fmpz_mod_poly_t f,
                                   flint_rand_t state, slong len);

void
fmpz_mod_poly_randtest_not_zero(fmpz_mod_poly_t f, 
                                flint_rand_t state, slong len);

void
fmpz_mod_poly_randtest_monic(fmpz_mod_poly_t f, flint_rand_t state, slong len);

void
fmpz_mod_poly_randtest_monic_irreducible(fmpz_mod_poly_t f,
                                         flint_rand_t state, slong len);

void
fmpz_mod_poly_randtest_trinomial(fmpz_mod_poly_t f, flint_rand_t state, slong len);

int
fmpz_mod_poly_randtest_trinomial_irreducible(fmpz_mod_poly_t f,
                                             flint_rand_t state, slong len,
                                             slong max_attempts);

void
fmpz_mod_poly_randtest_pentomial(fmpz_mod_poly_t f, flint_rand_t state, slong len);

int
fmpz_mod_poly_randtest_pentomial_irreducible(fmpz_mod_poly_t f,
                                             flint_rand_t state, slong len,
                                             slong max_attempts);

void
fmpz_mod_poly_randtest_sparse_irreducible(fmpz_mod_poly_t poly,
                                          flint_rand_t state, slong len);

/*  Attributes ***************************************************************/

#define fmpz_mod_poly_modulus(poly)  (&((poly)->p))

static __inline__ 
slong fmpz_mod_poly_degree(const fmpz_mod_poly_t poly)
{
    return poly->length - 1;
}

static __inline__ 
slong fmpz_mod_poly_length(const fmpz_mod_poly_t poly)
{
    return poly->length;
}

static __inline__
fmpz * fmpz_mod_poly_lead(const fmpz_mod_poly_t poly)
{
    if (poly->length)
        return poly->coeffs + (poly->length - 1);
    else
        return NULL;
}

/*  Assignment and basic manipulation ****************************************/

void fmpz_mod_poly_set(fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

void fmpz_mod_poly_swap(fmpz_mod_poly_t poly1, fmpz_mod_poly_t poly2);

void _fmpz_mod_poly_reverse(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_mod_poly_reverse(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, slong n);

static __inline__ 
void fmpz_mod_poly_zero(fmpz_mod_poly_t poly)
{
   _fmpz_mod_poly_set_length(poly, 0);
}

void fmpz_mod_poly_zero_coeffs(fmpz_mod_poly_t poly, slong i, slong j);

/*  Conversion ***************************************************************/

static __inline__ 
void fmpz_mod_poly_set_ui(fmpz_mod_poly_t f, ulong x)
{
    if (x == 0)
    {
        fmpz_mod_poly_zero(f);
    }
    else
    {
        fmpz_mod_poly_fit_length(f, 1);
        _fmpz_mod_poly_set_length(f, 1);
        fmpz_set_ui(f->coeffs, x);
        fmpz_mod(f->coeffs, f->coeffs, &(f->p));
        _fmpz_mod_poly_normalise(f);
    }
}

void fmpz_mod_poly_set_fmpz(fmpz_mod_poly_t poly, const fmpz_t c);

void fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f, const fmpz_poly_t g);

void fmpz_mod_poly_get_fmpz_poly(fmpz_poly_t f, const fmpz_mod_poly_t g);

/*  Comparison ***************************************************************/

static __inline__ 
int fmpz_mod_poly_equal(const fmpz_mod_poly_t poly1, 
                        const fmpz_mod_poly_t poly2)
{
    return fmpz_poly_equal((fmpz_poly_struct *) poly1, 
                           (fmpz_poly_struct *) poly2);
}

static __inline__ 
int fmpz_mod_poly_is_zero(const fmpz_mod_poly_t poly)
{
    return (poly->length == 0);
}

/*  Getting and setting coefficients *****************************************/

void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, slong n, const fmpz_t x);

void fmpz_mod_poly_set_coeff_ui(fmpz_mod_poly_t poly, slong n, ulong x);

static __inline__ 
void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly, slong n)
{
    if (n < poly->length)
        fmpz_set(x, poly->coeffs + n);
    else
        fmpz_zero(x);
}

static __inline__ void fmpz_mod_poly_set_coeff_mpz(fmpz_mod_poly_t poly, slong n,
    const mpz_t x)
{
    fmpz_t t;
    fmpz_init_set_readonly(t, x);
    fmpz_mod_poly_set_coeff_fmpz(poly, n, t);
    fmpz_clear_readonly(t);
}

static __inline__ void fmpz_mod_poly_get_coeff_mpz(mpz_t x, const fmpz_mod_poly_t poly, slong n)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mod_poly_get_coeff_fmpz(t, poly, n);
    fmpz_get_mpz(x, t);
    fmpz_clear(t);
}


/*  Shifting *****************************************************************/

void _fmpz_mod_poly_shift_left(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_mod_poly_shift_left(fmpz_mod_poly_t f, 
                              const fmpz_mod_poly_t g, slong n);

void _fmpz_mod_poly_shift_right(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_mod_poly_shift_right(fmpz_mod_poly_t f, 
                               const fmpz_mod_poly_t g, slong n);

/*  Addition and subtraction *************************************************/

void _fmpz_mod_poly_add(fmpz *res, const fmpz *poly1, slong len1, 
                                   const fmpz *poly2, slong len2, const fmpz_t p);

void fmpz_mod_poly_add(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

void _fmpz_mod_poly_sub(fmpz *res, const fmpz *poly1, slong len1, 
                                   const fmpz *poly2, slong len2, const fmpz_t p);

void fmpz_mod_poly_sub(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

void _fmpz_mod_poly_neg(fmpz *res, const fmpz *poly, slong len, const fmpz_t p);

void fmpz_mod_poly_neg(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

/*  Scalar multiplication ****************************************************/

void _fmpz_mod_poly_scalar_mul_fmpz(fmpz *res, const fmpz *poly, slong len, 
                                    const fmpz_t x, const fmpz_t p);

void fmpz_mod_poly_scalar_mul_fmpz(fmpz_mod_poly_t res, 
    const fmpz_mod_poly_t poly, const fmpz_t x);

/*  Multiplication ***********************************************************/

void _fmpz_mod_poly_mul(fmpz *res, const fmpz *poly1, slong len1, 
                                   const fmpz *poly2, slong len2, const fmpz_t p);

void fmpz_mod_poly_mul(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

void _fmpz_mod_poly_mullow(fmpz *res, const fmpz *poly1, slong len1, 
                                      const fmpz *poly2, slong len2, 
                                      const fmpz_t p, slong n);

void fmpz_mod_poly_mullow(fmpz_mod_poly_t res, 
    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n);

void _fmpz_mod_poly_sqr(fmpz *res, const fmpz *poly, slong len, const fmpz_t p);

void fmpz_mod_poly_sqr(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

void _fmpz_mod_poly_mulmod(fmpz * res, const fmpz * poly1, slong len1,
                           const fmpz * poly2, slong len2, const fmpz * f,
                           slong lenf, const fmpz_t p);

void fmpz_mod_poly_mulmod(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                    const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f);

void _fmpz_mod_poly_mulmod_preinv(fmpz * res, const fmpz * poly1, slong len1,
                    const fmpz * poly2, slong len2, const fmpz * f, slong lenf,
                    const fmpz* finv, slong lenfinv, const fmpz_t p);

void
fmpz_mod_poly_mulmod_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                         const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f,
                         const fmpz_mod_poly_t finv);

/*  Powering *****************************************************************/

void _fmpz_mod_poly_pow(fmpz *rop, const fmpz *op, slong len, ulong e, 
                        const fmpz_t p);

void fmpz_mod_poly_pow(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, ulong e);

void _fmpz_mod_poly_pow_trunc(fmpz * res, const fmpz * poly,
                              ulong e, slong trunc, const fmpz_t p);

void fmpz_mod_poly_pow_trunc(fmpz_mod_poly_t res,
                       const fmpz_mod_poly_t poly, ulong e, slong trunc);

void _fmpz_mod_poly_pow_trunc_binexp(fmpz * res, const fmpz * poly,
                                     ulong e, slong trunc, const fmpz_t p);

void fmpz_mod_poly_pow_trunc_binexp(fmpz_mod_poly_t res,
                              const fmpz_mod_poly_t poly, ulong e, slong trunc);

void
_fmpz_mod_poly_powmod_ui_binexp(fmpz * res, const fmpz * poly,
                                ulong e, const fmpz * f, slong lenf, const fmpz_t p);

void
fmpz_mod_poly_powmod_ui_binexp(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, ulong e,
                         const fmpz_mod_poly_t f);


void
_fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz * res, const fmpz * poly,
                              ulong e, const fmpz * f, slong lenf,
                              const fmpz * finv, slong lenfinv, const fmpz_t p);

void
fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, ulong e,
                         const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv);

void
_fmpz_mod_poly_powmod_fmpz_binexp(fmpz * res, const fmpz * poly,
                                  const fmpz_t e, const fmpz * f,
                                  slong lenf, const fmpz_t p);

void
fmpz_mod_poly_powmod_fmpz_binexp(fmpz_mod_poly_t res,
                           const fmpz_mod_poly_t poly, const fmpz_t e,
                           const fmpz_mod_poly_t f);

void
_fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz * res, const fmpz * poly,
                                  const fmpz_t e, const fmpz * f, slong lenf,
                                  const fmpz* finv, slong lenfinv,
                                  const fmpz_t p);

void
fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz_mod_poly_t res,
                           const fmpz_mod_poly_t poly, const fmpz_t e,
                           const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv);

void
_fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz * res, const fmpz_t e, const fmpz * f,
                                    slong lenf, const fmpz* finv, slong lenfinv,
                                    const fmpz_t p);

void
fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz_mod_poly_t res, const fmpz_t e,
                           const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv);

/*  Division *****************************************************************/

void _fmpz_mod_poly_divrem_basecase(fmpz * Q, fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_divrem_basecase(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

void _fmpz_mod_poly_div_basecase(fmpz * Q, fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_div_basecase(fmpz_mod_poly_t Q, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

void
_fmpz_mod_poly_div_newton_n_preinv (fmpz* Q, const fmpz* A, slong lenA,
                                    const fmpz* B, slong lenB, const fmpz* Binv,
                                    slong lenBinv, const fmpz_t mod);

void
fmpz_mod_poly_div_newton_n_preinv(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                           const fmpz_mod_poly_t B, const fmpz_mod_poly_t Binv);

void
_fmpz_mod_poly_divrem_newton_n_preinv (fmpz* Q, fmpz* R, const fmpz* A,
                            slong lenA, const fmpz* B, slong lenB,
                            const fmpz* Binv, slong lenBinv, const fmpz_t mod);

void
fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                               const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                               const fmpz_mod_poly_t Binv);

ulong fmpz_mod_poly_remove(fmpz_mod_poly_t f, const fmpz_mod_poly_t p);

void _fmpz_mod_poly_rem_basecase(fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_rem_basecase(fmpz_mod_poly_t R, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

void _fmpz_mod_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, fmpz * W, 
    const fmpz * A, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

void _fmpz_mod_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_divrem_divconquer(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                                     const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

static __inline__
void _fmpz_mod_poly_divrem(fmpz *Q, fmpz *R, 
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                           const fmpz_t invB, const fmpz_t p)
{
    _fmpz_mod_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, invB, p);
}

static __inline__ 
void fmpz_mod_poly_divrem(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_divrem_divconquer(Q, R, A, B);
}

void _fmpz_mod_poly_divrem_f(fmpz_t f, fmpz *Q, fmpz *R, 
                             const fmpz *A, slong lenA, 
                             const fmpz *B, slong lenB, const fmpz_t p);

void fmpz_mod_poly_divrem_f(fmpz_t f, fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                            const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

static __inline__ 
void _fmpz_mod_poly_rem(fmpz *R, 
                        const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                        const fmpz_t invB, const fmpz_t p)
{
    fmpz *Q = _fmpz_vec_init(lenA - lenB + 1);
    fmpz *T = _fmpz_vec_init(lenA);

    if (lenA < lenB)
    {
       _fmpz_vec_set(R, A, lenA);
       _fmpz_vec_zero(R + lenA, lenB - 1 - lenA);
    } else
    {
       _fmpz_mod_poly_divrem_divconquer(Q, T, A, lenA, B, lenB, invB, p);
       _fmpz_vec_set(R, T, lenB - 1);
    }

    _fmpz_vec_clear(T, lenA);
    _fmpz_vec_clear(Q, lenA - lenB + 1);
}

static __inline__ 
void fmpz_mod_poly_rem(fmpz_mod_poly_t R, 
                       const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_t Q;

    fmpz_mod_poly_init(Q, &(A->p));
    fmpz_mod_poly_divrem(Q, R, A, B);
    fmpz_mod_poly_clear(Q);
}

void _fmpz_mod_poly_div_newton_n_preinv (fmpz *Q, const fmpz* A, slong lenA,
                                         const fmpz* B, slong lenB, const fmpz* Binv,
                                         slong lenBinv, const fmpz_t p);

void fmpz_mod_poly_div_newton_n_preinv (fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                                        const fmpz_mod_poly_t B, const fmpz_mod_poly_t Binv);


void _fmpz_mod_poly_divrem_newton_n_preinv (fmpz* Q, fmpz* R, const fmpz* A,
                                            slong lenA, const fmpz* B, slong lenB,
                                            const fmpz* Binv, slong lenBinv, const fmpz_t p);

void fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                          const fmpz_mod_poly_t Binv);

/*  Power series inversion ***************************************************/

void 
_fmpz_mod_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, slong n, 
                                 const fmpz_t cinv, const fmpz_t p);

void fmpz_mod_poly_inv_series_newton(fmpz_mod_poly_t Qinv, 
    const fmpz_mod_poly_t Q, slong n);

/*  Greatest common divisor **************************************************/

void fmpz_mod_poly_make_monic(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

slong _fmpz_mod_poly_gcd_euclidean(fmpz *G, const fmpz *A, slong lenA, 
                                           const fmpz *B, slong lenB, 
                                           const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_gcd_euclidean(fmpz_mod_poly_t G, 
                                 const fmpz_mod_poly_t A,
                                 const fmpz_mod_poly_t B);

static __inline__ 
slong _fmpz_mod_poly_gcd(fmpz *G, const fmpz *A, slong lenA, 
                                 const fmpz *B, slong lenB, 
                                 const fmpz_t invB, const fmpz_t p)
{
    return _fmpz_mod_poly_gcd_euclidean(G, A, lenA, B, lenB, invB, p);
}

static __inline__ 
void fmpz_mod_poly_gcd(fmpz_mod_poly_t G, 
                       const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_gcd_euclidean(G, A, B);
}

slong _fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz *G, 
                                    const fmpz *A, slong lenA, 
                                    const fmpz *B, slong lenB, const fmpz_t p);

void fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G, 
                                   const fmpz_mod_poly_t A,
                                   const fmpz_mod_poly_t B);

static __inline__ 
slong _fmpz_mod_poly_gcd_f(fmpz_t f, fmpz *G, 
                          const fmpz *A, slong lenA, 
                          const fmpz *B, slong lenB, const fmpz_t p)
{
    return _fmpz_mod_poly_gcd_euclidean_f(f, G, A, lenA, B, lenB, p);
}

static __inline__ 
void fmpz_mod_poly_gcd_f(fmpz_t f, fmpz_mod_poly_t G, 
                         const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_gcd_euclidean_f(f, G, A, B);
}

slong _fmpz_mod_poly_xgcd_euclidean(fmpz *G, fmpz *S, fmpz *T, 
                                   const fmpz *A, slong lenA, 
                                   const fmpz *B, slong lenB, 
                                   const fmpz_t invB, const fmpz_t p);

void fmpz_mod_poly_xgcd_euclidean(fmpz_mod_poly_t G, 
                             fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

static __inline__ slong 
_fmpz_mod_poly_xgcd(fmpz *G, fmpz *S, fmpz *T, 
                    const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                    const fmpz_t invB, const fmpz_t p)
{
    return _fmpz_mod_poly_xgcd_euclidean(G, S, T, A, lenA, B, lenB, invB, p);
}

static __inline__ void 
fmpz_mod_poly_xgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                   const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_xgcd_euclidean(G, S, T, A, B);
}

slong _fmpz_mod_poly_gcdinv(fmpz *G, fmpz *S, 
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                           const fmpz_t p);

void fmpz_mod_poly_gcdinv(fmpz_mod_poly_t G, fmpz_mod_poly_t S, 
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

int _fmpz_mod_poly_invmod(fmpz *A, 
                          const fmpz *B, slong lenB, 
                          const fmpz *P, slong lenP, const fmpz_t p);

int fmpz_mod_poly_invmod(fmpz_mod_poly_t A, 
                         const fmpz_mod_poly_t B, const fmpz_mod_poly_t P);

/*  Derivative  **************************************************************/

void _fmpz_mod_poly_derivative(fmpz *res, const fmpz *poly, slong len, 
                               const fmpz_t p);
 
void fmpz_mod_poly_derivative(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

/*  Evaluation  **************************************************************/

void _fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz *poly, slong len, 
                                  const fmpz_t a, const fmpz_t p);

void fmpz_mod_poly_evaluate_fmpz(fmpz_t res, 
                                 const fmpz_mod_poly_t poly, const fmpz_t a);

fmpz_poly_struct ** _fmpz_mod_poly_tree_alloc(slong len);

void _fmpz_mod_poly_tree_free(fmpz_poly_struct ** tree, slong len);

void _fmpz_mod_poly_tree_build(fmpz_poly_struct ** tree, 
                             const fmpz * roots, slong len, const fmpz_t mod);

void _fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys, const fmpz * coeffs, 
                        slong len, const fmpz * xs, slong n, const fmpz_t mod);

void fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n);

void _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(fmpz * vs, 
              const fmpz * poly, slong plen, fmpz_poly_struct * const * tree, 
                                                 slong len, const fmpz_t mod);

void _fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, 
    const fmpz * poly, slong plen, const fmpz * xs, slong n, const fmpz_t mod);

void fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n);

void _fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys, const fmpz * coeffs, 
                        slong len, const fmpz * xs, slong n, const fmpz_t mod);

void fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n);



/*  Composition  *************************************************************/

void _fmpz_mod_poly_compose_horner(fmpz *res, const fmpz *poly1, slong len1, 
                                              const fmpz *poly2, slong len2, 
                                              const fmpz_t p);

void fmpz_mod_poly_compose_horner(fmpz_mod_poly_t res, 
                                  const fmpz_mod_poly_t poly1, 
                                  const fmpz_mod_poly_t poly2);

void _fmpz_mod_poly_compose_divconquer(fmpz *res, 
                                       const fmpz *poly1, slong len1, 
                                       const fmpz *poly2, slong len2, 
                                       const fmpz_t p);

void fmpz_mod_poly_compose_divconquer(fmpz_mod_poly_t res, 
                                  const fmpz_mod_poly_t poly1, 
                                  const fmpz_mod_poly_t poly2);

static __inline__
void _fmpz_mod_poly_compose(fmpz *res, const fmpz *poly1, slong len1, 
                                       const fmpz *poly2, slong len2, 
                                       const fmpz_t p)
{
    _fmpz_mod_poly_compose_divconquer(res, poly1, len1, poly2, len2, p);
}

static __inline__
void fmpz_mod_poly_compose(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, 
                                                const fmpz_mod_poly_t poly2)
{
    fmpz_mod_poly_compose_divconquer(res, poly1, poly2);
}

/* Modular composition  ******************************************************/

void
_fmpz_mod_poly_compose_mod(fmpz * res, const fmpz * f, slong lenf, const fmpz * g,
                                       const fmpz * h, slong lenh, const fmpz_t p);

void
fmpz_mod_poly_compose_mod(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                  const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3);

void
_fmpz_mod_poly_compose_mod_brent_kung(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, const fmpz * poly3, slong len3, const fmpz_t p);

void
fmpz_mod_poly_compose_mod_brent_kung(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                             const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3);

void
_fmpz_mod_poly_reduce_matrix_mod_poly (fmpz_mat_t A, const fmpz_mat_t B,
                                   const fmpz_mod_poly_t f);

void
_fmpz_mod_poly_precompute_matrix (fmpz_mat_t A, const fmpz * poly1,
                          const fmpz * poly2, slong len2, const fmpz * poly2inv,
                          slong len2inv, const fmpz_t p);

void
fmpz_mod_poly_precompute_matrix(fmpz_mat_t A, const fmpz_mod_poly_t poly1,
                   const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly2inv);

void
_fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz * res,
         const fmpz * poly1, slong len1, const fmpz_mat_t A, const fmpz * poly3,
         slong len3, const fmpz * poly3inv, slong len3inv, const fmpz_t p);

void
fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mat_t A,
                   const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv);

void
_fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz * res, const fmpz * poly1,
                 slong len1, const fmpz * poly2, const fmpz * poly3, slong len3,
                 const fmpz * poly3inv, slong len3inv, const fmpz_t p);

void
fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                   const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv);

void
_fmpz_mod_poly_compose_mod_horner(fmpz * res, const fmpz * f, slong lenf, const fmpz * g,
                                              const fmpz * h, slong lenh, const fmpz_t p);

void
fmpz_mod_poly_compose_mod_horner(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                         const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3);

/*  Radix conversion *********************************************************/

typedef struct {
    fmpz *V;
    fmpz *W;
    fmpz **Rpow;
    fmpz **Rinv;
    slong degR;
    slong k;
    fmpz invL;
} fmpz_mod_poly_radix_struct;

typedef fmpz_mod_poly_radix_struct fmpz_mod_poly_radix_t[1];

void _fmpz_mod_poly_radix_init(fmpz **Rpow, fmpz **Rinv, 
                    const fmpz *R, slong lenR, slong k, 
                    const fmpz_t invL, const fmpz_t p);

void fmpz_mod_poly_radix_init(fmpz_mod_poly_radix_t D, 
                              const fmpz_mod_poly_t R, slong degF);

void fmpz_mod_poly_radix_clear(fmpz_mod_poly_radix_t D);

void _fmpz_mod_poly_radix(fmpz **B, const fmpz *F, fmpz **Rpow, fmpz **Rinv, 
                          slong degR, slong k, slong i, fmpz *W, const fmpz_t p);

void fmpz_mod_poly_radix(fmpz_mod_poly_struct **B, const fmpz_mod_poly_t F, 
                         const fmpz_mod_poly_radix_t D);

/*  Input and output *********************************************************/

int _fmpz_mod_poly_fprint(FILE * file, const fmpz *poly, slong len, 
                          const fmpz_t p);

int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly);

int fmpz_mod_poly_fread(FILE * file, fmpz_mod_poly_t poly);

static __inline__ 
int fmpz_mod_poly_fprint_pretty(FILE * file, 
                                const fmpz_mod_poly_t poly, const char * x)
{
    return _fmpz_poly_fprint_pretty(file, poly->coeffs, poly->length, x);
}

static __inline__ 
int _fmpz_mod_poly_print(const fmpz *poly, slong len, const fmpz_t p)
{
    return _fmpz_mod_poly_fprint(stdout, poly, len, p);
}

static __inline__
int fmpz_mod_poly_print(const fmpz_mod_poly_t poly)
{
    return fmpz_mod_poly_fprint(stdout, poly);
}

static __inline__
int fmpz_mod_poly_print_pretty(const fmpz_mod_poly_t poly, const char * x)
{
    return fmpz_mod_poly_fprint_pretty(stdout, poly, x);
}

#ifdef __cplusplus
}
#endif

#include "fmpz_mod_poly_factor.h"

#endif

