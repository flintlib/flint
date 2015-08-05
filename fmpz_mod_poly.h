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
    Copyright (C) 2014 William Hart

******************************************************************************/

#ifndef FMPZ_MOD_POLY_H
#define FMPZ_MOD_POLY_H

#ifdef FMPZ_MOD_POLY_INLINES_C
#define FMPZ_MOD_POLY_INLINE FLINT_DLL
#else
#define FMPZ_MOD_POLY_INLINE static __inline__
#endif

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

#define FMPZ_MOD_POLY_HGCD_CUTOFF  128      /* HGCD: Basecase -> Recursion      */
#define FMPZ_MOD_POLY_GCD_CUTOFF  256       /* GCD:  Euclidean -> HGCD          */

#define FMPZ_MOD_POLY_INV_NEWTON_CUTOFF  64 /* Inv series newton: Basecase -> Newton */

/*  Type definitions *********************************************************/

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
    fmpz p;
} fmpz_mod_poly_struct;

typedef fmpz_mod_poly_struct fmpz_mod_poly_t[1];

typedef struct
{
   fmpz_t res;
   fmpz_t lc;
   slong len0;
   slong len1;
   slong off;
} fmpz_mod_poly_res_struct;

typedef fmpz_mod_poly_res_struct fmpz_mod_poly_res_t[1];

typedef struct
{
   fmpz_mod_poly_struct * pow;
   slong len;
} fmpz_mod_poly_frobenius_powers_2exp_struct;

typedef fmpz_mod_poly_frobenius_powers_2exp_struct fmpz_mod_poly_frobenius_powers_2exp_t[1];

typedef struct
{
   fmpz_mod_poly_struct * pow;
   slong len;
} fmpz_mod_poly_frobenius_powers_struct;

typedef fmpz_mod_poly_frobenius_powers_struct fmpz_mod_poly_frobenius_powers_t[1];

typedef struct
{
    fmpz_mat_struct A;
    fmpz_mod_poly_struct poly1;
    fmpz_mod_poly_struct poly2;
    fmpz_mod_poly_struct poly2inv;
}
fmpz_mod_poly_matrix_precompute_arg_t;

typedef struct
{
    fmpz_mat_struct A;
    fmpz_mod_poly_struct res;
    fmpz_mod_poly_struct poly1;
    fmpz_mod_poly_struct poly3;
    fmpz_mod_poly_struct poly3inv;
}
fmpz_mod_poly_compose_mod_precomp_preinv_arg_t;


/*  Initialisation and memory management *************************************/

FLINT_DLL void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_init2(fmpz_mod_poly_t poly, const fmpz_t p, slong alloc);

FLINT_DLL void fmpz_mod_poly_clear(fmpz_mod_poly_t poly);

FLINT_DLL void fmpz_mod_poly_realloc(fmpz_mod_poly_t poly, slong alloc);

FLINT_DLL void fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len);

/*  Normalisation and truncation *********************************************/

FLINT_DLL void _fmpz_mod_poly_normalise(fmpz_mod_poly_t poly);

FMPZ_MOD_POLY_INLINE 
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

FMPZ_MOD_POLY_INLINE 
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

FLINT_DLL void fmpz_mod_poly_set_trunc(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, slong n);

/*  Randomisation ************************************************************/

FLINT_DLL void fmpz_mod_poly_randtest(fmpz_mod_poly_t f, flint_rand_t state, slong len);

FLINT_DLL void fmpz_mod_poly_randtest_irreducible(fmpz_mod_poly_t f,
                                   flint_rand_t state, slong len);

FLINT_DLL void fmpz_mod_poly_randtest_not_zero(fmpz_mod_poly_t f, 
                                flint_rand_t state, slong len);

FLINT_DLL void fmpz_mod_poly_randtest_monic(fmpz_mod_poly_t f, flint_rand_t state, slong len);

FLINT_DLL void fmpz_mod_poly_randtest_monic_irreducible(fmpz_mod_poly_t f,
                                         flint_rand_t state, slong len);

FLINT_DLL void fmpz_mod_poly_randtest_trinomial(fmpz_mod_poly_t f, flint_rand_t state, slong len);

FLINT_DLL int fmpz_mod_poly_randtest_trinomial_irreducible(fmpz_mod_poly_t f,
                                             flint_rand_t state, slong len,
                                             slong max_attempts);

FLINT_DLL void fmpz_mod_poly_randtest_pentomial(fmpz_mod_poly_t f, flint_rand_t state, slong len);

FLINT_DLL int fmpz_mod_poly_randtest_pentomial_irreducible(fmpz_mod_poly_t f,
                                             flint_rand_t state, slong len,
                                             slong max_attempts);

FLINT_DLL void fmpz_mod_poly_randtest_sparse_irreducible(fmpz_mod_poly_t poly,
                                          flint_rand_t state, slong len);

/*  Attributes ***************************************************************/

#define fmpz_mod_poly_modulus(poly)  (&((poly)->p))

FMPZ_MOD_POLY_INLINE 
slong fmpz_mod_poly_degree(const fmpz_mod_poly_t poly)
{
    return poly->length - 1;
}

FMPZ_MOD_POLY_INLINE 
slong fmpz_mod_poly_length(const fmpz_mod_poly_t poly)
{
    return poly->length;
}

FMPZ_MOD_POLY_INLINE
fmpz * fmpz_mod_poly_lead(const fmpz_mod_poly_t poly)
{
    if (poly->length)
        return poly->coeffs + (poly->length - 1);
    else
        return NULL;
}

FMPZ_MOD_POLY_INLINE
int fmpz_mod_poly_is_one(const fmpz_mod_poly_t poly)
{
   return poly->length == 1 && fmpz_is_one(poly->coeffs + 0);
}

FMPZ_MOD_POLY_INLINE
int fmpz_mod_poly_is_x(const fmpz_mod_poly_t op)
{
    return (op->length) == 2 && (*(op->coeffs + 1) == WORD(1)) && (*(op->coeffs + 0) == WORD(0));
}

/*  Assignment and basic manipulation ****************************************/

FLINT_DLL void fmpz_mod_poly_set(fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

FLINT_DLL void fmpz_mod_poly_swap(fmpz_mod_poly_t poly1, fmpz_mod_poly_t poly2);

FLINT_DLL void _fmpz_mod_poly_reverse(fmpz * res, const fmpz * poly, slong len, slong n);

FLINT_DLL void fmpz_mod_poly_reverse(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, slong n);

FMPZ_MOD_POLY_INLINE 
void fmpz_mod_poly_zero(fmpz_mod_poly_t poly)
{
   _fmpz_mod_poly_set_length(poly, 0);
}

FLINT_DLL void fmpz_mod_poly_zero_coeffs(fmpz_mod_poly_t poly, slong i, slong j);

/*  Conversion ***************************************************************/

FMPZ_MOD_POLY_INLINE 
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

FLINT_DLL void fmpz_mod_poly_set_fmpz(fmpz_mod_poly_t poly, const fmpz_t c);

FLINT_DLL void fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f, const fmpz_poly_t g);

FLINT_DLL void fmpz_mod_poly_get_fmpz_poly(fmpz_poly_t f, const fmpz_mod_poly_t g);

/*  Comparison ***************************************************************/

FMPZ_MOD_POLY_INLINE 
int fmpz_mod_poly_equal(const fmpz_mod_poly_t poly1, 
                        const fmpz_mod_poly_t poly2)
{
    return fmpz_poly_equal((fmpz_poly_struct *) poly1, 
                           (fmpz_poly_struct *) poly2);
}

FMPZ_MOD_POLY_INLINE 
int fmpz_mod_poly_equal_trunc(const fmpz_mod_poly_t poly1, 
                        const fmpz_mod_poly_t poly2, slong n)
{
    return fmpz_poly_equal_trunc((fmpz_poly_struct *) poly1, 
                           (fmpz_poly_struct *) poly2, n);
}

FMPZ_MOD_POLY_INLINE 
int fmpz_mod_poly_is_zero(const fmpz_mod_poly_t poly)
{
    return (poly->length == 0);
}

/*  Getting and setting coefficients *****************************************/

FLINT_DLL void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, slong n, const fmpz_t x);

FLINT_DLL void fmpz_mod_poly_set_coeff_ui(fmpz_mod_poly_t poly, slong n, ulong x);

FMPZ_MOD_POLY_INLINE 
void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly, slong n)
{
    if (n < poly->length)
        fmpz_set(x, poly->coeffs + n);
    else
        fmpz_zero(x);
}

FMPZ_MOD_POLY_INLINE void fmpz_mod_poly_set_coeff_mpz(fmpz_mod_poly_t poly, slong n,
    const mpz_t x)
{
    fmpz_t t;
    fmpz_init_set_readonly(t, x);
    fmpz_mod_poly_set_coeff_fmpz(poly, n, t);
    fmpz_clear_readonly(t);
}

FMPZ_MOD_POLY_INLINE void fmpz_mod_poly_get_coeff_mpz(mpz_t x, const fmpz_mod_poly_t poly, slong n)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mod_poly_get_coeff_fmpz(t, poly, n);
    fmpz_get_mpz(x, t);
    fmpz_clear(t);
}


/*  Shifting *****************************************************************/

FLINT_DLL void _fmpz_mod_poly_shift_left(fmpz * res, const fmpz * poly, slong len, slong n);

FLINT_DLL void fmpz_mod_poly_shift_left(fmpz_mod_poly_t f, 
                              const fmpz_mod_poly_t g, slong n);

FLINT_DLL void _fmpz_mod_poly_shift_right(fmpz * res, const fmpz * poly, slong len, slong n);

FLINT_DLL void fmpz_mod_poly_shift_right(fmpz_mod_poly_t f, 
                               const fmpz_mod_poly_t g, slong n);

/*  Addition and subtraction *************************************************/

FLINT_DLL void _fmpz_mod_poly_add(fmpz *res, const fmpz *poly1, slong len1, 
                                   const fmpz *poly2, slong len2, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_add(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

FLINT_DLL void _fmpz_mod_poly_sub(fmpz *res, const fmpz *poly1, slong len1, 
                                   const fmpz *poly2, slong len2, const fmpz_t p);
								   
FLINT_DLL void fmpz_mod_poly_add_series(fmpz_mod_poly_t res, 
               const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n);

FLINT_DLL void _fmpz_mod_poly_sub(fmpz *res, const fmpz *poly1, slong len1, 
                                   const fmpz *poly2, slong len2, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_sub(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

FLINT_DLL void _fmpz_mod_poly_neg(fmpz *res, const fmpz *poly, slong len, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_sub_series(fmpz_mod_poly_t res, 
               const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n);

FLINT_DLL void _fmpz_mod_poly_neg(fmpz *res, const fmpz *poly, slong len, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_neg(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

/*  Scalar multiplication ****************************************************/

FLINT_DLL void _fmpz_mod_poly_scalar_mul_fmpz(fmpz *res, const fmpz *poly, slong len, 
                                    const fmpz_t x, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_scalar_mul_fmpz(fmpz_mod_poly_t res, 
    const fmpz_mod_poly_t poly, const fmpz_t x);

/*  Scalar division ****************************************************/

FLINT_DLL void _fmpz_mod_poly_scalar_div_fmpz(fmpz *res, const fmpz *poly, slong len, 
                                    const fmpz_t x, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_scalar_div_fmpz(fmpz_mod_poly_t res, 
                                   const fmpz_mod_poly_t poly, const fmpz_t x);

/*  Multiplication ***********************************************************/

FLINT_DLL void _fmpz_mod_poly_mul(fmpz *res, const fmpz *poly1, slong len1, 
                                   const fmpz *poly2, slong len2, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_mul(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2);

FLINT_DLL void _fmpz_mod_poly_mullow(fmpz *res, const fmpz *poly1, slong len1, 
                                      const fmpz *poly2, slong len2, 
                                      const fmpz_t p, slong n);

FLINT_DLL void fmpz_mod_poly_mullow(fmpz_mod_poly_t res, 
    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n);

FLINT_DLL void _fmpz_mod_poly_sqr(fmpz *res, const fmpz *poly, slong len, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_sqr(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

FLINT_DLL void _fmpz_mod_poly_mulmod(fmpz * res, const fmpz * poly1, slong len1,
                           const fmpz * poly2, slong len2, const fmpz * f,
                           slong lenf, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_mulmod(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                    const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f);

FLINT_DLL void _fmpz_mod_poly_mulmod_preinv(fmpz * res, const fmpz * poly1, slong len1,
                    const fmpz * poly2, slong len2, const fmpz * f, slong lenf,
                    const fmpz* finv, slong lenfinv, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_mulmod_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                         const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f,
                         const fmpz_mod_poly_t finv);

/*  Powering *****************************************************************/

FLINT_DLL void _fmpz_mod_poly_pow(fmpz *rop, const fmpz *op, slong len, ulong e, 
                        const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_pow(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, ulong e);

FLINT_DLL void _fmpz_mod_poly_pow_trunc(fmpz * res, const fmpz * poly,
                              ulong e, slong trunc, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_pow_trunc(fmpz_mod_poly_t res,
                       const fmpz_mod_poly_t poly, ulong e, slong trunc);

FLINT_DLL void _fmpz_mod_poly_pow_trunc_binexp(fmpz * res, const fmpz * poly,
                                     ulong e, slong trunc, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_pow_trunc_binexp(fmpz_mod_poly_t res,
                              const fmpz_mod_poly_t poly, ulong e, slong trunc);

FLINT_DLL void _fmpz_mod_poly_powmod_ui_binexp(fmpz * res, const fmpz * poly,
                                ulong e, const fmpz * f, slong lenf, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_powmod_ui_binexp(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, ulong e,
                         const fmpz_mod_poly_t f);

FLINT_DLL void _fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz * res, const fmpz * poly,
                              ulong e, const fmpz * f, slong lenf,
                              const fmpz * finv, slong lenfinv, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, ulong e,
                         const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv);

FLINT_DLL void _fmpz_mod_poly_powmod_fmpz_binexp(fmpz * res, const fmpz * poly,
                                  const fmpz_t e, const fmpz * f,
                                  slong lenf, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_powmod_fmpz_binexp(fmpz_mod_poly_t res,
                           const fmpz_mod_poly_t poly, const fmpz_t e,
                           const fmpz_mod_poly_t f);

FLINT_DLL void _fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz * res, const fmpz * poly,
                                  const fmpz_t e, const fmpz * f, slong lenf,
                                  const fmpz* finv, slong lenfinv,
                                  const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz_mod_poly_t res,
                          const fmpz_mod_poly_t poly, const fmpz_t e,
                          const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv);

FLINT_DLL void _fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz * res, const fmpz_t e, const fmpz * f,
                                   slong lenf, const fmpz* finv, slong lenfinv,
                                   const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz_mod_poly_t res, const fmpz_t e,
                          const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv);

FLINT_DLL void fmpz_mod_poly_frobenius_powers_2exp_precomp(fmpz_mod_poly_frobenius_powers_2exp_t pow, 
                 const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, ulong m);


FLINT_DLL void fmpz_mod_poly_frobenius_powers_2exp_clear(fmpz_mod_poly_frobenius_powers_2exp_t pow);

FLINT_DLL void fmpz_mod_poly_frobenius_power(fmpz_mod_poly_t res,
                            fmpz_mod_poly_frobenius_powers_2exp_t pow, 
                                             const fmpz_mod_poly_t f, ulong m);

FLINT_DLL void fmpz_mod_poly_frobenius_powers_precomp(fmpz_mod_poly_frobenius_powers_t pow, 
                  const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, ulong m);

FLINT_DLL void fmpz_mod_poly_frobenius_powers_clear(fmpz_mod_poly_frobenius_powers_t pow);

/*  Division *****************************************************************/

FLINT_DLL void _fmpz_mod_poly_divrem_basecase(fmpz * Q, fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_divrem_basecase(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FLINT_DLL void _fmpz_mod_poly_div_basecase(fmpz * Q, fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_div_basecase(fmpz_mod_poly_t Q, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FLINT_DLL void _fmpz_mod_poly_div_newton_n_preinv (fmpz* Q, const fmpz* A, slong lenA,
                                    const fmpz* B, slong lenB, const fmpz* Binv,
                                    slong lenBinv, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_div_newton_n_preinv(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                           const fmpz_mod_poly_t B, const fmpz_mod_poly_t Binv);

FLINT_DLL void _fmpz_mod_poly_divrem_newton_n_preinv (fmpz* Q, fmpz* R, const fmpz* A,
                            slong lenA, const fmpz* B, slong lenB,
                            const fmpz* Binv, slong lenBinv, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                               const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                               const fmpz_mod_poly_t Binv);

FLINT_DLL ulong fmpz_mod_poly_remove(fmpz_mod_poly_t f, const fmpz_mod_poly_t p);

FLINT_DLL void _fmpz_mod_poly_rem_basecase(fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_rem_basecase(fmpz_mod_poly_t R, 
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FLINT_DLL void _fmpz_mod_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, fmpz * W, 
    const fmpz * A, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

FLINT_DLL void _fmpz_mod_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
    const fmpz * A, slong lenA, const fmpz * B, slong lenB, 
    const fmpz_t invB, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_divrem_divconquer(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                                     const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FMPZ_MOD_POLY_INLINE
void _fmpz_mod_poly_divrem(fmpz *Q, fmpz *R, 
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                           const fmpz_t invB, const fmpz_t p)
{
    _fmpz_mod_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, invB, p);
}

FMPZ_MOD_POLY_INLINE 
void fmpz_mod_poly_divrem(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_divrem_divconquer(Q, R, A, B);
}

FLINT_DLL void _fmpz_mod_poly_divrem_f(fmpz_t f, fmpz *Q, fmpz *R, 
                             const fmpz *A, slong lenA, 
                             const fmpz *B, slong lenB, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_divrem_f(fmpz_t f, fmpz_mod_poly_t Q, fmpz_mod_poly_t R, 
                            const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FMPZ_MOD_POLY_INLINE 
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

FMPZ_MOD_POLY_INLINE 
void fmpz_mod_poly_rem(fmpz_mod_poly_t R, 
                       const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_t Q;

    fmpz_mod_poly_init(Q, &(A->p));
    fmpz_mod_poly_divrem(Q, R, A, B);
    fmpz_mod_poly_clear(Q);
}

FMPZ_MOD_POLY_INLINE 
void fmpz_mod_poly_rem_f(fmpz_t f, fmpz_mod_poly_t R, 
                       const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_t Q;

    fmpz_mod_poly_init(Q, &(A->p));
    fmpz_mod_poly_divrem_f(f, Q, R, A, B);
    fmpz_mod_poly_clear(Q);
}

FLINT_DLL void _fmpz_mod_poly_div_newton_n_preinv (fmpz *Q, const fmpz* A, slong lenA,
                                         const fmpz* B, slong lenB, const fmpz* Binv,
                                         slong lenBinv, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_div_newton_n_preinv (fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                                        const fmpz_mod_poly_t B, const fmpz_mod_poly_t Binv);


FLINT_DLL void _fmpz_mod_poly_divrem_newton_n_preinv (fmpz* Q, fmpz* R, const fmpz* A,
                                            slong lenA, const fmpz* B, slong lenB,
                                            const fmpz* Binv, slong lenBinv, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                          const fmpz_mod_poly_t Binv);

/*  Power series inversion ***************************************************/

FLINT_DLL void _fmpz_mod_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, slong n, 
                                 const fmpz_t cinv, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_inv_series_newton(fmpz_mod_poly_t Qinv, 
    const fmpz_mod_poly_t Q, slong n);

FLINT_DLL void fmpz_mod_poly_inv_series_newton_f(fmpz_t f, fmpz_mod_poly_t Qinv, 
    const fmpz_mod_poly_t Q, slong n);

FMPZ_MOD_POLY_INLINE void 
_fmpz_mod_poly_inv_series(fmpz * Qinv, const fmpz * Q, slong n, 
                                 const fmpz_t cinv, const fmpz_t p)
{
   _fmpz_mod_poly_inv_series_newton(Qinv, Q, n, cinv, p);
}

FMPZ_MOD_POLY_INLINE void 
fmpz_mod_poly_inv_series(fmpz_mod_poly_t Qinv, 
    const fmpz_mod_poly_t Q, slong n)
{
   fmpz_mod_poly_inv_series_newton(Qinv, Q, n);
}

FMPZ_MOD_POLY_INLINE void 
fmpz_mod_poly_inv_series_f(fmpz_t f, fmpz_mod_poly_t Qinv, 
    const fmpz_mod_poly_t Q, slong n)
{
   fmpz_mod_poly_inv_series_newton_f(f, Qinv, Q, n);
}

/*  Power series division ***************************************************/

FLINT_DLL void _fmpz_mod_poly_div_series(fmpz * Q, const fmpz * A, slong Alen,
                      const fmpz * B, slong Blen, const fmpz_t p, slong n);

FLINT_DLL void fmpz_mod_poly_div_series(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, 
                                         const fmpz_mod_poly_t B, slong n);

/*  Greatest common divisor **************************************************/

FLINT_DLL void fmpz_mod_poly_make_monic(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

FLINT_DLL void fmpz_mod_poly_make_monic_f(fmpz_t f, fmpz_mod_poly_t res, 
                                                   const fmpz_mod_poly_t poly);

FLINT_DLL slong _fmpz_mod_poly_gcd_euclidean(fmpz *G, const fmpz *A, slong lenA, 
                                           const fmpz *B, slong lenB, 
                                           const fmpz_t invB, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_gcd_euclidean(fmpz_mod_poly_t G, 
                                 const fmpz_mod_poly_t A,
                                 const fmpz_mod_poly_t B);

FLINT_DLL slong _fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz *G, 
                                    const fmpz *A, slong lenA, 
                                    const fmpz *B, slong lenB, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G, 
                                   const fmpz_mod_poly_t A,
                                   const fmpz_mod_poly_t B);

FMPZ_MOD_POLY_INLINE 
slong _fmpz_mod_poly_gcd_f(fmpz_t f, fmpz *G, 
                          const fmpz *A, slong lenA, 
                          const fmpz *B, slong lenB, const fmpz_t p)
{
    return _fmpz_mod_poly_gcd_euclidean_f(f, G, A, lenA, B, lenB, p);
}

FMPZ_MOD_POLY_INLINE 
void fmpz_mod_poly_gcd_f(fmpz_t f, fmpz_mod_poly_t G, 
                         const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_gcd_euclidean_f(f, G, A, B);
}

FLINT_DLL slong _fmpz_mod_poly_hgcd_recursive(fmpz **M, slong *lenM, 
    fmpz *A, slong *lenA, fmpz *B, slong *lenB, 
    const fmpz *a, slong lena, const fmpz *b, slong lenb, 
    fmpz *P, const fmpz_t mod, int flag, fmpz_mod_poly_res_t res);

FLINT_DLL slong _fmpz_mod_poly_hgcd(fmpz **M, slong *lenM, 
                     fmpz *A, slong *lenA, fmpz *B, slong *lenB, 
                     const fmpz *a, slong lena, const fmpz *b, slong lenb, 
                     const fmpz_t mod);

FLINT_DLL slong _fmpz_mod_poly_gcd_hgcd(fmpz *G, const fmpz *A, slong lenA, 
                                  const fmpz *B, slong lenB, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_gcd_hgcd(fmpz_mod_poly_t G, 
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FMPZ_MOD_POLY_INLINE 
slong _fmpz_mod_poly_gcd(fmpz *G, const fmpz *A, slong lenA, 
                                 const fmpz *B, slong lenB, 
                                 const fmpz_t invB, const fmpz_t p)
{
    if (FLINT_MIN(lenA, lenB) < FMPZ_MOD_POLY_GCD_CUTOFF)
       return _fmpz_mod_poly_gcd_euclidean(G, A, lenA, B, lenB, invB, p);
    else
       return _fmpz_mod_poly_gcd_hgcd(G, A, lenA, B, lenB, p);
}

FMPZ_MOD_POLY_INLINE 
void fmpz_mod_poly_gcd(fmpz_mod_poly_t G, 
                       const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    if (FLINT_MIN(A->length, B->length) < FMPZ_MOD_POLY_GCD_CUTOFF)
       fmpz_mod_poly_gcd_euclidean(G, A, B);
    else
       fmpz_mod_poly_gcd_hgcd(G, A, B);
}

FLINT_DLL slong _fmpz_mod_poly_xgcd_euclidean(fmpz *G, fmpz *S, fmpz *T, 
                                   const fmpz *A, slong lenA, 
                                   const fmpz *B, slong lenB, 
                                   const fmpz_t invB, const fmpz_t p);

FLINT_DLL slong _fmpz_mod_poly_xgcd_euclidean_f(fmpz_t f, fmpz *G, fmpz *S, fmpz *T, 
                                   const fmpz *A, slong lenA, 
                                   const fmpz *B, slong lenB, 
                                   const fmpz_t invB, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_xgcd_euclidean(fmpz_mod_poly_t G, 
                             fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FLINT_DLL void fmpz_mod_poly_xgcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G, 
                             fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FLINT_DLL slong _fmpz_mod_poly_xgcd_hgcd(fmpz *G, fmpz *S, fmpz *T, 
                          const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                          const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_xgcd_hgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S, 
          fmpz_mod_poly_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FMPZ_MOD_POLY_INLINE slong 
_fmpz_mod_poly_xgcd(fmpz *G, fmpz *S, fmpz *T, 
                    const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                    const fmpz_t invB, const fmpz_t p)
{
    if (FLINT_MIN(lenA, lenB) < FMPZ_MOD_POLY_GCD_CUTOFF)
       return _fmpz_mod_poly_xgcd_euclidean(G, S, T, A, lenA, B, lenB, invB, p);
    else
       return _fmpz_mod_poly_xgcd_hgcd(G, S, T, A, lenA, B, lenB, p);
}

FMPZ_MOD_POLY_INLINE slong 
_fmpz_mod_poly_xgcd_f(fmpz_t f, fmpz *G, fmpz *S, fmpz *T, 
                    const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                    const fmpz_t invB, const fmpz_t p)
{
    return _fmpz_mod_poly_xgcd_euclidean_f(f, G, S, T, A, lenA, B, lenB, invB, p);
}

FMPZ_MOD_POLY_INLINE void 
fmpz_mod_poly_xgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                   const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    if (FLINT_MIN(A->length, B->length) < FMPZ_MOD_POLY_GCD_CUTOFF)
       fmpz_mod_poly_xgcd_euclidean(G, S, T, A, B);
    else
       fmpz_mod_poly_xgcd_hgcd(G, S, T, A, B);
}

FMPZ_MOD_POLY_INLINE void 
fmpz_mod_poly_xgcd_f(fmpz_t f, fmpz_mod_poly_t G, fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                   const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_poly_xgcd_euclidean_f(f, G, S, T, A, B);
}

FLINT_DLL slong _fmpz_mod_poly_gcdinv(fmpz *G, fmpz *S, 
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                           const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_gcdinv(fmpz_mod_poly_t G, fmpz_mod_poly_t S, 
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FLINT_DLL slong _fmpz_mod_poly_gcdinv_f(fmpz_t f, fmpz *G, fmpz *S, 
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB, 
                           const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_gcdinv_f(fmpz_t f, fmpz_mod_poly_t G, fmpz_mod_poly_t S, 
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B);

FLINT_DLL int _fmpz_mod_poly_invmod(fmpz *A, 
                          const fmpz *B, slong lenB, 
                          const fmpz *P, slong lenP, const fmpz_t p);

FLINT_DLL int _fmpz_mod_poly_invmod_f(fmpz_t f, fmpz *A, 
                          const fmpz *B, slong lenB, 
                          const fmpz *P, slong lenP, const fmpz_t p);

FLINT_DLL int fmpz_mod_poly_invmod(fmpz_mod_poly_t A, 
                         const fmpz_mod_poly_t B, const fmpz_mod_poly_t P);

FLINT_DLL int fmpz_mod_poly_invmod_f(fmpz_t f, fmpz_mod_poly_t A, 
                         const fmpz_mod_poly_t B, const fmpz_mod_poly_t P);

/*  Resultant  ***************************************************************/

FLINT_DLL void _fmpz_mod_poly_resultant_euclidean(fmpz_t res, 
                                    const fmpz *poly1, slong len1, 
                              const fmpz *poly2, slong len2, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_resultant_euclidean(fmpz_t r, const fmpz_mod_poly_t f, 
                                                      const fmpz_mod_poly_t g);

FLINT_DLL void _fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz *A, slong lenA, 
                                  const fmpz *B, slong lenB, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz_mod_poly_t A, 
                                                      const fmpz_mod_poly_t B);

FMPZ_MOD_POLY_INLINE void 
_fmpz_mod_poly_resultant(fmpz_t res, const fmpz *poly1, slong len1, 
                     const fmpz *poly2, slong len2, const fmpz_t mod)
{
    if (len1 < FMPZ_MOD_POLY_GCD_CUTOFF)
        _fmpz_mod_poly_resultant_euclidean(res, poly1, len1, poly2, len2, mod);
    else
        _fmpz_mod_poly_resultant_hgcd(res, poly1, len1, poly2, len2, mod);
}

FMPZ_MOD_POLY_INLINE void 
fmpz_mod_poly_resultant(fmpz_t res, const fmpz_mod_poly_t f, 
                                                       const fmpz_mod_poly_t g)
{
    if (FLINT_MAX(f->length, g->length) < FMPZ_MOD_POLY_GCD_CUTOFF)
       fmpz_mod_poly_resultant_euclidean(res, f, g);
    else
       fmpz_mod_poly_resultant_hgcd(res, f, g);
}

/*  Discriminant  ************************************************************/

FLINT_DLL void _fmpz_mod_poly_discriminant(fmpz_t d, const fmpz *poly, 
                                                  slong len, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_discriminant(fmpz_t d, const fmpz_mod_poly_t f);

/*  Derivative  **************************************************************/

FLINT_DLL void _fmpz_mod_poly_derivative(fmpz *res, const fmpz *poly, slong len, 
                               const fmpz_t p);
 
FLINT_DLL void fmpz_mod_poly_derivative(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly);

/*  Evaluation  **************************************************************/

FLINT_DLL void _fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz *poly, slong len, 
                                  const fmpz_t a, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_evaluate_fmpz(fmpz_t res, 
                                 const fmpz_mod_poly_t poly, const fmpz_t a);

FLINT_DLL fmpz_poly_struct ** _fmpz_mod_poly_tree_alloc(slong len);

FLINT_DLL void _fmpz_mod_poly_tree_free(fmpz_poly_struct ** tree, slong len);

FLINT_DLL void _fmpz_mod_poly_tree_build(fmpz_poly_struct ** tree, 
                             const fmpz * roots, slong len, const fmpz_t mod);

FLINT_DLL void _fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys, const fmpz * coeffs, 
                        slong len, const fmpz * xs, slong n, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n);

FLINT_DLL void _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(fmpz * vs, 
              const fmpz * poly, slong plen, fmpz_poly_struct * const * tree, 
                                                 slong len, const fmpz_t mod);

FLINT_DLL void _fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, 
    const fmpz * poly, slong plen, const fmpz * xs, slong n, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n);

FLINT_DLL void _fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys, const fmpz * coeffs, 
                        slong len, const fmpz * xs, slong n, const fmpz_t mod);

FLINT_DLL void fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n);



/*  Composition  *************************************************************/

FLINT_DLL void _fmpz_mod_poly_compose_horner(fmpz *res, const fmpz *poly1, slong len1, 
                                              const fmpz *poly2, slong len2, 
                                              const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_compose_horner(fmpz_mod_poly_t res, 
                                  const fmpz_mod_poly_t poly1, 
                                  const fmpz_mod_poly_t poly2);

FLINT_DLL void _fmpz_mod_poly_compose_divconquer(fmpz *res, 
                                       const fmpz *poly1, slong len1, 
                                       const fmpz *poly2, slong len2, 
                                       const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_compose_divconquer(fmpz_mod_poly_t res, 
                                  const fmpz_mod_poly_t poly1, 
                                  const fmpz_mod_poly_t poly2);

FMPZ_MOD_POLY_INLINE
void _fmpz_mod_poly_compose(fmpz *res, const fmpz *poly1, slong len1, 
                                       const fmpz *poly2, slong len2, 
                                       const fmpz_t p)
{
    _fmpz_mod_poly_compose_divconquer(res, poly1, len1, poly2, len2, p);
}

FMPZ_MOD_POLY_INLINE
void fmpz_mod_poly_compose(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, 
                                                const fmpz_mod_poly_t poly2)
{
    fmpz_mod_poly_compose_divconquer(res, poly1, poly2);
}

/* Modular composition  ******************************************************/

FLINT_DLL void _fmpz_mod_poly_compose_mod(fmpz * res, const fmpz * f, slong lenf, const fmpz * g,
                                       const fmpz * h, slong lenh, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_compose_mod(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                  const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3);

FLINT_DLL void _fmpz_mod_poly_compose_mod_brent_kung(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, const fmpz * poly3, slong len3, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_compose_mod_brent_kung(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                             const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3);

FLINT_DLL void _fmpz_mod_poly_reduce_matrix_mod_poly (fmpz_mat_t A, const fmpz_mat_t B,
                                   const fmpz_mod_poly_t f);

FLINT_DLL void _fmpz_mod_poly_precompute_matrix (fmpz_mat_t A, const fmpz * poly1,
                          const fmpz * poly2, slong len2, const fmpz * poly2inv,
                          slong len2inv, const fmpz_t p);

FLINT_DLL void * _fmpz_mod_poly_precompute_matrix_worker(void * arg_ptr);

FLINT_DLL void fmpz_mod_poly_precompute_matrix(fmpz_mat_t A, const fmpz_mod_poly_t poly1,
                   const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly2inv);

FLINT_DLL void _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz * res,
         const fmpz * poly1, slong len1, const fmpz_mat_t A, const fmpz * poly3,
         slong len3, const fmpz * poly3inv, slong len3inv, const fmpz_t p);

FLINT_DLL void * _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker(void * arg_ptr);

FLINT_DLL void fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mat_t A,
                   const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv);

FLINT_DLL void _fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz * res, const fmpz * poly1,
                 slong len1, const fmpz * poly2, const fmpz * poly3, slong len3,
                 const fmpz * poly3inv, slong len3inv, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                   const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv);

FLINT_DLL void _fmpz_mod_poly_compose_mod_horner(fmpz * res, const fmpz * f, slong lenf, const fmpz * g,
                                              const fmpz * h, slong lenh, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_compose_mod_horner(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                         const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3);

FLINT_DLL void _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(fmpz_mod_poly_struct * res,
                 const fmpz_mod_poly_struct * polys, slong len1, slong l,
                 const fmpz * poly, slong len, const fmpz * polyinv,
                 slong leninv, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(fmpz_mod_poly_struct * res,
                    const fmpz_mod_poly_struct * polys, slong len1, slong n,
                    const fmpz_mod_poly_t poly, const fmpz_mod_poly_t polyinv);

FLINT_DLL void _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(fmpz_mod_poly_struct * res,
                                                 const fmpz_mod_poly_struct *
                                                 polys, slong lenpolys,
                                                 slong l, const fmpz * poly,
                                                 slong len,
                                                 const fmpz * polyinv,
                                                 slong leninv, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(fmpz_mod_poly_struct * res,
                                                const fmpz_mod_poly_struct *
                                                polys, slong len1, slong n,
                                                const fmpz_mod_poly_t poly,
                                                const fmpz_mod_poly_t polyinv);

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

FLINT_DLL void _fmpz_mod_poly_radix_init(fmpz **Rpow, fmpz **Rinv, 
                    const fmpz *R, slong lenR, slong k, 
                    const fmpz_t invL, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_radix_init(fmpz_mod_poly_radix_t D, 
                              const fmpz_mod_poly_t R, slong degF);

FLINT_DLL void fmpz_mod_poly_radix_clear(fmpz_mod_poly_radix_t D);

FLINT_DLL void _fmpz_mod_poly_radix(fmpz **B, const fmpz *F, fmpz **Rpow, fmpz **Rinv, 
                          slong degR, slong k, slong i, fmpz *W, const fmpz_t p);

FLINT_DLL void fmpz_mod_poly_radix(fmpz_mod_poly_struct **B, const fmpz_mod_poly_t F, 
                         const fmpz_mod_poly_radix_t D);

/*  Input and output *********************************************************/

FLINT_DLL int _fmpz_mod_poly_fprint(FILE * file, const fmpz *poly, slong len, 
                          const fmpz_t p);

FLINT_DLL int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly);

FLINT_DLL int fmpz_mod_poly_fread(FILE * file, fmpz_mod_poly_t poly);

FMPZ_MOD_POLY_INLINE 
int fmpz_mod_poly_fprint_pretty(FILE * file, 
                                const fmpz_mod_poly_t poly, const char * x)
{
    return _fmpz_poly_fprint_pretty(file, poly->coeffs, poly->length, x);
}

FMPZ_MOD_POLY_INLINE 
int _fmpz_mod_poly_print(const fmpz *poly, slong len, const fmpz_t p)
{
    return _fmpz_mod_poly_fprint(stdout, poly, len, p);
}

FMPZ_MOD_POLY_INLINE
int fmpz_mod_poly_print(const fmpz_mod_poly_t poly)
{
    return fmpz_mod_poly_fprint(stdout, poly);
}

FMPZ_MOD_POLY_INLINE
int fmpz_mod_poly_print_pretty(const fmpz_mod_poly_t poly, const char * x)
{
    return fmpz_mod_poly_fprint_pretty(stdout, poly, x);
}

#ifdef __cplusplus
}
#endif

#include "fmpz_mod_poly_factor.h"

#endif

