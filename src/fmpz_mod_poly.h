/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_POLY_H
#define FMPZ_MOD_POLY_H

#ifdef FMPZ_MOD_POLY_INLINES_C
#define FMPZ_MOD_POLY_INLINE
#else
#define FMPZ_MOD_POLY_INLINE static inline
#endif

#include "thread_pool.h"
#include "nmod_types.h"
#include "fmpz_mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FMPZ_MOD_POLY_HGCD_CUTOFF  128      /* HGCD: Basecase -> Recursion      */
#define FMPZ_MOD_POLY_GCD_CUTOFF  256       /* GCD:  Euclidean -> HGCD          */

#define FMPZ_MOD_POLY_INV_NEWTON_CUTOFF  64 /* Inv series newton: Basecase -> Newton */
#define FMPZ_MOD_POLY_DIV_DIVCONQUER_CUTOFF

#define FMPZ_MOD_POLY_EVALUATE_FMPZ_VEC  32 /* Evaluate fmpz_vec: Iter -> Fast  */

/*  Type definitions *********************************************************/

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
    fmpz_mat_struct * A;
    fmpz_mod_poly_struct * poly1;
    fmpz_mod_poly_struct * poly2;
    fmpz_mod_poly_struct * poly2inv;
    const fmpz_mod_ctx_struct * ctx;
}
fmpz_mod_poly_matrix_precompute_arg_t;

typedef struct
{
    fmpz_mat_struct * A;
    fmpz_mod_poly_struct * res;
    fmpz_mod_poly_struct * poly1;
    fmpz_mod_poly_struct * poly3;
    fmpz_mod_poly_struct * poly3inv;
    const fmpz_mod_ctx_struct * ctx;
}
fmpz_mod_poly_compose_mod_precomp_preinv_arg_t;


/*  Initialisation and memory management *************************************/

FMPZ_MOD_POLY_INLINE
void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    poly->coeffs = NULL;
    poly->alloc  = 0;
    poly->length = 0;
}

void fmpz_mod_poly_init2(fmpz_mod_poly_t poly, slong alloc,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_clear(fmpz_mod_poly_t poly,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_realloc(fmpz_mod_poly_t poly, slong alloc,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len);

FMPZ_MOD_POLY_INLINE
void fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len,
                                                      const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_fit_length(poly, len);
}


/*  Normalisation and truncation *********************************************/

FMPZ_MOD_POLY_INLINE
void _fmpz_mod_poly_normalise(fmpz_mod_poly_t poly)
{
    slong i;

    for (i = poly->length - 1; (i >= 0) && !poly->coeffs[i]; i--) ;
    poly->length = i + 1;
}

void _fmpz_mod_poly_set_length(fmpz_mod_poly_t poly, slong len);

void fmpz_mod_poly_truncate(fmpz_mod_poly_t poly, slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_set_trunc(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_canonical(const fmpz_mod_poly_t A, const fmpz_mod_ctx_t ctx);


/*  Randomisation ************************************************************/

void fmpz_mod_poly_randtest(fmpz_mod_poly_t f, flint_rand_t state,
                                          slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_randtest_irreducible(fmpz_mod_poly_t f,
                      flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_randtest_not_zero(fmpz_mod_poly_t f,
                      flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_randtest_monic(fmpz_mod_poly_t f,
                      flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_randtest_monic_irreducible(fmpz_mod_poly_t f,
                      flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_randtest_monic_primitive(fmpz_mod_poly_t f,
                      flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_randtest_trinomial(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_randtest_trinomial_irreducible(fmpz_mod_poly_t f,
                            flint_rand_t state, slong len, slong max_attempts,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_randtest_pentomial(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_randtest_pentomial_irreducible(fmpz_mod_poly_t f,
                            flint_rand_t state, slong len, slong max_attempts,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_randtest_sparse_irreducible(fmpz_mod_poly_t poly,
                      flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx);

/*  Attributes ***************************************************************/

slong fmpz_mod_poly_degree(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);
slong fmpz_mod_poly_length(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

fmpz * fmpz_mod_poly_lead(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_monic(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);
int fmpz_mod_poly_is_one(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);
int fmpz_mod_poly_is_gen(const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx);
int fmpz_mod_poly_is_unit(const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx);

/*  Assignment and basic manipulation ****************************************/

void fmpz_mod_poly_set(fmpz_mod_poly_t poly1,
                        const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_INLINE
void fmpz_mod_poly_swap(fmpz_mod_poly_t poly1,
                               fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)
{
    FLINT_SWAP(fmpz_mod_poly_struct, *poly1, *poly2);
}

void _fmpz_mod_poly_reverse(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_mod_poly_reverse(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_INLINE
void fmpz_mod_poly_zero(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
   _fmpz_mod_poly_set_length(poly, 0);
}

void fmpz_mod_poly_one(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_gen(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_zero_coeffs(fmpz_mod_poly_t poly, slong i, slong j, const fmpz_mod_ctx_t ctx);

ulong fmpz_mod_poly_deflation(const fmpz_mod_poly_t input, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_deflate(fmpz_mod_poly_t result, const fmpz_mod_poly_t input, ulong deflation, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_inflate(fmpz_mod_poly_t result, const fmpz_mod_poly_t input, ulong inflation, const fmpz_mod_ctx_t ctx);

/*  Conversion ***************************************************************/

void fmpz_mod_poly_set_ui(fmpz_mod_poly_t f, ulong x, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_set_nmod_poly(fmpz_mod_poly_t f, const nmod_poly_t g);
void fmpz_mod_poly_set_fmpz(fmpz_mod_poly_t poly, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f, const fmpz_poly_t g, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_get_nmod_poly(nmod_poly_t f, const fmpz_mod_poly_t g);
void fmpz_mod_poly_get_fmpz_poly(fmpz_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx);

/*  Comparison ***************************************************************/

int fmpz_mod_poly_equal(const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx);
int fmpz_mod_poly_equal_trunc(const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_is_zero(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

/*  Getting and setting coefficients *****************************************/

void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, slong n, const fmpz_t x, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_set_coeff_ui(fmpz_mod_poly_t poly, slong n, ulong x, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_set_coeff_si(fmpz_mod_poly_t poly, slong n, slong x, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly,
                                             slong n, const fmpz_mod_ctx_t ctx);

/*  Shifting *****************************************************************/

void _fmpz_mod_poly_shift_left(fmpz * res, const fmpz * poly, slong len, slong n);
void fmpz_mod_poly_shift_left(fmpz_mod_poly_t f, const fmpz_mod_poly_t g, slong n, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_shift_right(fmpz * res, const fmpz * poly, slong len, slong n);

void fmpz_mod_poly_shift_right(fmpz_mod_poly_t f, const fmpz_mod_poly_t g, slong n, const fmpz_mod_ctx_t ctx);

/*  Addition and subtraction *************************************************/

void _fmpz_mod_poly_neg(fmpz *res, const fmpz *poly, slong len, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_neg(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_add(fmpz *res, const fmpz *poly1, slong len1, const fmpz *poly2, slong len2, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_add(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_sub(fmpz *res, const fmpz *poly1, slong len1, const fmpz *poly2, slong len2, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_sub(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_add_series(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_sub_series(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx);

/*  Scalar multiplication ****************************************************/

void _fmpz_mod_poly_scalar_mul_fmpz(fmpz *res, const fmpz *poly, slong len, const fmpz_t x, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_scalar_mul_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_poly_scalar_mul_ui(fmpz *res, const fmpz *poly, slong len, ulong x, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_scalar_mul_ui(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong x, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_scalar_addmul_fmpz(fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_t x, const fmpz_mod_ctx_t ctx);

/*  Scalar division ****************************************************/

void _fmpz_mod_poly_scalar_div_fmpz(fmpz *res, const fmpz *poly, slong len, const fmpz_t x, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_scalar_div_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx);

/*  Multiplication ***********************************************************/

void _fmpz_mod_poly_mul(fmpz *res, const fmpz *poly1, slong len1, const fmpz *poly2, slong len2, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_mul(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_mullow(fmpz *res, const fmpz *poly1, slong len1, const fmpz *poly2, slong len2, slong n, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_mullow(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_sqr(fmpz *res, const fmpz *poly, slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_mulhigh(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong start, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_sqr(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_mulmod(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz * f, slong lenf, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_mulmod(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_mulmod_preinv(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz * f, slong lenf, const fmpz* finv, slong lenfinv, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_mulmod_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, const fmpz_mod_ctx_t ctx);

/*  Powering *****************************************************************/

void _fmpz_mod_poly_pow(fmpz *rop, const fmpz *op, slong len, ulong e, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_pow(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, ulong e, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_pow_trunc(fmpz * res, const fmpz * poly, ulong e, slong trunc, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_pow_trunc(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong e, slong trunc, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_poly_pow_trunc_binexp(fmpz * res, const fmpz * poly, ulong e, slong trunc, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_pow_trunc_binexp(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong e, slong trunc, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_powmod_ui_binexp(fmpz * res, const fmpz * poly, ulong e, const fmpz * f, slong lenf, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_powmod_ui_binexp(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong e, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz * res, const fmpz * poly,
                              ulong e, const fmpz * f, slong lenf,
                              const fmpz * finv, slong lenfinv, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, ulong e,
                         const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_powmod_fmpz_binexp(fmpz * res, const fmpz * poly,
                                  const fmpz_t e, const fmpz * f,
                                  slong lenf, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_powmod_fmpz_binexp(fmpz_mod_poly_t res,
                            const fmpz_mod_poly_t poly, const fmpz_t e,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz * res, const fmpz * poly,
                              const fmpz_t e, const fmpz * f, slong lenf,
                              const fmpz* finv, slong lenfinv, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz_mod_poly_t res,
                          const fmpz_mod_poly_t poly, const fmpz_t e,
                          const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz * res, const fmpz_t e,
                  const fmpz * f, slong lenf, const fmpz* finv, slong lenfinv,
                                                               const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz_mod_poly_t res,
         const fmpz_t e, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_powmod_linear_fmpz_preinv(fmpz_mod_poly_t res,
                      const fmpz_t a, const fmpz_t e, const fmpz_mod_poly_t f,
                         const fmpz_mod_poly_t finv, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_powers_mod_preinv_naive(fmpz ** res,
                    const fmpz * f, slong flen, slong n, const fmpz * g,
                 slong glen, const fmpz * ginv, slong ginvlen, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_powers_mod_naive(fmpz_mod_poly_struct * res,
                    const fmpz_mod_poly_t f, slong n, const fmpz_mod_poly_t g,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_powers_mod_preinv_threaded_pool(fmpz ** res,
       const fmpz * f, slong flen, slong n, const fmpz * g, slong glen,
                          const fmpz * ginv, slong ginvlen, const fmpz_mod_ctx_t ctx,
	                      thread_pool_handle * threads, slong num_threads);

void fmpz_mod_poly_powers_mod_bsgs(fmpz_mod_poly_struct * res,
                    const fmpz_mod_poly_t f, slong n, const fmpz_mod_poly_t g,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_frobenius_powers_2exp_precomp(
            fmpz_mod_poly_frobenius_powers_2exp_t pow, const fmpz_mod_poly_t f,
                const fmpz_mod_poly_t finv, ulong m, const fmpz_mod_ctx_t ctx);


void fmpz_mod_poly_frobenius_powers_2exp_clear(
          fmpz_mod_poly_frobenius_powers_2exp_t pow, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_frobenius_power(fmpz_mod_poly_t res,
                            fmpz_mod_poly_frobenius_powers_2exp_t pow,
                   const fmpz_mod_poly_t f, ulong m, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_frobenius_powers_precomp(
                fmpz_mod_poly_frobenius_powers_t pow, const fmpz_mod_poly_t f,
                const fmpz_mod_poly_t finv, ulong m, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_frobenius_powers_clear(
               fmpz_mod_poly_frobenius_powers_t pow, const fmpz_mod_ctx_t ctx);

/*  Division *****************************************************************/

void _fmpz_mod_poly_divrem_basecase(fmpz * Q, fmpz * R,
                    const fmpz * A, slong lenA, const fmpz * B, slong lenB,
                                            const fmpz_t invB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_divrem_basecase(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_div_newton_n_preinv (fmpz* Q, const fmpz* A,
                      slong lenA, const fmpz* B, slong lenB, const fmpz* Binv,
                                    slong lenBinv, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_div_newton_n_preinv(fmpz_mod_poly_t Q,
                         const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                         const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_divrem_newton_n_preinv (fmpz* Q, fmpz* R,
             const fmpz* A, slong lenA, const fmpz* B, slong lenB,
                            const fmpz* Binv, slong lenBinv, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q,
          fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                         const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx);

ulong fmpz_mod_poly_remove(fmpz_mod_poly_t f,
                            const fmpz_mod_poly_t p, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_rem_basecase(fmpz * R,
                        const fmpz * A, slong lenA, const fmpz * B, slong lenB,
                                            const fmpz_t invB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_rem_basecase(fmpz_mod_poly_t R,
                            const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_divrem(fmpz *Q, fmpz *R,
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                           const fmpz_t invB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_divrem(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_div(fmpz *Q,
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                           const fmpz_t invB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_div(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                              const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_divrem_f(fmpz_t f, fmpz *Q, fmpz *R,
                             const fmpz *A, slong lenA,
                             const fmpz *B, slong lenB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_divrem_f(fmpz_t f, fmpz_mod_poly_t Q,
          fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_rem(fmpz *R, const fmpz *A, slong lenA, const fmpz *B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_rem(fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_rem_f(fmpz_t f, fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx);

/* Divisibility testing ******************************************************/

int _fmpz_mod_poly_divides_classical(fmpz * Q, const fmpz * A,
             slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_divides_classical(fmpz_mod_poly_t Q,
                        const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

int _fmpz_mod_poly_divides(fmpz * Q, const fmpz * A, slong lenA,
                         const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_divides(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                            const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx);

/*  Power series inversion ***************************************************/

void _fmpz_mod_poly_inv_series(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_inv_series(fmpz_mod_poly_t Qinv, const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx);

void
fmpz_mod_poly_inv_series_f(fmpz_t f, fmpz_mod_poly_t Qinv,
                    const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_div_series(fmpz * Q, const fmpz * A, slong Alen,
                      const fmpz * B, slong Blen, slong n, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_div_series(fmpz_mod_poly_t Q,
                    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, slong n,
                                                     const fmpz_mod_ctx_t ctx);

/*  Greatest common divisor **************************************************/

void fmpz_mod_poly_make_monic(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_make_monic_f(fmpz_t f, fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_gcd(fmpz *G, const fmpz *A, slong lenA,
                                           const fmpz *B, slong lenB,
                                           const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_gcd(fmpz_mod_poly_t G,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz *G,
                                    const fmpz *A, slong lenA,
                                    const fmpz *B, slong lenB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_INLINE
slong _fmpz_mod_poly_gcd_f(fmpz_t f, fmpz *G,
                          const fmpz *A, slong lenA,
                          const fmpz *B, slong lenB, const fmpz_mod_ctx_t ctx)
{
    return _fmpz_mod_poly_gcd_euclidean_f(f, G, A, lenA, B, lenB, ctx);
}

FMPZ_MOD_POLY_INLINE
void fmpz_mod_poly_gcd_f(fmpz_t f, fmpz_mod_poly_t G, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_gcd_euclidean_f(f, G, A, B, ctx);
}

slong _fmpz_mod_poly_hgcd(fmpz **M, slong *lenM,
                     fmpz *A, slong *lenA, fmpz *B, slong *lenB,
                     const fmpz *a, slong lena, const fmpz *b, slong lenb,
                     const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_xgcd_euclidean_f(fmpz_t f, fmpz *G, fmpz *S, fmpz *T,
                                   const fmpz *A, slong lenA,
                                   const fmpz *B, slong lenB,
                                   const fmpz_t invB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_xgcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G,
                             fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_xgcd(fmpz *G, fmpz *S, fmpz *T,
                          const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                          const fmpz_t invB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_xgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S,
          fmpz_mod_poly_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_INLINE slong
_fmpz_mod_poly_xgcd_f(fmpz_t f, fmpz *G, fmpz *S, fmpz *T,
                    const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                    const fmpz_t invB, const fmpz_mod_ctx_t ctx)
{
    return _fmpz_mod_poly_xgcd_euclidean_f(f, G, S, T, A, lenA, B, lenB, invB, ctx);
}

FMPZ_MOD_POLY_INLINE void
fmpz_mod_poly_xgcd_f(fmpz_t f, fmpz_mod_poly_t G, fmpz_mod_poly_t S,
          fmpz_mod_poly_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_xgcd_euclidean_f(f, G, S, T, A, B, ctx);
}

slong _fmpz_mod_poly_gcdinv_euclidean_f(fmpz_t f, fmpz *G, fmpz *S,
                    const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                                            const fmpz_t invA, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_gcdinv_euclidean_f(fmpz_t f, fmpz_mod_poly_t G,
          fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_gcdinv_euclidean(fmpz *G, fmpz *S,
                  const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                                            const fmpz_t invA, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_gcdinv_euclidean(fmpz_mod_poly_t G,
          fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_gcdinv(fmpz *G, fmpz *S,
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                           const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_gcdinv(fmpz_mod_poly_t G, fmpz_mod_poly_t S,
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_gcdinv_f(fmpz_t f, fmpz *G, fmpz *S,
                           const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                           const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_gcdinv_f(fmpz_t f, fmpz_mod_poly_t G,
          fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx);

int _fmpz_mod_poly_invmod(fmpz *A,
                          const fmpz *B, slong lenB,
                          const fmpz *P, slong lenP, const fmpz_mod_ctx_t ctx);

int _fmpz_mod_poly_invmod_f(fmpz_t f, fmpz *A,
                          const fmpz *B, slong lenB,
                          const fmpz *P, slong lenP, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_invmod(fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                            const fmpz_mod_poly_t P, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_invmod_f(fmpz_t f, fmpz_mod_poly_t A,
                         const fmpz_mod_poly_t B, const fmpz_mod_poly_t P,
                                                     const fmpz_mod_ctx_t ctx);

/* Square root ***************************************************************/

void _fmpz_mod_poly_invsqrt_series(fmpz * g,
                                  const fmpz * h, slong hlen, slong n, const fmpz_mod_ctx_t mod);

void fmpz_mod_poly_invsqrt_series(fmpz_mod_poly_t g,
                         const fmpz_mod_poly_t h, slong n, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_sqrt_series(fmpz * g,
                                  const fmpz * h, slong hlen, slong n, const fmpz_mod_ctx_t mod);

void fmpz_mod_poly_sqrt_series(fmpz_mod_poly_t g,
                         const fmpz_mod_poly_t h, slong n, const fmpz_mod_ctx_t ctx);

int _fmpz_mod_poly_sqrt(fmpz * s,
                                const fmpz * p, slong len, const fmpz_mod_ctx_t mod);

int fmpz_mod_poly_sqrt(fmpz_mod_poly_t b,
                                  const fmpz_mod_poly_t a, const fmpz_mod_ctx_t ctx);

/*  Minpoly  *****************************************************************/

slong _fmpz_mod_poly_minpoly_bm(fmpz * poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_minpoly_bm(fmpz_mod_poly_t poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_minpoly_hgcd(fmpz * poly, const fmpz* seq, slong len, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_minpoly_hgcd(fmpz_mod_poly_t poly, const fmpz* seq, slong len, const fmpz_mod_ctx_t ctx);

slong _fmpz_mod_poly_minpoly(fmpz* poly, const fmpz* seq, slong len, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_minpoly(fmpz_mod_poly_t poly, const fmpz* seq, slong len, const fmpz_mod_ctx_t ctx);

/*  Resultant  ***************************************************************/

void _fmpz_mod_poly_resultant_euclidean(fmpz_t res,
                                    const fmpz *poly1, slong len1,
                              const fmpz *poly2, slong len2, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_resultant_euclidean(fmpz_t r,
                            const fmpz_mod_poly_t f, const fmpz_mod_poly_t g,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz *A, slong lenA,
                                  const fmpz *B, slong lenB, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz_mod_poly_t A,
                            const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_INLINE void
_fmpz_mod_poly_resultant(fmpz_t res, const fmpz *poly1, slong len1,
                     const fmpz *poly2, slong len2, const fmpz_mod_ctx_t ctx)
{
    if (len1 < FMPZ_MOD_POLY_GCD_CUTOFF)
        _fmpz_mod_poly_resultant_euclidean(res, poly1, len1, poly2, len2, ctx);
    else
        _fmpz_mod_poly_resultant_hgcd(res, poly1, len1, poly2, len2, ctx);
}

FMPZ_MOD_POLY_INLINE void
fmpz_mod_poly_resultant(fmpz_t res, const fmpz_mod_poly_t f,
                             const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)
{
    if (FLINT_MAX(f->length, g->length) < FMPZ_MOD_POLY_GCD_CUTOFF)
       fmpz_mod_poly_resultant_euclidean(res, f, g, ctx);
    else
       fmpz_mod_poly_resultant_hgcd(res, f, g, ctx);
}

/*  Discriminant  ************************************************************/

void _fmpz_mod_poly_discriminant(fmpz_t d, const fmpz *poly, slong len, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_discriminant(fmpz_t d, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

/*  Derivative  **************************************************************/

void _fmpz_mod_poly_derivative(fmpz *res, const fmpz *poly, slong len, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_derivative(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

/*  Evaluation  **************************************************************/

void _fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz *poly, slong len, const fmpz_t a, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz_mod_poly_t poly, const fmpz_t a, const fmpz_mod_ctx_t ctx);

fmpz_poly_struct ** _fmpz_mod_poly_tree_alloc(slong len);
void _fmpz_mod_poly_tree_free(fmpz_poly_struct ** tree, slong len);
void _fmpz_mod_poly_tree_build(fmpz_poly_struct ** tree, const fmpz * roots, slong len, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys, const fmpz * coeffs, slong len, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys, const fmpz_mod_poly_t poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(fmpz * vs, const fmpz * poly, slong plen, fmpz_poly_struct * const * tree, slong len, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, const fmpz * poly, slong plen, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, const fmpz_mod_poly_t poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys, const fmpz * coeffs, slong len, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx);
void fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys, const fmpz_mod_poly_t poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx);

/*  Composition  *************************************************************/

void _fmpz_mod_poly_compose(fmpz *res, const fmpz *poly1, slong len1,
                                              const fmpz *poly2, slong len2,
                                              const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_compose(fmpz_mod_poly_t res,
                    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                                     const fmpz_mod_ctx_t ctx);

/* Modular composition  ******************************************************/

void _fmpz_mod_poly_compose_mod(fmpz * res, const fmpz * f, slong lenf, const fmpz * g,
                                       const fmpz * h, slong lenh, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_compose_mod(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                        const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_compose_mod_brent_kung(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, const fmpz * poly3, slong len3, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_compose_mod_brent_kung(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                        const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_reduce_matrix_mod_poly (fmpz_mat_t A,
        const fmpz_mat_t B, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_precompute_matrix(fmpz_mat_t A, const fmpz * poly1,
                          const fmpz * poly2, slong len2, const fmpz * poly2inv,
                          slong len2inv, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_precompute_matrix_worker(void * arg_ptr);

void fmpz_mod_poly_precompute_matrix(fmpz_mat_t A,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                     const fmpz_mod_poly_t poly2inv, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz * res,
         const fmpz * poly1, slong len1, const fmpz_mat_t A, const fmpz * poly3,
         slong len3, const fmpz * poly3inv, slong len3inv, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker(void * arg_ptr);

void fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mat_t A,
                   const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz * res, const fmpz * poly1,
                 slong len1, const fmpz * poly2, const fmpz * poly3, slong len3,
                 const fmpz * poly3inv, slong len3inv, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                   const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_compose_mod_horner(fmpz * res, const fmpz * f, slong lenf, const fmpz * g,
                                              const fmpz * h, slong lenh, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_compose_mod_horner(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                        const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(fmpz_mod_poly_struct * res,
                 const fmpz_mod_poly_struct * polys, slong len1, slong l,
                 const fmpz * g, slong glen, const fmpz * poly, slong len,
		 const fmpz * polyinv, slong leninv, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(
      fmpz_mod_poly_struct * res, const fmpz_mod_poly_struct * polys,
      slong len1,slong n, const fmpz_mod_poly_t g, const fmpz_mod_poly_t poly,
                      const fmpz_mod_poly_t polyinv, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(fmpz_mod_poly_struct * res,
             const fmpz_mod_poly_struct * polys, slong lenpolys, slong l,
             const fmpz * g, slong glen, const fmpz * poly, slong len,
             const fmpz * polyinv, slong leninv, const fmpz_mod_ctx_t ctx,
                              thread_pool_handle * threads, slong num_threads);

void fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(fmpz_mod_poly_struct * res,
            const fmpz_mod_poly_struct * polys, slong len1, slong n,
            const fmpz_mod_poly_t g, const fmpz_mod_poly_t poly,
            const fmpz_mod_poly_t polyinv, const fmpz_mod_ctx_t ctx,
                              thread_pool_handle * threads, slong num_threads);

void
fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(fmpz_mod_poly_struct * res,
                    const fmpz_mod_poly_struct * polys, slong len1, slong n,
                    const fmpz_mod_poly_t g, const fmpz_mod_poly_t poly,
                      const fmpz_mod_poly_t polyinv, const fmpz_mod_ctx_t ctx);

/* Norms *********************************************************************/

slong fmpz_mod_poly_hamming_weight(const fmpz_mod_poly_t A, const fmpz_mod_ctx_t ctx);

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
                    const fmpz_t invL, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_radix_init(fmpz_mod_poly_radix_t D,
                const fmpz_mod_poly_t R, slong degF, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_radix_clear(fmpz_mod_poly_radix_t D);

void _fmpz_mod_poly_radix(fmpz **B, const fmpz *F, fmpz **Rpow, fmpz **Rinv,
                          slong degR, slong k, slong i, fmpz *W, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_radix(fmpz_mod_poly_struct **B, const fmpz_mod_poly_t F,
                      const fmpz_mod_poly_radix_t D, const fmpz_mod_ctx_t ctx);

/*  Input and output *********************************************************/

char * fmpz_mod_poly_get_str(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);
char * fmpz_mod_poly_get_str_pretty(const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx);

#ifdef FLINT_HAVE_FILE
int _fmpz_mod_poly_fprint(FILE * file, const fmpz *poly, slong len, const fmpz_t p);
int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);
int fmpz_mod_poly_fprint_pretty(FILE * file, const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_fread(FILE * file, fmpz_mod_poly_t poly, fmpz_mod_ctx_t ctx);
#endif

int _fmpz_mod_poly_print(const fmpz *poly, slong len, const fmpz_t p);
int fmpz_mod_poly_print(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);
int fmpz_mod_poly_print_pretty(const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx);

/* Products *****************************************************************/

void _fmpz_mod_poly_product_roots_fmpz_vec(fmpz * poly,
                                   const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_product_roots_fmpz_vec(fmpz_mod_poly_t poly,
                           const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx);

int fmpz_mod_poly_find_distinct_nonzero_roots(fmpz * roots,
                            const fmpz_mod_poly_t P, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_split_rabin(fmpz_mod_poly_t a, fmpz_mod_poly_t b,
               const fmpz_mod_poly_t f, const fmpz_t halfp, fmpz_mod_poly_t t,
         fmpz_mod_poly_t t2, flint_rand_t randstate, const fmpz_mod_ctx_t ctx);

/* Characteristic polynomial and minimal polynomial */

void fmpz_mod_mat_charpoly_berkowitz(fmpz_mod_poly_t p,
                             const fmpz_mod_mat_t M, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_INLINE void fmpz_mod_mat_charpoly(fmpz_mod_poly_t p,
                              const fmpz_mod_mat_t M, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_mat_charpoly_berkowitz(p, M, ctx);
}


void fmpz_mod_mat_minpoly(fmpz_mod_poly_t p, const fmpz_mod_mat_t M,
                                                     const fmpz_mod_ctx_t ctx);

/* Berlekamp-Massey Algorithm - see fmpz_mod_poly/berlekamp_massey.c for more info ********/
typedef struct {
    slong npoints;
    fmpz_mod_poly_t R0, R1;
    fmpz_mod_poly_t V0, V1;
    fmpz_mod_poly_t qt, rt;
    fmpz_mod_poly_t points;
} fmpz_mod_berlekamp_massey_struct;
typedef fmpz_mod_berlekamp_massey_struct fmpz_mod_berlekamp_massey_t[1];

void fmpz_mod_berlekamp_massey_init(
                      fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_berlekamp_massey_start_over(
                      fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_berlekamp_massey_clear(fmpz_mod_berlekamp_massey_t B,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_berlekamp_massey_print(
                const fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_berlekamp_massey_add_points(
                   fmpz_mod_berlekamp_massey_t B, const fmpz * a, slong count,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_berlekamp_massey_add_zeros(
                                   fmpz_mod_berlekamp_massey_t B, slong count,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_berlekamp_massey_add_point(
                                fmpz_mod_berlekamp_massey_t B, const fmpz_t a,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_berlekamp_massey_add_point_ui(
                                      fmpz_mod_berlekamp_massey_t B, ulong a,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_berlekamp_massey_reduce(
                      fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_INLINE const fmpz * fmpz_mod_berlekamp_massey_points(
                        const fmpz_mod_berlekamp_massey_t B)
{
    return B->points->coeffs;
}

FMPZ_MOD_POLY_INLINE slong fmpz_mod_berlekamp_massey_point_count(
                        const fmpz_mod_berlekamp_massey_t B)
{
    return B->points->length;
}

FMPZ_MOD_POLY_INLINE const fmpz_mod_poly_struct * fmpz_mod_berlekamp_massey_V_poly(
                        const fmpz_mod_berlekamp_massey_t B)
{
    return B->V1;
}

FMPZ_MOD_POLY_INLINE const fmpz_mod_poly_struct * fmpz_mod_berlekamp_massey_R_poly(
                        const fmpz_mod_berlekamp_massey_t B)
{
    return B->R1;
}

/* Inlines *******************************************************************/

void fmpz_mod_poly_add_si(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, slong c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_sub_si(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, slong c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_si_sub(fmpz_mod_poly_t res, slong c,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_add_fmpz(fmpz_mod_poly_t res,
         const fmpz_mod_poly_t poly, const fmpz_t c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_sub_fmpz(fmpz_mod_poly_t res,
         const fmpz_mod_poly_t poly, const fmpz_t c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_fmpz_sub(fmpz_mod_poly_t res, const fmpz_t c,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

/* Declare dead functions *******************************************************************/

#define fmpz_mod_poly_set_coeff_mpz _Pragma("GCC error \"'fmpz_mod_poly_set_coeff_mpz' is deprecated. Use 'fmpz_mod_poly_set_coeff_fmpz' instead.\"")
#define fmpz_mod_poly_get_coeff_mpz _Pragma("GCC error \"'fmpz_mod_poly_get_coeff_mpz' is deprecated. Use 'fmpz_mod_poly_get_coeff_fmpz' instead.\"")

#ifdef __cplusplus
}
#endif

#endif
