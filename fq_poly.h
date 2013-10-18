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

    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#ifndef FQ_POLY_H
#define FQ_POLY_H

#undef ulong                /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include "fq.h"
#include "fq_vec.h"

/*  Type definitions *********************************************************/

typedef struct
{
    fq_struct *coeffs;
    long alloc;
    long length;
}
fq_poly_struct;

typedef fq_poly_struct fq_poly_t[1];

/*  Memory management ********************************************************/

void fq_poly_init(fq_poly_t poly, const fq_ctx_t ctx);

void fq_poly_init2(fq_poly_t poly, long alloc, const fq_ctx_t ctx);

void fq_poly_realloc(fq_poly_t poly, long alloc, const fq_ctx_t ctx);

void fq_poly_truncate(fq_poly_t poly, long len, const fq_ctx_t ctx);

void fq_poly_fit_length(fq_poly_t poly, long len, const fq_ctx_t ctx);

void fq_poly_clear(fq_poly_t poly, const fq_ctx_t ctx);

void _fq_poly_normalise(fq_poly_t poly, const fq_ctx_t ctx);

void _fq_poly_normalise2(fq_struct *poly, long *length, const fq_ctx_t ctx);

static __inline__ void _fq_poly_set_length(fq_poly_t poly, long len, const fq_ctx_t ctx)
{
    if (poly->length > len)
    {
        long i;

        for (i = len; i < poly->length; i++)
            fq_zero(poly->coeffs + i, ctx);
    }
    poly->length = len;
}

#define FQ_VEC_NORM(vec, i, ctx)                    \
do {                                                \
    while ((i) && fq_is_zero((vec) + (i) - 1, ctx)) \
        (i)--;                                      \
} while (0)

/*  Polynomial parameters  ***************************************************/

static __inline__ long fq_poly_length(const fq_poly_t poly, const fq_ctx_t ctx)
{
    return poly->length;
}

static __inline__ long fq_poly_degree(const fq_poly_t poly, const fq_ctx_t ctx)
{
    return poly->length - 1;
}

static __inline__ fq_struct * fq_poly_lead(const fq_poly_t poly, const fq_ctx_t ctx)
{
    return poly->length > 0 ? poly->coeffs + (poly->length - 1) : NULL;
}

/*  Randomisation  ***********************************************************/

void fq_poly_randtest(fq_poly_t f, flint_rand_t state, 
                      long len, const fq_ctx_t ctx);

void fq_poly_randtest_not_zero(fq_poly_t f, flint_rand_t state, 
                               long len, const fq_ctx_t ctx);

void fq_poly_randtest_monic (fq_poly_t f, flint_rand_t state,
                             long len, const fq_ctx_t ctx);

void fq_poly_randtest_irreducible (fq_poly_t f, flint_rand_t state,
                                   long len, const fq_ctx_t ctx);

/*  Factoring ****************************************************************/
typedef struct
{
    fq_poly_struct *poly;
    slong *exp;
    slong num;
    slong alloc;
} fq_poly_factor_struct;

typedef fq_poly_factor_struct fq_poly_factor_t[1];


void fq_poly_factor_init(fq_poly_factor_t fac, const fq_ctx_t ctx);

void fq_poly_factor_clear(fq_poly_factor_t fac, const fq_ctx_t ctx);

void fq_poly_factor_realloc(fq_poly_factor_t fac, slong alloc,
                            const fq_ctx_t ctx);

void fq_poly_factor_fit_length(fq_poly_factor_t fac, slong len,
                               const fq_ctx_t ctx);

void fq_poly_factor_set(fq_poly_factor_t res, const fq_poly_factor_t fac,
                        const fq_ctx_t ctx);

void fq_poly_factor_insert(fq_poly_factor_t fac, const fq_poly_t poly, 
                           slong exp, const fq_ctx_t ctx);

void fq_poly_factor_print(const fq_poly_factor_t fac, const fq_ctx_t ctx);

void
fq_poly_factor_print_pretty(const fq_poly_factor_t fac, const char * var,
                            const fq_ctx_t ctx);

void fq_poly_factor_concat(fq_poly_factor_t res, const fq_poly_factor_t fac,
                           const fq_ctx_t ctx);

void fq_poly_factor_pow(fq_poly_factor_t fac, slong exp, const fq_ctx_t ctx);

int
_fq_poly_is_squarefree(const fq_struct * f, slong len, const fq_ctx_t ctx);

int
fq_poly_is_squarefree(const fq_poly_t f, const fq_ctx_t ctx);

void
fq_poly_factor_squarefree(fq_poly_factor_t res, const fq_poly_t f,
                          const fq_ctx_t ctx);

int
fq_poly_is_irreducible(const fq_poly_t f, const fq_ctx_t ctx);

int 
fq_poly_is_irreducible_ddf(const fq_poly_t f, const fq_ctx_t ctx);

void
fq_poly_factor_distinct_deg(fq_poly_factor_t res, const fq_poly_t poly, 
                            slong * const *degs, const fq_ctx_t ctx);

int
fq_poly_factor_equal_deg_prob(fq_poly_t factor, flint_rand_t state,
                              const fq_poly_t pol, slong d,
                              const fq_ctx_t ctx);

void
fq_poly_factor_equal_deg(fq_poly_factor_t factors, const fq_poly_t pol, 
                         slong d, const fq_ctx_t ctx);

void
fq_poly_factor_cantor_zassenhaus(fq_poly_factor_t res, const fq_poly_t f,
                                 const fq_ctx_t ctx);

void
fq_poly_factor_kaltofen_shoup(fq_poly_factor_t res, const fq_poly_t poly,
                              const fq_ctx_t ctx);

void
fq_poly_factor_berlekamp(fq_poly_factor_t factors, const fq_poly_t f,
                         const fq_ctx_t ctx);

void
fq_poly_factor_with_berlekamp(fq_poly_factor_t result, fq_t leading_coeff,
                              const fq_poly_t input, const fq_ctx_t ctx);

void
fq_poly_factor_with_cantor_zassenhaus(fq_poly_factor_t result, fq_t leading_coeff,
                                      const fq_poly_t input, const fq_ctx_t ctx);

void
fq_poly_factor_with_kaltofen_shoup(fq_poly_factor_t result, fq_t leading_coeff,
                                   const fq_poly_t input, const fq_ctx_t ctx);

void
fq_poly_factor(fq_poly_factor_t result, fq_t leading_coeff,
               const fq_poly_t input, const fq_ctx_t ctx);


/*  Assignment and basic manipulation  ***************************************/

void _fq_poly_set(fq_struct *rop, const fq_struct *op, long len, const fq_ctx_t ctx);

void fq_poly_set(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx);

void fq_poly_set_fq(fq_poly_t poly, const fq_t c, const fq_ctx_t ctx);

void fq_poly_swap(fq_poly_t op1, fq_poly_t op2, const fq_ctx_t ctx);

static __inline__ void _fq_poly_zero(fq_struct *rop, long len, const fq_ctx_t ctx)
{
    long i;

    for (i = 0; i < len; i++)
        fq_zero(rop + i, ctx);
}

static __inline__ void fq_poly_zero(fq_poly_t poly, const fq_ctx_t ctx)
{
    _fq_poly_set_length(poly, 0, ctx);
}

void fq_poly_one(fq_poly_t poly, const fq_ctx_t ctx);

void fq_poly_gen(fq_poly_t f, const fq_ctx_t ctx);

void _fq_poly_make_monic(fq_struct *rop, const fq_struct *op, long length, const fq_ctx_t ctx);

void fq_poly_make_monic(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx);

void
_fq_poly_reverse(fq_struct * res, const fq_struct * poly, slong len, slong n,
                 const fq_ctx_t ctx);

void
fq_poly_reverse(fq_poly_t res, const fq_poly_t poly, slong n,
                const fq_ctx_t ctx);

ulong
fq_poly_deflation(const fq_poly_t input, const fq_ctx_t ctx);

void
fq_poly_deflate(fq_poly_t result, const fq_poly_t input, ulong deflation,
                const fq_ctx_t ctx);

void
fq_poly_inflate(fq_poly_t result, const fq_poly_t input, ulong inflation,
                const fq_ctx_t ctx);

/*  Getting and setting coefficients  ****************************************/

void fq_poly_get_coeff(fq_t x, const fq_poly_t poly, long n, const fq_ctx_t ctx);

void fq_poly_set_coeff(fq_poly_t poly, long n, const fq_t x, const fq_ctx_t ctx);

static __inline__
void fq_poly_set_coeff_fmpz(fq_poly_t poly, long n, const fmpz_t x, const fq_ctx_t ctx)
{
    if (fmpz_is_zero(x))
    {
        fq_poly_zero(poly, ctx);
    }
    else
    {
        fq_poly_fit_length(poly, 1, ctx);
        fq_set_fmpz(poly->coeffs, x, ctx);
        _fq_poly_set_length(poly, 1, ctx);
    }
}

static __inline__ int 
fq_poly_is_gen(const fq_poly_t poly,  const fq_ctx_t ctx)
{
    return ((poly->length == 2) &&
            fq_is_zero(poly->coeffs, ctx) &&
            fq_is_one(poly->coeffs + 1, ctx));
}

/*  Comparison  **************************************************************/

int fq_poly_equal(const fq_poly_t poly1, const fq_poly_t poly2, const fq_ctx_t ctx);

static __inline__ int fq_poly_is_zero(const fq_poly_t poly, const fq_ctx_t ctx)
{
    return (poly->length == 0);
}

static __inline__ int fq_poly_is_one(const fq_poly_t op, const fq_ctx_t ctx)
{
    return (op->length == 1) && (fq_is_one(op->coeffs + 0, ctx));
}

static __inline__ int fq_poly_is_unit(const fq_poly_t op, const fq_ctx_t ctx)
{
    return (op->length == 1) && (!(fq_is_zero(op->coeffs + 0, ctx)));
}

static __inline__ int fq_poly_equal_fq(const fq_poly_t poly, const fq_t c, const fq_ctx_t ctx)
{
    return ((poly->length == 0) && fq_is_zero(c, ctx)) ||
        ((poly->length == 1) && fq_equal(poly->coeffs, c, ctx));
}

/*  Addition and subtraction  ************************************************/

void _fq_poly_add(fq_struct *res, 
                  const fq_struct *poly1, long len1, 
                  const fq_struct *poly2, long len2, 
                  const fq_ctx_t ctx);

void fq_poly_add(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, 
                 const fq_ctx_t ctx);

void _fq_poly_sub(fq_struct *res, 
                  const fq_struct *poly1, long len1, 
                  const fq_struct *poly2, long len2, 
                  const fq_ctx_t ctx);

void fq_poly_sub(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, 
                 const fq_ctx_t ctx);

void _fq_poly_neg(fq_struct *rop, const fq_struct *op, long len, 
                  const fq_ctx_t ctx);

void fq_poly_neg(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx);

/*  Scalar multiplication and division  **************************************/

void _fq_poly_scalar_mul_fq(fq_struct *rop, 
    const fq_struct *op, long len, const fq_t x, const fq_ctx_t ctx);

void fq_poly_scalar_mul_fq(fq_poly_t rop, 
    const fq_poly_t op, const fq_t x, const fq_ctx_t ctx);

void _fq_poly_scalar_addmul_fq(fq_struct *rop, 
    const fq_struct *op, long len, const fq_t x, const fq_ctx_t ctx);

void fq_poly_scalar_addmul_fq(fq_poly_t rop, 
    const fq_poly_t op, const fq_t x, const fq_ctx_t ctx);

void _fq_poly_scalar_submul_fq(fq_struct *rop, 
    const fq_struct *op, long len, const fq_t x, const fq_ctx_t ctx);

void fq_poly_scalar_submul_fq(fq_poly_t rop, 
    const fq_poly_t op, const fq_t x, const fq_ctx_t ctx);

/*  Multiplication  **********************************************************/

void _fq_poly_mul_classical(fq_struct *rop, 
                            const fq_struct *op1, long len1, 
                            const fq_struct *op2, long len2, 
                            const fq_ctx_t ctx);

void fq_poly_mul_classical(fq_poly_t rop, 
                           const fq_poly_t op1, const fq_poly_t op2, 
                           const fq_ctx_t ctx);

void _fq_poly_mul_reorder(fq_struct *rop, 
                          const fq_struct *op1, long len1, 
                          const fq_struct *op2, long len2, 
                          const fq_ctx_t ctx);

void fq_poly_mul_reorder(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2,
                         const fq_ctx_t ctx);

void _fq_poly_mul_KS(fq_struct *rop, const fq_struct *op1, long len1, 
                                     const fq_struct *op2, long len2, 
                                     const fq_ctx_t ctx);

void fq_poly_mul_KS(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, 
                    const fq_ctx_t ctx);

void _fq_poly_mul(fq_struct *rop, const fq_struct *op1, long len1, 
                                  const fq_struct *op2, long len2, 
                                  const fq_ctx_t ctx);

void fq_poly_mul(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2,
                 const fq_ctx_t ctx);

void _fq_poly_mullow_classical(fq_struct *rop, 
                               const fq_struct *op1, long len1, 
                               const fq_struct *op2, long len2, long n, 
                               const fq_ctx_t ctx);

void fq_poly_mullow_classical(fq_poly_t rop, 
    const fq_poly_t op1, const fq_poly_t op2, long n, const fq_ctx_t ctx);

void _fq_poly_mullow_KS(fq_struct *rop, 
                        const fq_struct *op1, long len1, 
                        const fq_struct *op2, long len2, long n, 
                        const fq_ctx_t ctx);

void fq_poly_mullow_KS(fq_poly_t rop, 
                       const fq_poly_t op1, const fq_poly_t op2, long n, 
                       const fq_ctx_t ctx);

void _fq_poly_mullow(fq_struct *rop, 
                     const fq_struct *op1, long len1, 
                     const fq_struct *op2, long len2, long n, 
                     const fq_ctx_t ctx);

void fq_poly_mullow(fq_poly_t rop, 
                    const fq_poly_t op1, const fq_poly_t op2, long n, 
                    const fq_ctx_t ctx);

void _fq_poly_mulmod(fq_struct * res, 
                     const fq_struct * poly1, slong len1,
                     const fq_struct * poly2, slong len2, 
                     const fq_struct * f, slong lenf, 
                     const fq_ctx_t ctx);

void fq_poly_mulmod(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2,
                    const fq_poly_t f, const fq_ctx_t ctx);

void _fq_poly_mulmod_preinv(fq_struct *res, const fq_struct *poly1, slong len1,
                            const fq_struct  *poly2, slong len2,
                            const fq_struct *f, slong lenf,
                            const fq_struct *finv, slong lenfinv,
                            const fq_ctx_t ctx);
void
fq_poly_mulmod_preinv(fq_poly_t res, const fq_poly_t poly1,
                      const fq_poly_t poly2, const fq_poly_t f,
                      const fq_poly_t finv, const fq_ctx_t ctx);

/* Squaring ******************************************************************/

void _fq_poly_sqr_classical(fq_struct *rop, 
                            const fq_struct *op, long len, 
                            const fq_ctx_t ctx);

void fq_poly_sqr_classical(fq_poly_t rop, 
                           const fq_poly_t op, const fq_ctx_t ctx);

void _fq_poly_sqr_reorder(fq_struct *rop, 
                          const fq_struct *op, long len, const fq_ctx_t ctx);

void fq_poly_sqr_reorder(fq_poly_t rop, 
                         const fq_poly_t op, const fq_ctx_t ctx);

void _fq_poly_sqr_KS(fq_struct *rop, const fq_struct *op, long len, 
                                     const fq_ctx_t ctx);

void fq_poly_sqr_KS(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx);

void _fq_poly_sqr(fq_struct *rop, const fq_struct *op, long len, 
                                  const fq_ctx_t ctx);

void fq_poly_sqr(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx);

/*  Powering  ****************************************************************/

void _fq_poly_pow(fq_struct *rop, const fq_struct *op, long len, ulong e, 
                                  const fq_ctx_t ctx);

void fq_poly_pow(fq_poly_t rop, const fq_poly_t op, ulong e, 
                 const fq_ctx_t ctx);

void
_fq_poly_powmod_fmpz_binexp(fq_struct * res, const fq_struct * poly,
                            const fmpz_t e, const fq_struct * f, slong lenf, 
                            const fq_ctx_t ctx);

void
fq_poly_powmod_fmpz_binexp(fq_poly_t res, const fq_poly_t poly, const fmpz_t e,
                           const fq_poly_t f, const fq_ctx_t ctx);

void
_fq_poly_powmod_fmpz_binexp_preinv(fq_struct * res, const fq_struct * poly,
                                   const fmpz_t e, const fq_struct * f, slong lenf, 
                                   const fq_struct *finv, slong lenfinv,
                                   const fq_ctx_t ctx);

void
fq_poly_powmod_fmpz_binexp_preinv(fq_poly_t res, const fq_poly_t poly, const fmpz_t e,
                                  const fq_poly_t f, const fq_poly_t finv,
                                  const fq_ctx_t ctx);

void
_fq_poly_powmod_ui_binexp(fq_struct * res, const fq_struct * poly, ulong e, 
                          const fq_struct * f, slong lenf, const fq_ctx_t ctx);

void
fq_poly_powmod_ui_binexp(fq_poly_t res, const fq_poly_t poly, ulong e,
                         const fq_poly_t f, const fq_ctx_t ctx);

void
_fq_poly_powmod_ui_binexp_preinv(fq_struct * res, const fq_struct * poly, ulong e, 
                                 const fq_struct * f, slong lenf,
                                 const fq_struct * finv, slong lenfinv,
                                 const fq_ctx_t ctx);

void
fq_poly_powmod_ui_binexp_preinv(fq_poly_t res, const fq_poly_t poly, ulong e,
                                const fq_poly_t f, const fq_poly_t finv,
                                const fq_ctx_t ctx);

/*  Shifting  ****************************************************************/

void _fq_poly_shift_left(fq_struct *rop, const fq_struct *op, long len, long n, const fq_ctx_t ctx);

void fq_poly_shift_left(fq_poly_t rop, const fq_poly_t op, long n, const fq_ctx_t ctx);

void _fq_poly_shift_right(fq_struct *rop, const fq_struct *op, long len, long n, const fq_ctx_t ctx);

void fq_poly_shift_right(fq_poly_t rop, const fq_poly_t op, long n, const fq_ctx_t ctx);

/*  Norms  *******************************************************************/

long _fq_poly_hamming_weight(const fq_struct *op, long len, const fq_ctx_t ctx);

long fq_poly_hamming_weight(const fq_poly_t op, const fq_ctx_t ctx);

/*  Greatest common divisor  *************************************************/

void fq_poly_gcd_euclidean(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, 
                 const fq_ctx_t ctx);

long _fq_poly_gcd_euclidean(fq_struct* G,const fq_struct* A, long lenA, 
                            const fq_struct* B, long lenB, const fq_ctx_t ctx);

static __inline__
long _fq_poly_gcd(fq_struct* G, const fq_struct* A, long lenA, 
                  const fq_struct* B, long lenB, const fq_ctx_t ctx)
{
    return _fq_poly_gcd_euclidean(G, A, lenA, B, lenB, ctx);
}

static __inline__
void fq_poly_gcd(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, 
                 const fq_ctx_t ctx)
{
    fq_poly_gcd_euclidean(rop,op1,op2,ctx);
}


/*  Euclidean division  ******************************************************/

#define FQ_POLY_DIVREM_DIVCONQUER_CUTOFF  16

ulong fq_poly_remove(fq_poly_t f, const fq_poly_t g, const fq_ctx_t ctx);

void _fq_poly_div_basecase(fq_struct *Q, fq_struct *R, 
                           const fq_struct *A, slong lenA,
                           const fq_struct *B, slong lenB, 
                           const fq_t invB, const fq_ctx_t ctx);

void fq_poly_div_basecase(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B,
                          const fq_ctx_t ctx);

void _fq_poly_divrem_basecase(fq_struct *Q, fq_struct *R, 
    const fq_struct *A, long lenA, const fq_struct *B, long lenB, 
    const fq_t invB, const fq_ctx_t ctx);

void fq_poly_divrem_basecase(fq_poly_t Q, fq_poly_t R, 
                             const fq_poly_t A, const fq_poly_t B, 
                             const fq_ctx_t ctx);

void
_fq_poly_divrem_divconquer_recursive(fq_struct * Q, fq_struct * BQ, fq_struct * W, 
                                     const fq_struct * A,
                                     const fq_struct * B, slong lenB, 
                                     const fq_t invB,
                                     const fq_ctx_t ctx);



void _fq_poly_divrem_divconquer(fq_struct *Q, fq_struct *R, 
                                const fq_struct *A, slong lenA,
                                const fq_struct *B, slong lenB, 
                                const fq_t invB,
                                const fq_ctx_t ctx);

void
fq_poly_divrem_divconquer(fq_poly_t Q, fq_poly_t R,
                          const fq_poly_t A, const fq_poly_t B,
                          const fq_ctx_t ctx);

static __inline__ 
void _fq_poly_divrem(fq_struct *Q, fq_struct *R, 
    const fq_struct *A, long lenA, const fq_struct *B, long lenB, 
    const fq_t invB, const fq_ctx_t ctx)
{
    _fq_poly_divrem_basecase(Q, R, A, lenA, B, lenB, invB, ctx);
}

static __inline__ 
void fq_poly_divrem(fq_poly_t Q, fq_poly_t R, 
                    const fq_poly_t A, const fq_poly_t B, 
                    const fq_ctx_t ctx)
{
    fq_poly_divrem_basecase(Q, R, A, B, ctx);
}

static __inline__ 
void _fq_poly_rem(fq_struct *R, const fq_struct *A, long lenA,
                  const fq_struct *B, long lenB, const fq_t invB,
                  const fq_ctx_t ctx)
{
    fq_struct *Q = _fq_vec_init(lenA + lenB, ctx); /*TODO: smaller bound */
    _fq_poly_divrem(Q, R, A, lenA, B, lenB, invB, ctx);
    _fq_vec_clear(Q, lenA + lenB, ctx);
}



static __inline__ 
void fq_poly_rem(fq_poly_t R, 
                    const fq_poly_t A, const fq_poly_t B, 
                    const fq_ctx_t ctx)
{
    fq_poly_t Q;
    fq_poly_init2(Q,A->length+B->length, ctx);/*TDOO: smaller bound*/
    fq_poly_divrem_basecase(Q, R, A, B, ctx);
    fq_poly_clear(Q, ctx);
}

void 
_fq_poly_inv_series_newton(fq_struct * Qinv, const fq_struct * Q, slong n, 
                           const fq_t cinv, const fq_ctx_t ctx);

void fq_poly_inv_series_newton(fq_poly_t Qinv, const fq_poly_t Q, slong n,
                               const fq_ctx_t ctx);

void _fq_poly_div_newton_preinv (fq_struct *Q, const fq_struct *A, slong lenA,
                                 const fq_struct* B, slong lenB,
                                 const fq_struct* Binv, slong lenBinv,
                                 const fq_ctx_t ctx);

void fq_poly_div_newton_preinv (fq_poly_t Q, const fq_poly_t A,
                                const fq_poly_t B, const fq_poly_t Binv,
                                const fq_ctx_t ctx);

void
_fq_poly_divrem_newton_preinv (fq_struct* Q, fq_struct* R,
                               const fq_struct* A, slong lenA,
                               const fq_struct* B, slong lenB,
                               const fq_struct* Binv, slong lenBinv, 
                               const fq_ctx_t ctx);

void
fq_poly_divrem_newton_preinv(fq_poly_t Q, fq_poly_t R,
                             const fq_poly_t A, const fq_poly_t B,
                             const fq_poly_t Binv, const fq_ctx_t ctx);

/*  Divisibility testing  ***************************************************/

int _fq_poly_divides(fq_struct *Q, 
                     const fq_struct *A, long lenA, 
                     const fq_struct *B, long lenB, const fq_t invB, 
                     const fq_ctx_t ctx);

int fq_poly_divides(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, 
                                 const fq_ctx_t ctx);

/*  Derivative  **************************************************************/

void _fq_poly_derivative(fq_struct *rop, const fq_struct *op, long len, 
                                         const fq_ctx_t ctx);

void fq_poly_derivative(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx);

/*  Evaluation  **************************************************************/

void _fq_poly_evaluate_fq(fq_t rop, const fq_struct *op, long len, 
                                    const fq_t a, const fq_ctx_t ctx);

void fq_poly_evaluate_fq(fq_t res, const fq_poly_t f, const fq_t a, 
                         const fq_ctx_t ctx);

/*  Composition  *************************************************************/

void _fq_poly_compose_divconquer(fq_struct *rop, 
                                 const fq_struct *op1, long len1, 
                                 const fq_struct *op2, long len2, 
                                 const fq_ctx_t ctx);

void fq_poly_compose_divconquer(fq_poly_t rop, 
                                const fq_poly_t op1, const fq_poly_t op2, 
                                const fq_ctx_t ctx);

void _fq_poly_compose_horner(fq_struct *rop, const fq_struct *op1, long len1, 
                                             const fq_struct *op2, long len2, 
                                             const fq_ctx_t ctx);

void fq_poly_compose_horner(fq_poly_t rop, 
                            const fq_poly_t op1, const fq_poly_t op2, 
                            const fq_ctx_t ctx);

void _fq_poly_compose(fq_struct *rop, const fq_struct *op1, long len1, 
                                      const fq_struct *op2, long len2, 
                                      const fq_ctx_t ctx);

void fq_poly_compose(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, 
                     const fq_ctx_t ctx);

void
_fq_poly_compose_mod(fq_struct * res, 
                     const fq_struct * f, slong lenf, 
                     const fq_struct * g,
                     const fq_struct * h, slong lenh, 
                     const fq_ctx_t ctx);

void
fq_poly_compose_mod(fq_poly_t res, const fq_poly_t poly1,
                    const fq_poly_t poly2, const fq_poly_t poly3,
                    const fq_ctx_t ctx);

void
_fq_poly_compose_mod_horner(fq_struct * res,
                            const fq_struct * f, slong lenf, 
                            const fq_struct * g,
                            const fq_struct * h, slong lenh, 
                            const fq_ctx_t ctx);

void
fq_poly_compose_mod_horner(fq_poly_t res,
                           const fq_poly_t poly1,
                           const fq_poly_t poly2,
                           const fq_poly_t poly3,
                           const fq_ctx_t ctx);

void
fq_poly_compose_mod_brent_kung(fq_poly_t res, const fq_poly_t poly1,
                               const fq_poly_t poly2, const fq_poly_t poly3,
                               const fq_ctx_t ctx);

void
_fq_poly_compose_mod_brent_kung(fq_struct * res,
                                const fq_struct * poly1, slong len1,
                                const fq_struct * poly2,
                                const fq_struct * poly3, slong len3,
                                const fq_ctx_t ctx);

void
_fq_poly_compose_mod_brent_kung_preinv(fq_struct * res,
                                       const fq_struct * poly1, slong len1,
                                       const fq_struct * poly2,
                                       const fq_struct * poly3, slong len3,
                                       const fq_struct * poly3inv, slong len3inv,
                                       const fq_ctx_t ctx);
void
fq_poly_compose_mod_brent_kung_preinv(fq_poly_t res, const fq_poly_t poly1,
                                      const fq_poly_t poly2, const fq_poly_t poly3,
                                      const fq_poly_t poly3inv, const fq_ctx_t ctx);


/*  Input and output  ********************************************************/

int _fq_poly_fprint_pretty(FILE *file, const fq_struct *poly, long len, 
                            const char *x, const fq_ctx_t ctx);

int fq_poly_fprint_pretty(FILE * file, const fq_poly_t poly, const char *x, 
                          const fq_ctx_t ctx);

int _fq_poly_fprint(FILE * file, const fq_struct *poly, long len, 
                    const fq_ctx_t ctx);

int fq_poly_fprint(FILE * file, const fq_poly_t poly, 
                   const fq_ctx_t ctx);

static __inline__ 
int _fq_poly_print(const fq_struct *poly, slong len, const fq_ctx_t ctx)
{
    return _fq_poly_fprint(stdout, poly, len, ctx);
}

static __inline__
int fq_poly_print(const fq_poly_t poly, const fq_ctx_t ctx)
{
    return fq_poly_fprint(stdout, poly, ctx);
}


static __inline__ 
int _fq_poly_print_pretty(const fq_struct *poly, long len, 
                          const char *x, const fq_ctx_t ctx)
{
    return _fq_poly_fprint_pretty(stdout, poly, len, x, ctx);
}

static __inline__ 
int fq_poly_print_pretty(const fq_poly_t poly, const char *x, 
                         const fq_ctx_t ctx)
{
    return fq_poly_fprint_pretty(stdout, poly, x, ctx);
}

#endif

