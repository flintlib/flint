/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_MPOLY_H
#define FMPQ_MPOLY_H

#ifdef FMPQ_MPOLY_INLINES_C
#define FMPQ_MPOLY_INLINE FLINT_DLL
#else
#define FMPQ_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "fmpz_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Context object ************************************************************/

typedef struct
{
    fmpz_mpoly_ctx_t zctx;
} fmpq_mpoly_ctx_struct;

typedef fmpq_mpoly_ctx_struct fmpq_mpoly_ctx_t[1];

FMPQ_MPOLY_INLINE
void fmpq_mpoly_ctx_init(fmpq_mpoly_ctx_t ctx, 
                                            slong nvars, const ordering_t ord)
{
    fmpz_mpoly_ctx_init(ctx->zctx, nvars, ord);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_ctx_init_rand(fmpq_mpoly_ctx_t ctx,
                                           flint_rand_t state, slong max_nvars)
{
    fmpz_mpoly_ctx_init_rand(ctx->zctx, state, max_nvars);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_ctx_clear(fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_ctx_clear(ctx->zctx);
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_ctx_nvars(fmpq_mpoly_ctx_t ctx)
{
    return ctx->zctx->minfo->nvars;
}

FMPQ_MPOLY_INLINE
ordering_t fmpq_mpoly_ctx_ord(fmpq_mpoly_ctx_t ctx)
{
    return ctx->zctx->minfo->ord;
}

/* Polynomials over Q ********************************************************/

/*
    A polynomial f is represented as
        content * zpoly,
    where zpoly should have positive leading coefficient and trivial content.
    If f is zero, then the representation should have
        content = 0 and zpoly = 0
*/

typedef struct
{                       /* non zero case:                   |  zero case: */
    fmpq_t content;     /* positive or negative content     |  (or zero)  */
    fmpz_mpoly_t zpoly; /* contentless poly, lc is positive |  (or zero)  */
} fmpq_mpoly_struct;

typedef fmpq_mpoly_struct fmpq_mpoly_t[1];


FMPQ_MPOLY_INLINE
void fmpq_mpoly_denominator(fmpz_t d, const fmpq_mpoly_t poly,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_set(d, fmpq_denref(poly->content));
}

/* geobuckets ****************************************************************/
typedef struct fmpq_mpoly_geobucket
{
    fmpq_mpoly_struct polys[FLINT_BITS/2];
    slong length;
} fmpq_mpoly_geobucket_struct;

typedef fmpq_mpoly_geobucket_struct fmpq_mpoly_geobucket_t[1];

FLINT_DLL void fmpq_mpoly_geobucket_init(fmpq_mpoly_geobucket_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_clear(fmpq_mpoly_geobucket_t B,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_empty(fmpq_mpoly_t p,
                         fmpq_mpoly_geobucket_t B, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_print(fmpq_mpoly_geobucket_t B,
                                  const char ** x, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_fit_length(fmpq_mpoly_geobucket_t B,
                                          slong i, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_geobucket_fix(fmpq_mpoly_geobucket_t B, slong i,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_set(fmpq_mpoly_geobucket_t B,
                                   fmpq_mpoly_t p, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_add(fmpq_mpoly_geobucket_t B,
                                   fmpq_mpoly_t p, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_sub(fmpq_mpoly_geobucket_t B,
                                   fmpq_mpoly_t p, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_set_fmpz(fmpq_mpoly_geobucket_t B,
                                         fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_set_fmpq(fmpq_mpoly_geobucket_t B,
                                         fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_gen(fmpq_mpoly_geobucket_t B, slong var,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_add_inplace(fmpq_mpoly_geobucket_t B1,
                        fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_sub_inplace(fmpq_mpoly_geobucket_t B1,
                        fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_neg_inplace(fmpq_mpoly_geobucket_t B1,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_mul_inplace(fmpq_mpoly_geobucket_t B1,
                        fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_pow_ui_inplace(fmpq_mpoly_geobucket_t B1,
                                          ulong k, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_geobucket_pow_fmpz_inplace(fmpq_mpoly_geobucket_t B1,
                                   const fmpz_t k, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_geobucket_divides_inplace(fmpq_mpoly_geobucket_t B1,
                        fmpq_mpoly_geobucket_t B2, const fmpq_mpoly_ctx_t ctx);

/*  Memory management ********************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_init(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_init(poly->content);
    fmpz_mpoly_init(poly->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_init2(fmpq_mpoly_t poly, slong alloc, 
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpq_init(poly->content);
    fmpz_mpoly_init2(poly->zpoly, alloc, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_realloc(fmpq_mpoly_t poly, slong alloc, 
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_realloc(poly->zpoly, alloc, ctx->zctx);
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_length(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_length(poly->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_fit_length(fmpq_mpoly_t poly, slong len, 
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_fit_length(poly->zpoly, len, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_clear(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_clear(poly->zpoly, ctx->zctx);
    fmpq_clear(poly->content);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_truncate(fmpq_mpoly_t poly, slong newlen, 
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_truncate(poly->zpoly, newlen, ctx->zctx);
}

/*
   if poly->bits < bits, set poly->bits = bits and reallocate poly->exps
*/
FMPQ_MPOLY_INLINE
void fmpq_mpoly_fit_bits(fmpq_mpoly_t poly,
                                        slong bits, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_fit_bits(poly->zpoly, bits, ctx->zctx);
}

/*  Basic manipulation *******************************************************/

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_fmpq(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
   return fmpz_mpoly_is_fmpz(poly->zpoly, ctx->zctx);
}

FLINT_DLL void fmpq_mpoly_get_fmpq(fmpq_t x, const fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_fmpq(fmpq_mpoly_t poly,
                                   const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_fmpz(fmpq_mpoly_t poly,
                                   const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_ui(fmpq_mpoly_t poly,
                                          ulong c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_si(fmpq_mpoly_t poly,
                                          slong c, const fmpq_mpoly_ctx_t ctx);


FLINT_DLL void fmpq_mpoly_canonicalise(fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
int fmpq_mpoly_degrees_fit_si(const fmpq_mpoly_t poly,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    return poly->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_degrees_fit_si(poly->zpoly->exps,
                                       poly->zpoly->length, poly->zpoly->bits,
                                                             ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_degrees_si(slong * degs, const fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, poly->zpoly->exps, poly->zpoly->length,
                                          poly->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_degrees_fmpz(fmpz ** degs, const fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, poly->zpoly->exps, poly->zpoly->length,
                                          poly->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_degree_si(const fmpq_mpoly_t poly, slong var,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_degree_si(poly->zpoly->exps, poly->zpoly->length,
                                     poly->zpoly->bits, var, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_degree_fmpz(fmpz_t deg, const fmpq_mpoly_t poly, slong var,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    mpoly_degree_fmpz(deg, poly->zpoly->exps, poly->zpoly->length,
                                     poly->zpoly->bits, var, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_totaldegree_fits_si(const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_totaldegree_fits_si(A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
slong fmpq_mpoly_totaldegree_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    return mpoly_totaldegree_si(A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_totaldegree_fmpz(fmpz_t td, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    mpoly_totaldegree_fmpz(td, A->zpoly->exps, A->zpoly->length,
                                             A->zpoly->bits, ctx->zctx->minfo);
}


FMPQ_MPOLY_INLINE
void fmpq_mpoly_gen(fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_one(poly->content);
    fmpz_mpoly_gen(poly->zpoly, i, ctx->zctx);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_gen(const fmpq_mpoly_t poly,
                                          slong k, const fmpq_mpoly_ctx_t ctx)
{
    return fmpq_is_one(poly->content)
          && fmpz_mpoly_is_gen(poly->zpoly, k, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_swap(fmpq_mpoly_t poly1, 
                                fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_mpoly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_zero(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_zero(poly->content);
    fmpz_mpoly_zero(poly->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_one(fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_one(poly->content);
    fmpz_mpoly_one(poly->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_zero(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_is_zero(poly->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_is_one(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
   return fmpq_is_one(poly->content)
          && fmpz_mpoly_is_one(poly->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_sort_terms(fmpq_mpoly_t poly1, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_sort_terms(poly1->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_combine_like_terms(fmpq_mpoly_t poly1, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_combine_like_terms(poly1->zpoly, ctx->zctx);
    fmpq_mpoly_canonicalise(poly1, ctx);
}

/* coefficients of monomials *************************************************/

/* get/set a coeff of poly1 corresponding to the monomial poly2 */
/* these two functions throw if poly2 is not a monomial */

FLINT_DLL void fmpq_mpoly_get_coeff_fmpq_monomial(fmpq_t c,
                        const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_coeff_fmpq_monomial(fmpq_mpoly_t poly1,
                                    const fmpq_t c, const fmpq_mpoly_t poly2,
                                                   const fmpq_mpoly_ctx_t ctx);

/* get/set a coeff of poly1 corresponding to an exponent vector exp */

FLINT_DLL void _fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t poly,
                 const fmpq_t c, const fmpz * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t poly,
                const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_coeff_fmpq_ui(fmpq_mpoly_t poly,
                const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t poly,
                                 const fmpz * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t poly,
                               fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_get_coeff_fmpq_ui(fmpq_t c, const fmpq_mpoly_t poly,
                                const ulong * exp, const fmpq_mpoly_ctx_t ctx);

/* pushing *******************************************************************/

FLINT_DLL void _fmpq_mpoly_emplacebackterm_fmpq_fmpz(fmpq_mpoly_t poly,
                     fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pushterm_fmpq_fmpz(fmpq_mpoly_t poly,
               const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pushterm_fmpz_fmpz(fmpq_mpoly_t poly,
               const fmpz_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pushterm_ui_fmpz(fmpq_mpoly_t poly,
                      ulong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pushterm_si_fmpz(fmpq_mpoly_t poly,
                      slong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void _fmpq_mpoly_emplacebackterm_fmpq_ui(fmpq_mpoly_t poly,
                       fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pushterm_fmpq_ui(fmpq_mpoly_t poly,
                const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pushterm_fmpz_ui(fmpq_mpoly_t poly,
                const fmpz_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pushterm_ui_ui(fmpq_mpoly_t poly,
                       ulong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pushterm_si_ui(fmpq_mpoly_t poly,
                       slong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx);

/* getters and setters for nth term ******************************************/

/* get/set the coefficient of the nth term into/from c */
/* these two functions throw if n is out of range */

FLINT_DLL void fmpq_mpoly_get_termcoeff_fmpq(fmpq_t c, const fmpq_mpoly_t poly,
                                          slong n, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_set_termcoeff_fmpq(fmpq_mpoly_t poly,
                          slong n, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

/* does the exponent vector of the nth term fit? */
FMPQ_MPOLY_INLINE
int fmpq_mpoly_termexp_fits_si(const fmpq_mpoly_t poly,
                                           slong n, const fmpq_mpoly_ctx_t ctx)
{
    return poly->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_termexp_fits_si(poly->zpoly->exps,
                                       poly->zpoly->bits, n, ctx->zctx->minfo);
}

FMPQ_MPOLY_INLINE
int fmpq_mpoly_termexp_fits_ui(const fmpq_mpoly_t poly,
                                           slong n, const fmpq_mpoly_ctx_t ctx)
{
    return poly->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_termexp_fits_ui(poly->zpoly->exps,
                                       poly->zpoly->bits, n, ctx->zctx->minfo);
}

/* get/set the exponent vector of the nth term into/from exps */
FMPQ_MPOLY_INLINE
void fmpq_mpoly_get_termexp_ui(ulong * exps, const fmpq_mpoly_t poly,
                                          slong n, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_get_termexp_ui(exps, poly->zpoly, n, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_get_termexp_fmpz(fmpz ** exps, const fmpq_mpoly_t poly,
                                          slong n, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_get_termexp_fmpz(exps, poly->zpoly, n, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_set_termexp_ui(fmpq_mpoly_t poly,
                      slong n, const ulong * exps, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_termexp_ui(poly->zpoly, n, exps, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_set_termexp_fmpz(fmpq_mpoly_t poly,
                      slong n, fmpz * const * exps, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_termexp_fmpz(poly->zpoly, n, exps, ctx->zctx);
}


/* Set and negate ************************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_set(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpq_set(poly1->content, poly2->content);
    fmpz_mpoly_set(poly1->zpoly, poly2->zpoly, ctx->zctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_neg(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpq_neg(poly1->content, poly2->content);
    fmpz_mpoly_set(poly1->zpoly, poly2->zpoly, ctx->zctx);
}


/* Comparison ****************************************************************/

FMPQ_MPOLY_INLINE
int fmpq_mpoly_equal(const fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    return fmpq_equal(poly1->content, poly2->content)
        && fmpz_mpoly_equal(poly1->zpoly, poly2->zpoly, ctx->zctx);
}

FLINT_DLL int fmpq_mpoly_equal_fmpq(const fmpq_mpoly_t poly, const fmpq_t c,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_equal_fmpz(const fmpq_mpoly_t poly, const fmpz_t c,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_equal_ui(const fmpq_mpoly_t poly, ulong        c,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_equal_si(const fmpq_mpoly_t poly, slong        c,
                                                   const fmpq_mpoly_ctx_t ctx);

/* Basic arithmetic **********************************************************/

FLINT_DLL void fmpq_mpoly_add_fmpq(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_add_fmpz(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_add_ui(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, ulong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_add_si(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, slong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub_fmpq(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub_fmpz(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub_ui(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, ulong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub_si(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, slong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_add(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                         const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_sub(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                         const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx);

/* Scalar operations *********************************************************/

FLINT_DLL void fmpq_mpoly_scalar_mul_fmpq(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_mul_fmpz(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_mul_ui(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, ulong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_mul_si(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, slong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_div_fmpq(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_div_fmpz(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, const fmpz_t c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_div_ui(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, ulong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_scalar_div_si(fmpq_mpoly_t poly1,
         const fmpq_mpoly_t poly2, slong        c, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_make_monic_inplace(fmpq_mpoly_t poly1,
                                                   const fmpq_mpoly_ctx_t ctx);

/* Multiplication ************************************************************/

FLINT_DLL void fmpq_mpoly_mul(fmpq_mpoly_t poly1,
                 const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3, 
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_mul_threaded(fmpq_mpoly_t poly1,
                 const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

FLINT_DLL void fmpq_mpoly_pow_fmpz(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                 const fmpz_t pow, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_pow_si(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                                        slong pow, const fmpq_mpoly_ctx_t ctx);

/* Calculus ******************************************************************/

FLINT_DLL void fmpq_mpoly_derivative(fmpq_mpoly_t poly1,
              const fmpq_mpoly_t poly2, slong var, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_integral(fmpq_mpoly_t poly1,
              const fmpq_mpoly_t poly2, slong var, const fmpq_mpoly_ctx_t ctx);

/* Divisibility **************************************************************/

FLINT_DLL int fmpq_mpoly_divides(fmpq_mpoly_t poly1,
                  const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

/* Division ******************************************************************/

FLINT_DLL void fmpq_mpoly_div(fmpq_mpoly_t q,
                     const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_divrem(fmpq_mpoly_t q, fmpq_mpoly_t r,
                  const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                                   const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_divrem_ideal(fmpq_mpoly_struct ** q, fmpq_mpoly_t r,
    const fmpq_mpoly_t poly2, fmpq_mpoly_struct * const * poly3, slong len,
                                                   const fmpq_mpoly_ctx_t ctx);

/* Evaluation ****************************************************************/
/*
not implemented yet
void fmpq_mpoly_evaluate_all_fmpq(fmpq_t ev, fmpq_mpoly_t poly,
                                            fmpq ** val, fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_evaluate_one_fmpq(fmpq_mpoly_t poly1, fmpq_mpoly_t poly2,
                                  slong var, fmpq_t val, fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_compose(fmpq_mpoly_t res, fmpq_mpoly_t poly1,
    fmpq_mpoly_struct ** polys2, fmpq_mpoly_ctx_t ctx1, fmpq_mpoly_ctx_t ctx2);
*/

/* GCD ***********************************************************************/

FLINT_DLL int fmpq_mpoly_gcd(fmpq_mpoly_t poly1, const fmpq_mpoly_t poly2,
                         const fmpq_mpoly_t poly3, const fmpq_mpoly_ctx_t ctx);
/*
not implemented yet
FLINT_DLL void fmpq_mpoly_resultant(fmpq_mpoly_t poly1,
                const fmpq_mpoly_t poly2, const fmpq_mpoly_t poly3,
                                        slong var, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL void fmpq_mpoly_discriminant(fmpq_mpoly_t poly1,
              const fmpq_mpoly_t poly2, slong var, const fmpq_mpoly_ctx_t ctx);
*/

/* Input/output **************************************************************/

FLINT_DLL int fmpq_mpoly_set_str_pretty(fmpq_mpoly_t poly, const char * str,
                                  const char ** x, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL char * fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t poly,
                                  const char ** x, const fmpq_mpoly_ctx_t ctx);

FLINT_DLL int fmpq_mpoly_fprint_pretty(FILE * file, 
         const fmpq_mpoly_t poly, const char ** x, const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_INLINE
int fmpq_mpoly_print_pretty(const fmpq_mpoly_t poly,
                                   const char ** x, const fmpq_mpoly_ctx_t ctx)
{
    return fmpq_mpoly_fprint_pretty(stdout, poly, x, ctx);
}

/* Random generation *********************************************************/

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bounds(fmpq_mpoly_t poly, flint_rand_t state,
                   slong length, mp_bitcnt_t coeff_bits, ulong * exp_bounds,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bounds(poly->zpoly, state,
                                     length, coeff_bits, exp_bounds, ctx->zctx);
    do {
        fmpq_randtest(poly->content, state, coeff_bits + 1);
    } while (fmpq_is_zero(poly->content));
    fmpq_mpoly_canonicalise(poly, ctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bound(fmpq_mpoly_t poly, flint_rand_t state,
                   slong length, mp_bitcnt_t coeff_bits, ulong exp_bound,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bound(poly->zpoly, state,
                                     length, coeff_bits, exp_bound, ctx->zctx);
    do {
        fmpq_randtest(poly->content, state, coeff_bits + 1);
    } while (fmpq_is_zero(poly->content));
    fmpq_mpoly_canonicalise(poly, ctx);
}

FMPQ_MPOLY_INLINE
void fmpq_mpoly_randtest_bits(fmpq_mpoly_t poly, flint_rand_t state,
                   slong length, mp_bitcnt_t coeff_bits, mp_bitcnt_t exp_bits,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_randtest_bits(poly->zpoly, state,
                                      length, coeff_bits, exp_bits, ctx->zctx);
    do {
        fmpq_randtest(poly->content, state, coeff_bits + 1);
    } while (fmpq_is_zero(poly->content));
    fmpq_mpoly_canonicalise(poly, ctx);
}



/******************************************************************************

   Internal consistency checks

******************************************************************************/

FLINT_DLL void fmpq_mpoly_assert_canonical(const fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx);


/*
   test that r is a valid remainder upon division by g over Q
   this means that no term of r is divisible by lt(g)
*/
FMPQ_MPOLY_INLINE
void fmpq_mpoly_remainder_test(const fmpq_mpoly_t r, const fmpq_mpoly_t g,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_remainder_strongtest(r->zpoly, g->zpoly, ctx->zctx);
}



#ifdef __cplusplus
}
#endif

#endif
