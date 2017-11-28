/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MPOLY_H
#define NMOD_MPOLY_H

#ifdef NMOD_MPOLY_INLINES_C
#define NMOD_MPOLY_INLINE FLINT_DLL
#else
#define NMOD_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif


/*
    all of the data we need to do arithmetic mod n
    and low level arithmetic modulo polynomials mod n
*/
typedef struct
{
    nmod_t mod;
    slong deg;    /* degree of finite field over its prime field */
    mp_limb_t * minpoly; /* minpoly[0] is constant term, minpoly[deg] should be 1 */
} nmodf_ctx_struct;

typedef nmodf_ctx_struct nmodf_ctx_t[1];

NMOD_MPOLY_INLINE
void nmodf_ctx_init(nmodf_ctx_t ctx, ulong modulus)
{
    ctx->mod.n = modulus;
    ctx->mod.ninv = n_preinvert_limb(modulus);
    count_leading_zeros(ctx->mod.norm, modulus);

    ctx->deg = 1;

    ctx->minpoly = (ulong *) flint_malloc(2*sizeof(ulong));
    ctx->minpoly[0] = 0;
    ctx->minpoly[1] = 1;
}

NMOD_MPOLY_INLINE
void nmodf_ctx_clear(nmodf_ctx_t ctx)
{
    if (ctx->minpoly != NULL)
        flint_free(ctx->minpoly);
}


/*  Type definitions *********************************************************/

typedef struct
{
    nmodf_ctx_t ffinfo;
    slong n;        /* number of elements in exponent vector (including deg) */
    ordering_t ord; /* polynomial ordering */
} nmod_mpoly_ctx_struct;

typedef nmod_mpoly_ctx_struct nmod_mpoly_ctx_t[1];

typedef struct
{
   mp_limb_t * coeffs;
   ulong * exps;  
   slong alloc;
   slong length;
   slong bits;     /* number of bits per exponent */
} nmod_mpoly_struct;

typedef nmod_mpoly_struct nmod_mpoly_t[1];

/* sparse univariates with multivariate coefficients */
typedef struct
{
   nmod_mpoly_struct * coeffs; /* multivariate coefficients */
   ulong * exps;
   slong alloc;
   slong length;
   slong var; /* univariate variable number */
} nmod_mpoly_univar_struct;

typedef nmod_mpoly_univar_struct nmod_mpoly_univar_t[1];

/* geobuckets ****************************************************************/
typedef struct nmod_mpoly_geobucket
{
    nmod_mpoly_struct polys[FLINT_BITS/2];
    slong length;
} nmod_mpoly_geobucket_struct;

typedef nmod_mpoly_geobucket_struct nmod_mpoly_geobucket_t[1];

FLINT_DLL void nmod_mpoly_geobucket_init(nmod_mpoly_geobucket_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_clear(nmod_mpoly_geobucket_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_empty(nmod_mpoly_t p,
                         nmod_mpoly_geobucket_t B, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_print(nmod_mpoly_geobucket_t B,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_fit_length(nmod_mpoly_geobucket_t B,
                                          slong i, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_set(nmod_mpoly_geobucket_t B,
                                   nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_add(nmod_mpoly_geobucket_t B,
                                   nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_sub(nmod_mpoly_geobucket_t B,
                                   nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_set_ui(nmod_mpoly_geobucket_t B,
                                          ulong c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_gen(nmod_mpoly_geobucket_t B, slong var,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_add_inplace(nmod_mpoly_geobucket_t B1,
                        nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_sub_inplace(nmod_mpoly_geobucket_t B1,
                        nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_neg_inplace(nmod_mpoly_geobucket_t B1,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_mul_inplace(nmod_mpoly_geobucket_t B1,
                        nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_pow_inplace(nmod_mpoly_geobucket_t B1,
                                          slong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_geobucket_divides_inplace(nmod_mpoly_geobucket_t B1,
                        nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx);

/* Context object ************************************************************/

FLINT_DLL void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx, 
                             slong nvars, const ordering_t ord, ulong modulus);

FLINT_DLL void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx);

/*  Memory management ********************************************************/

FLINT_DLL void nmod_mpoly_init(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_init2(nmod_mpoly_t poly, slong alloc, 
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_realloc(ulong ** poly, ulong ** exps,
                                            slong * alloc, slong len, slong N, const nmodf_ctx_t fctx);

FLINT_DLL void nmod_mpoly_realloc(nmod_mpoly_t poly, slong alloc, 
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_fit_length(ulong ** poly,
                             ulong ** exps, slong * alloc, slong len, slong N, const nmodf_ctx_t fctx);

FLINT_DLL void nmod_mpoly_fit_length(nmod_mpoly_t poly, slong len, 
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_clear(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void _nmod_mpoly_set_length(nmod_mpoly_t poly, slong newlen, 
                                                   const nmod_mpoly_ctx_t ctx)
{
    poly->length = newlen;
}

NMOD_MPOLY_INLINE
void nmod_mpoly_truncate(nmod_mpoly_t poly, slong newlen, 
                                                   const nmod_mpoly_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        poly->length = newlen;
    }  
}

/*
   if poly->bits < bits, set poly->bits = bits and reallocate poly->exps
*/
NMOD_MPOLY_INLINE
void nmod_mpoly_fit_bits(nmod_mpoly_t poly,
                                        slong bits, const nmod_mpoly_ctx_t ctx)
{
   slong N;
   ulong * t;

   FLINT_ASSERT(bits <= FLINT_BITS);

   if (poly->bits < bits)
   {
      if (poly->alloc != 0)
      {
         N = words_per_exp(ctx->n, bits);
         t = flint_malloc(N*poly->alloc*sizeof(ulong));
         mpoly_unpack_monomials(t, bits, poly->exps,
                                             poly->bits, poly->length, ctx->n);
         flint_free(poly->exps);
         poly->exps = t;
      }

      poly->bits = bits;
   }
}

/*  Basic manipulation *******************************************************/

FLINT_DLL void nmod_mpoly_degrees(slong * degs, const nmod_mpoly_t poly,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL slong nmod_mpoly_degree(const nmod_mpoly_t poly, slong var,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_max_degrees(ulong * max_degs, const ulong * exps,
                    slong len, slong bits, slong n, int deg, int rev, slong N);

FLINT_DLL void nmod_mpoly_max_degrees(ulong * max_degs,
                          const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_gen(ulong * poly, ulong * exps, slong i,
                               slong bits, slong n, int deg, int rev, slong N);

FLINT_DLL void nmod_mpoly_gen(nmod_mpoly_t poly, slong i,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_ui(nmod_mpoly_t poly,
                                          ulong c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_si(nmod_mpoly_t poly,
                                          slong c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_nmod(nmod_mpoly_t poly,
                                   const nmod_t c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_equal_ui(const nmod_mpoly_t poly,
                                          ulong c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_equal_si(const nmod_mpoly_t poly,
                                          slong c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_equal_nmod(const nmod_mpoly_t poly,
                                   const nmod_t c, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
void nmod_mpoly_swap(nmod_mpoly_t poly1, 
                                nmod_mpoly_t poly2, const nmod_mpoly_ctx_t ctx)
{
   nmod_mpoly_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}

NMOD_MPOLY_INLINE
void nmod_mpoly_zero(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
   _nmod_mpoly_set_length(poly, 0, ctx);
}

NMOD_MPOLY_INLINE
void nmod_mpoly_one(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_set_ui(poly, UWORD(1), ctx);
}

NMOD_MPOLY_INLINE
int nmod_mpoly_is_zero(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
   return poly->length == 0;
}

NMOD_MPOLY_INLINE
int nmod_mpoly_is_one(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
   return nmod_mpoly_equal_ui(poly, 1, ctx);
}

FLINT_DLL int nmod_mpoly_is_gen(const nmod_mpoly_t poly,
                                          slong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL ulong nmod_mpoly_get_coeff_ui(nmod_t x,
                 const nmod_mpoly_t poly, slong n, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_coeff_ui(nmod_mpoly_t poly, 
                          slong n, ulong x, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_get_monomial(ulong * exps, const nmod_mpoly_t poly, 
                                          slong n, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_monomial(nmod_mpoly_t poly, 
                      slong n, const ulong * exps, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_term_ui(nmod_mpoly_t poly,
                       ulong const * exp, ulong c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_get_term_ui(nmod_t c, const nmod_mpoly_t poly,
                                ulong const * exp, const nmod_mpoly_ctx_t ctx);

/* Set and negate ************************************************************/

FLINT_DLL void nmod_mpoly_set(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_neg(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                   const nmod_mpoly_ctx_t ctx);

/* Comparison ****************************************************************/

FLINT_DLL int nmod_mpoly_equal(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                   const nmod_mpoly_ctx_t ctx);

/* Basic arithmetic **********************************************************/

FLINT_DLL void nmod_mpoly_add_ui(nmod_mpoly_t poly1,
                const nmod_mpoly_t poly2, ulong c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_sub_ui(nmod_mpoly_t poly1,
                const nmod_mpoly_t poly2, ulong c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_add(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                         const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_sub(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                         const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx);

/* Scalar operations *********************************************************/

FLINT_DLL void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t poly1,
                const nmod_mpoly_t poly2, ulong c, const nmod_mpoly_ctx_t ctx);

/* Multiplication ************************************************************/

FLINT_DLL void nmod_mpoly_mul_johnson(nmod_mpoly_t poly1,
                 const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, 
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t poly1,
                 const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

FLINT_DLL void nmod_mpoly_pow(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                          slong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_pow_rmul(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                          slong k, const nmod_mpoly_ctx_t ctx);

/* Calculus ******************************************************************/

FLINT_DLL void nmod_mpoly_derivative(nmod_mpoly_t poly1,
              const nmod_mpoly_t poly2, slong var, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_integral(nmod_mpoly_t poly1,
              const nmod_mpoly_t poly2, slong var, const nmod_mpoly_ctx_t ctx);

/* Divisibility **************************************************************/

FLINT_DLL int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t poly1,
                  const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

/* Division ******************************************************************/

FLINT_DLL void nmod_mpoly_div_monagan_pearce(nmod_mpoly_t q,
                     const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_divrem_monagan_pearce(nmod_mpoly_t q, nmod_mpoly_t r,
                  const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

/* Evaluation ****************************************************************/

void nmod_mpoly_evaluate_all_tree_nmod(nmod_t ev, nmod_mpoly_t poly,
                                           ulong ** val, nmod_mpoly_ctx_t ctx);

void nmod_mpoly_evaluate_one_nmod(nmod_mpoly_t poly1, nmod_mpoly_t poly2,
                                  slong var, nmod_t val, nmod_mpoly_ctx_t ctx);

void nmod_mpoly_compose(nmod_mpoly_t res, nmod_mpoly_t poly1,
    nmod_mpoly_struct ** polys2, nmod_mpoly_ctx_t ctx1, nmod_mpoly_ctx_t ctx2);

/* Univariates ***************************************************************/

FLINT_DLL void nmod_mpoly_univar_init(nmod_mpoly_univar_t poly,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_univar_clear(nmod_mpoly_univar_t poly,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_univar_swap(nmod_mpoly_univar_t poly1,
                        nmod_mpoly_univar_t poly2, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_univar_fit_length(nmod_mpoly_univar_t poly,
                                     slong length, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_univar_print_pretty(const nmod_mpoly_univar_t poly,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_from_univar(nmod_mpoly_t poly1,
                  const nmod_mpoly_univar_t poly2, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_to_univar(nmod_mpoly_univar_t poly1,
              const nmod_mpoly_t poly2, slong var, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_univar_equal(nmod_mpoly_univar_t poly1,
                  const nmod_mpoly_univar_t poly2, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_univar_add(nmod_mpoly_univar_t poly1,
            const nmod_mpoly_univar_t poly2, const nmod_mpoly_univar_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_univar_mul(nmod_mpoly_univar_t poly1,
            const nmod_mpoly_univar_t poly2, const nmod_mpoly_univar_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_to_nmod_poly(nmod_poly_t poly1, slong * shift1,
              const nmod_mpoly_t poly2, slong var, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_from_nmod_poly(nmod_mpoly_t poly1,
        const nmod_poly_t poly2, slong shift2, slong var,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_univar_derivative(nmod_mpoly_univar_t poly1,
                  const nmod_mpoly_univar_t poly2, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_univar_prem(nmod_mpoly_univar_t polyA,
            const nmod_mpoly_univar_t polyB, nmod_mpoly_univar_t polyC,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_univar_pgcd(nmod_mpoly_univar_t poly1,
            const nmod_mpoly_univar_t polyP, const nmod_mpoly_univar_t polyQ,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_univar_pgcd_ducos(nmod_mpoly_univar_t poly1,
            const nmod_mpoly_univar_t polyP, const nmod_mpoly_univar_t polyQ,
                                                   const nmod_mpoly_ctx_t ctx);

/* GCD ***********************************************************************/

FLINT_DLL void nmod_mpoly_term_content(nmod_mpoly_t poly1,
                         const nmod_mpoly_t poly2, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_gcd_prs(nmod_mpoly_t poly1, nmod_mpoly_t poly2,
                               nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_gcd_is_unit(const nmod_mpoly_t a, const nmod_mpoly_t b,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_resultant(nmod_mpoly_t poly1,
                const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                        slong var, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_discriminant(nmod_mpoly_t poly1,
              const nmod_mpoly_t poly2, slong var, const nmod_mpoly_ctx_t ctx);

/* Reduction *****************************************************************/

FLINT_DLL slong
_nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** polyq, 
       ulong ** polyr, ulong ** expr, slong * allocr, const ulong * poly2,
          const ulong * exp2, slong len2, nmod_mpoly_struct * const * poly3,
                        ulong * const * exp3, slong len, slong N, slong bits,
                       const nmod_mpoly_ctx_t ctx, ulong maskhi, ulong masklo);

FLINT_DLL void
nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** q, nmod_mpoly_t r,
    const nmod_mpoly_t poly2, nmod_mpoly_struct * const * poly3, slong len,
                                                   const nmod_mpoly_ctx_t ctx);

/* Input/output **************************************************************/

FLINT_DLL int nmod_mpoly_set_str_pretty(nmod_mpoly_t poly, const char * str,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

FLINT_DLL char * nmod_mpoly_get_str_pretty(const nmod_mpoly_t poly,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int _nmod_mpoly_fprint_pretty(FILE * file, const ulong * poly, 
                           const ulong * exps, slong len, const char ** x,
             slong bits, slong n, int deg, int rev, slong N, const nmodf_ctx_t fctx);

FLINT_DLL int nmod_mpoly_fprint_pretty(FILE * file, 
         const nmod_mpoly_t poly, const char ** x, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
int _nmod_mpoly_print_pretty(const ulong * poly, 
                const ulong * exps, slong len, const char ** x,
              slong bits, slong n, int deg, int rev, slong N, const nmodf_ctx_t fctx)
{
   return _nmod_mpoly_fprint_pretty(stdout, poly, exps, len,
                                    x, bits, n, deg, rev, N, fctx);
}

NMOD_MPOLY_INLINE
int nmod_mpoly_print_pretty(const nmod_mpoly_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
{
   return nmod_mpoly_fprint_pretty(stdout, poly, x, ctx);
}

/* Random generation *********************************************************/

void nmod_mpoly_randtest(nmod_mpoly_t poly, flint_rand_t state,
                    slong length, slong exp_bound, const nmod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/






/******************************************************************************

   Internal consistency checks

******************************************************************************/

/*
   test that the terms in poly are in the correct order
*/
NMOD_MPOLY_INLINE
void nmod_mpoly_test(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
   slong i, N;
   ulong maskhi, masklo;

   masks_from_bits_ord(maskhi, masklo, poly->bits, ctx->ord);
   N = words_per_exp(ctx->n, poly->bits);

   if (!mpoly_monomials_test(poly->exps, poly->length, N, maskhi, masklo))
      flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

    for (i = 0; i < poly->length; i++)
    {
        if (poly->coeffs[i] == 0)
            flint_throw(FLINT_ERROR, "Polynomial has a zero coefficient");
    }

}


/*
   test that r is a valid remainder upon division by g over Q
   this means that no term of r is divisible by lt(g)
*/
NMOD_MPOLY_INLINE
void nmod_mpoly_remainder_strongtest(const nmod_mpoly_t r, const nmod_mpoly_t g,
                                                    const nmod_mpoly_ctx_t ctx)
{
   slong i, N, bits;
   ulong mask = 0;
   ulong * rexp, * gexp;

   bits = FLINT_MAX(r->bits, g->bits);
   N = words_per_exp(ctx->n, bits);

   if (g->length == 0 )
      flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

   if (r->length == 0 )
      return;

   rexp = (ulong *) flint_malloc(N*r->length*sizeof(ulong));
   gexp = (ulong *) flint_malloc(N*1        *sizeof(ulong));
   mpoly_unpack_monomials(rexp, bits, r->exps, r->bits, r->length, ctx->n);
   mpoly_unpack_monomials(gexp, bits, g->exps, g->bits, 1,         ctx->n);

   /* mask with high bit set in each field of exponent vector */
   for (i = 0; i < FLINT_BITS/bits; i++)
      mask = (mask << bits) + (UWORD(1) << (bits - 1));

   for (i = 0; i < r->length; i++)
      if (mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask))
      {
         flint_printf("nmod_mpoly_remainder_strongtest FAILED i = %wd\n", i);
         flint_printf("rem ");nmod_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
         flint_printf("den ");nmod_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
         flint_abort();
      }

   flint_free(rexp);
   flint_free(gexp);
}



#ifdef __cplusplus
}
#endif

#endif
