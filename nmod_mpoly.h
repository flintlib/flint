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

#include "fq_nmod_poly.h"


#ifdef __cplusplus
 extern "C" {
#endif


/* all of the data we need to do arithmetic mod n */

typedef struct
{
    nmod_t mod;
    ulong * extras; /* more information is sure to needed later */
} nmodf_ctx_struct; /* keep this here to enfore proper memory management of ctx struct */

typedef nmodf_ctx_struct nmodf_ctx_t[1];

NMOD_MPOLY_INLINE
void nmodf_ctx_init(nmodf_ctx_t ctx, ulong modulus)
{
    ctx->mod.n = modulus;
    ctx->mod.ninv = n_preinvert_limb(modulus);
    count_leading_zeros(ctx->mod.norm, modulus);

    ctx->extras = (ulong *) flint_malloc(2*sizeof(ulong));
}

NMOD_MPOLY_INLINE
void nmodf_ctx_reset(nmodf_ctx_t ctx, ulong modulus)
{
    ctx->mod.n = modulus;
    ctx->mod.ninv = n_preinvert_limb(modulus);
    count_leading_zeros(ctx->mod.norm, modulus);
}

NMOD_MPOLY_INLINE
void nmodf_ctx_clear(nmodf_ctx_t ctx)
{
    flint_free(ctx->extras);
}

/*  Type definitions *********************************************************/

typedef struct
{
    nmodf_ctx_t ffinfo;
    mpoly_ctx_t minfo;
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


/*
    A dense mpoly is stored as a flat array of coeffcients.
    Suppose deg_bounds = {a, b, c}. The coefficient of the monomial with 
    exponents {i, j, k} is stored at the coefficient of index
        c + k*(b + j*(a + i*0))    
*/
typedef struct
{
    slong nvars;
    slong degb_alloc;
    slong * deg_bounds;
    slong coeff_alloc;
    mp_limb_t * coeffs;
} nmod_mpolyd_struct;

typedef nmod_mpolyd_struct nmod_mpolyd_t[1];

typedef struct
{
    slong nvars;
    slong * perm;
} nmod_mpolyd_ctx_struct;

typedef nmod_mpolyd_ctx_struct nmod_mpolyd_ctx_t[1];


typedef struct
{
    slong nvars;
    slong degb_alloc;
    slong * deg_bounds;
    slong coeff_alloc;
    fq_nmod_struct * coeffs;
} fq_nmod_mpolyd_struct;
typedef fq_nmod_mpolyd_struct fq_nmod_mpolyd_t[1];

typedef struct
{
    slong nvars;
    slong * perm;
    fq_nmod_ctx_t fqctx;
} fq_nmod_mpolyd_ctx_struct;
typedef fq_nmod_mpolyd_ctx_struct fq_nmod_mpolyd_ctx_t[1];




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

FLINT_DLL void nmod_mpoly_geobucket_pow_si_inplace(nmod_mpoly_geobucket_t B1,
                                          slong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_geobucket_pow_fmpz_inplace(nmod_mpoly_geobucket_t B1,
                                   const fmpz_t k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_geobucket_divides_inplace(nmod_mpoly_geobucket_t B1,
                        nmod_mpoly_geobucket_t B2, const nmod_mpoly_ctx_t ctx);

/* Context object ************************************************************/

FLINT_DLL void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx, 
                         slong nvars, const ordering_t ord, mp_limb_t modulus);

FLINT_DLL void nmod_mpoly_ctx_init_rand(nmod_mpoly_ctx_t ctx, flint_rand_t state,
                                           slong max_nvars, mp_limb_t modulus);

FLINT_DLL void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx);

/*  Memory management ********************************************************/

FLINT_DLL void nmod_mpoly_init(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_init2(nmod_mpoly_t poly, slong alloc, 
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_init3(nmod_mpoly_t poly, slong alloc, mp_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_realloc(ulong ** poly, ulong ** exps,
                                            slong * alloc, slong len, slong N);

FLINT_DLL void nmod_mpoly_realloc(nmod_mpoly_t poly, slong alloc, 
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_fit_length(mp_limb_t ** poly,
                             ulong ** exps, slong * alloc, slong len, slong N);

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

   if (poly->bits < bits)
   {
      if (poly->alloc != 0)
      {
         N = mpoly_words_per_exp(bits, ctx->minfo);
         t = flint_malloc(N*poly->alloc*sizeof(ulong));
         mpoly_repack_monomials(t, bits, poly->exps,
                                         poly->bits, poly->length, ctx->minfo);
         flint_free(poly->exps);
         poly->exps = t;
      }

      poly->bits = bits;
   }
}

/*  Basic manipulation *******************************************************/

FLINT_DLL void _nmod_mpoly_set_term_ui_fmpz(nmod_mpoly_t poly,
                        ulong c, const fmpz * exp, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_term_ui_fmpz(nmod_mpoly_t poly,
                      ulong c, fmpz * const * exp, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_set_term_ui_ui(nmod_mpoly_t poly,
                       ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx);

FLINT_DLL ulong _nmod_mpoly_get_term_ui_fmpz(const nmod_mpoly_t poly,
                                 const fmpz * exp, const nmod_mpoly_ctx_t ctx);

FLINT_DLL ulong nmod_mpoly_get_term_ui_fmpz(const nmod_mpoly_t poly,
                               fmpz * const * exp, const nmod_mpoly_ctx_t ctx);

FLINT_DLL ulong nmod_mpoly_get_term_ui_ui(const nmod_mpoly_t poly,
                                const ulong * exp, const nmod_mpoly_ctx_t ctx);



FLINT_DLL void nmod_mpoly_degrees_si(slong * degs, const nmod_mpoly_t poly,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL slong nmod_mpoly_degree_si(const nmod_mpoly_t poly, slong var,
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

FLINT_DLL int nmod_mpoly_is_nmod(const nmod_mpoly_t poly,
                                                   const nmod_mpoly_ctx_t ctx);

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

/* Set and negate ************************************************************/

FLINT_DLL void nmod_mpoly_set(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_neg(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                   const nmod_mpoly_ctx_t ctx);

/* Comparison ****************************************************************/

FLINT_DLL int nmod_mpoly_equal(const nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
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

FLINT_DLL void nmod_mpoly_make_monic(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                    const nmod_mpoly_ctx_t ctx);

/* Multiplication ************************************************************/

FLINT_DLL void nmod_mpoly_mul(nmod_mpoly_t poly1,
                 const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, 
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_mul_johnson(nmod_mpoly_t poly1,
                 const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, 
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t poly1,
                 const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_mul_dense(nmod_mpoly_t P,
                        const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_mul_array(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                         const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx);

slong _nmod_mpoly_mul_johnson(mp_limb_t ** coeff1, ulong ** exp1, slong * alloc,
                 const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                 const mp_limb_t * coeff3, const ulong * exp3, slong len3,
      mp_bitcnt_t bits, slong N, const ulong * cmpmask, const nmodf_ctx_t fctx);

/* Powering ******************************************************************/

FLINT_DLL void nmod_mpoly_pow_si(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                          slong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_pow_rmul(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                          slong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_pow_fmpz(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                 const fmpz_t pow, const nmod_mpoly_ctx_t ctx);

/* Calculus ******************************************************************/

FLINT_DLL void nmod_mpoly_derivative(nmod_mpoly_t poly1,
              const nmod_mpoly_t poly2, slong var, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_integral(nmod_mpoly_t poly1,
              const nmod_mpoly_t poly2, slong var, const nmod_mpoly_ctx_t ctx);

/* Divisibility **************************************************************/

FLINT_DLL int nmod_mpoly_divides(nmod_mpoly_t Q,
                        const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t poly1,
                  const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_divides_dense(nmod_mpoly_t Q,
                        const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

slong _nmod_mpoly_divides_monagan_pearce(
                     mp_limb_t ** coeff1,      ulong ** exp1, slong * alloc,
                const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                const mp_limb_t * coeff3, const ulong * exp3, slong len3,
     mp_bitcnt_t bits, slong N, const ulong * cmpmask, const nmodf_ctx_t fctx);

/* Division ******************************************************************/

FLINT_DLL void nmod_mpoly_div_monagan_pearce(nmod_mpoly_t q,
                     const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_divrem_monagan_pearce(nmod_mpoly_t q, nmod_mpoly_t r,
                  const nmod_mpoly_t poly2, const nmod_mpoly_t poly3,
                                                   const nmod_mpoly_ctx_t ctx);

/* Evaluation ****************************************************************/
/*
not implemented yet
void nmod_mpoly_evaluate_all_tree_nmod(nmod_t ev, nmod_mpoly_t poly,
                                           ulong ** val, nmod_mpoly_ctx_t ctx);

void nmod_mpoly_evaluate_one_nmod(nmod_mpoly_t poly1, nmod_mpoly_t poly2,
                                  slong var, nmod_t val, nmod_mpoly_ctx_t ctx);

void nmod_mpoly_compose(nmod_mpoly_t res, nmod_mpoly_t poly1,
    nmod_mpoly_struct ** polys2, nmod_mpoly_ctx_t ctx1, nmod_mpoly_ctx_t ctx2);
*/
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

/* terms *********************************************************************/

FLINT_DLL void _nmod_mpoly_radix_sort1(nmod_mpoly_t A, slong left, slong right,
                              mp_bitcnt_t pos, ulong cmpmask, ulong totalmask);

FLINT_DLL void _nmod_mpoly_radix_sort(nmod_mpoly_t A, slong left, slong right,
                                    mp_bitcnt_t pos, slong N, ulong * cmpmask);

FLINT_DLL void nmod_mpoly_sort_terms(nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_combine_terms(nmod_mpoly_t A,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_emplacebackterm_ui_ui(nmod_mpoly_t poly,
                   mp_limb_t c, const ulong * exp, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_pushterm_ui_ui(nmod_mpoly_t poly,
                       ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void _nmod_mpoly_emplacebackterm_ui_ffmpz(nmod_mpoly_t poly,
                    mp_limb_t c, const fmpz * exp, const nmod_mpoly_ctx_t ctx);

/* dense helpers *************************************************************/

FLINT_DLL void nmod_mpolyd_ctx_init(nmod_mpolyd_ctx_t dctx, slong nvars);

FLINT_DLL void nmod_mpolyd_ctx_clear(nmod_mpolyd_ctx_t dctx);

NMOD_MPOLY_INLINE void nmod_mpolyd_swap(nmod_mpolyd_t poly1,
                                nmod_mpolyd_t poly2)
{
   nmod_mpolyd_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}

FLINT_DLL int nmod_mpolyd_set_degbounds(nmod_mpolyd_t A, slong * bounds);

FLINT_DLL int nmod_mpolyd_set_degbounds_perm(nmod_mpolyd_t A,
                                 const nmod_mpolyd_ctx_t dctx, slong * bounds);

FLINT_DLL void nmod_mpoly_convert_to_nmod_mpolyd(
                                  nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_convert_to_nmod_mpolyd_degbound(
                                  nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_convert_from_nmod_mpolyd(
                                 nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx,
                          const nmod_mpolyd_t B, const nmod_mpolyd_ctx_t dctx);

FLINT_DLL void nmod_mpolyd_init(nmod_mpolyd_t poly, slong nvars);

FLINT_DLL void nmod_mpolyd_fit_length(nmod_mpolyd_t poly, slong len);

FLINT_DLL void nmod_mpolyd_zero(nmod_mpolyd_t poly);

FLINT_DLL void nmod_mpolyd_set_nvars(nmod_mpolyd_t poly, slong nvars);

FLINT_DLL void nmod_mpolyd_set_ui(nmod_mpolyd_t A, mp_limb_t v);

FLINT_DLL void nmod_mpolyd_set(nmod_mpolyd_t A, const nmod_mpolyd_t B);

FLINT_DLL void nmod_mpolyd_clear(nmod_mpolyd_t poly);

FLINT_DLL void nmod_mpolyd_print(nmod_mpolyd_t poly);

FLINT_DLL slong nmod_mpolyd_length(const nmod_mpolyd_t A);

FLINT_DLL slong nmod_mpolyd_last_degree(const nmod_mpolyd_t A,
                                                       const nmodf_ctx_t fctx);

FLINT_DLL void nmod_mpolyd_div_last_poly(nmod_mpolyd_t A, nmod_poly_t b,
                                                       const nmodf_ctx_t fctx);

FLINT_DLL void nmod_mpolyd_last_lc(nmod_poly_t lc, const nmod_mpolyd_t A,
                                                       const nmodf_ctx_t fctx);

FLINT_DLL slong nmod_mpolyd_leadmon(slong * exps, const nmod_mpolyd_t A);

FLINT_DLL void nmod_mpolyd_mul_last_poly(nmod_mpolyd_t M,
                       nmod_mpolyd_t A, nmod_poly_t b, const nmodf_ctx_t fctx);

FLINT_DLL void nmod_mpolyd_mul_scalar(nmod_mpolyd_t A, mp_limb_t b,
                                                       const nmodf_ctx_t fctx);


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

FLINT_DLL int nmod_mpolyd_ctx_set_for_gcd(nmod_mpolyd_ctx_t dctx,
                            const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyd_last_content(nmod_poly_t cont, const nmod_mpolyd_t A,
                                                       const nmodf_ctx_t fctx);

FLINT_DLL void nmod_mpolyd_gcd_brown_univar(nmod_mpolyd_t G,
                                 nmod_mpolyd_t Abar,       nmod_mpolyd_t Bbar,
                           const nmod_mpolyd_t A   , const nmod_mpolyd_t B,
                                                       const nmodf_ctx_t fctx);

FLINT_DLL int nmod_mpolyd_gcd_brown_smprime(nmod_mpolyd_t G,
                          nmod_mpolyd_t Abar, nmod_mpolyd_t Bbar,
                          nmod_mpolyd_t A, nmod_mpolyd_t B,
                                                       const nmodf_ctx_t fctx);

FLINT_DLL int nmod_mpolyd_gcd_brown_lgprime(nmod_mpolyd_t G,
                                 nmod_mpolyd_t Abar,  nmod_mpolyd_t Bbar,
                                 nmod_mpolyd_t A, nmod_mpolyd_t B,
                                                       const nmodf_ctx_t fctx);

FLINT_DLL int fq_nmod_mpolyd_gcd_brown_smprime(fq_nmod_mpolyd_t G,
                                fq_nmod_mpolyd_t Abar, fq_nmod_mpolyd_t Bbar,
                                fq_nmod_mpolyd_t A, fq_nmod_mpolyd_t B,
                                              const fq_nmod_mpolyd_ctx_t dctx);

FLINT_DLL int nmod_mpoly_gcd_brown(nmod_mpoly_t G,
                               const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A,
                             const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int _nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx,
                                         int keepbits, flint_rand_t randstate);

/*
    nmod_mpolyu_t
    sparse univariates with nmod_mpoly_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
   nmod_mpoly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   mp_bitcnt_t bits;    /* default bits to construct coeffs */
} nmod_mpolyu_struct;
typedef nmod_mpolyu_struct nmod_mpolyu_t[1];

FLINT_DLL void nmod_mpolyu_init(nmod_mpolyu_t A, mp_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_clear(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

FLINT_DLL void nmod_mpolyu_swap(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                  const nmod_mpoly_ctx_t uctx);

FLINT_DLL void nmod_mpolyu_zero(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

FLINT_DLL int nmod_mpolyu_is_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

FLINT_DLL void nmod_mpolyu_print_pretty(const nmod_mpolyu_t poly,
                                  const char ** x, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_fit_length(nmod_mpolyu_t A, slong length,
                                                  const nmod_mpoly_ctx_t uctx);

FLINT_DLL void nmod_mpolyu_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

FLINT_DLL nmod_mpoly_struct * _nmod_mpolyu_get_coeff(nmod_mpolyu_t A,
                                       ulong pow, const nmod_mpoly_ctx_t uctx);

FLINT_DLL void nmod_mpolyu_shift_right(nmod_mpolyu_t A, ulong s);

FLINT_DLL void nmod_mpolyu_shift_left(nmod_mpolyu_t A, ulong s);

FLINT_DLL void nmod_mpolyu_scalar_mul_nmod(nmod_mpolyu_t A, mp_limb_t c,
                                                  const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_set(nmod_mpolyu_t A, const nmod_mpolyu_t B,
                                                  const nmod_mpoly_ctx_t uctx);

FLINT_DLL void nmod_mpolyu_cvtto_poly(nmod_poly_t a, nmod_mpolyu_t A,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_cvtfrom_poly(nmod_mpolyu_t A, nmod_poly_t a,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_cvtfrom_poly_notmain(nmod_mpolyu_t A, nmod_poly_t a,
                                        slong var, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_to_mpolyu_perm(nmod_mpolyu_t A,
                                      const nmod_mpoly_t B, const slong * perm,
                      const nmod_mpoly_ctx_t uctx, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_from_mpolyu_perm(nmod_mpoly_t A,
                       const nmod_mpolyu_t B, int keepbits, const slong * perm,
                      const nmod_mpoly_ctx_t uctx, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpolyu_divides(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_divexact_mpoly(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                   nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx);


FLINT_DLL void nmod_mpolyu_mul_mpoly(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                   nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx);




/*
    nmod_mpolyn_t
    multivariates with nmod_poly_t coefficients
*/
typedef struct
{
   nmod_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} nmod_mpolyn_struct;
typedef nmod_mpolyn_struct nmod_mpolyn_t[1];

FLINT_DLL void nmod_mpolyn_init(nmod_mpolyn_t A, mp_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_clear(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_swap(nmod_mpolyn_t A, nmod_mpolyn_t B);

FLINT_DLL void nmod_mpolyn_zero(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpolyn_is_zero(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_print_pretty(const nmod_mpolyn_t A, const char ** x_in,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_fit_length(nmod_mpolyn_t A, slong length,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_set_length(nmod_mpolyn_t A, slong newlen,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_fit_bits(nmod_mpolyn_t A, slong bits,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_set(nmod_mpolyn_t A, const nmod_mpolyn_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_set_mpoly(nmod_mpolyn_t A, const nmod_mpoly_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_mul_poly(nmod_mpolyn_t A, const nmod_mpolyn_t B, 
                              const nmod_poly_t c, const nmod_mpoly_ctx_t ctx);

/*
    nmod_mpolyun_t
    sparse univariates with nmod_mpolyn_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
    nmod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    mp_bitcnt_t bits;   /* default bits to construct coeffs */
} nmod_mpolyun_struct;
typedef nmod_mpolyun_struct nmod_mpolyun_t[1];

FLINT_DLL void nmod_mpolyun_init(nmod_mpolyun_t A, mp_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_clear(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_swap(nmod_mpolyun_t A, nmod_mpolyun_t B);

FLINT_DLL void nmod_mpolyun_zero(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_print_pretty(const nmod_mpolyun_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_fit_length(nmod_mpolyun_t A, slong length,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_shift_right(nmod_mpolyun_t A, ulong s);

FLINT_DLL void nmod_mpolyun_shift_left(nmod_mpolyun_t A, ulong s);

FLINT_DLL slong nmod_mpolyun_lastdeg(nmod_mpolyun_t A,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_set(nmod_mpolyun_t A, const nmod_mpolyun_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_cvtto_mpolyun(nmod_mpolyun_t A, nmod_mpolyu_t B,
                                          slong k, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_cvtfrom_mpolyun(nmod_mpolyu_t A, nmod_mpolyun_t B,
                                        slong var, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_set_mpolyu(nmod_mpolyun_t A, const nmod_mpolyu_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_eval_last(nmod_mpolyu_t B, nmod_mpolyun_t A,
                                  mp_limb_t alpha, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_mul_poly(nmod_mpolyun_t A, const nmod_mpolyun_t B,
                              const nmod_poly_t c, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_content_last(nmod_poly_t a, nmod_mpolyun_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_divexact_last(nmod_mpolyun_t A, nmod_poly_t b,
                                                   const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE nmod_poly_struct *
nmod_mpolyn_leadcoeff_ref(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

NMOD_MPOLY_INLINE nmod_poly_struct *
nmod_mpolyun_leadcoeff_ref(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpolyn_leadcoeff_ref(A->coeffs + 0, ctx);
}


NMOD_MPOLY_INLINE mp_limb_t
nmod_mpoly_leadcoeff(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs[0];
}

NMOD_MPOLY_INLINE mp_limb_t
nmod_mpolyu_leadcoeff(nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return nmod_mpoly_leadcoeff(A->coeffs + 0, ctx);
}

FLINT_DLL void nmod_mpolyu_init(nmod_mpolyu_t A, mp_bitcnt_t bits,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_fit_length(nmod_mpolyu_t A, slong length,
                                                  const nmod_mpoly_ctx_t uctx);

FLINT_DLL void nmod_mpolyu_clear(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx);

FLINT_DLL void nmod_mpolyu_setform(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyu_scalar_mul_nmod(nmod_mpolyu_t A, mp_limb_t c,
                                                   const nmod_mpoly_ctx_t ctx);

typedef enum {
    nmod_gcds_success,
    nmod_gcds_form_main_degree_too_high,
    nmod_gcds_form_wrong,
    nmod_gcds_no_solution,
    nmod_gcds_scales_not_found,
    nmod_gcds_eval_point_not_found,
    nmod_gcds_eval_gcd_deg_too_high
} nmod_gcds_ret_t;

FLINT_DLL nmod_gcds_ret_t nmod_mpolyu_gcds_zippel(nmod_mpolyu_t G,
                           nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpolyu_t f,
                                        slong var, const nmod_mpoly_ctx_t ctx,
                                     flint_rand_t randstate, slong * degbound);

FLINT_DLL int nmod_mpolyu_gcdp_zippel(nmod_mpolyu_t G,
                             nmod_mpolyu_t A, nmod_mpolyu_t B, slong var,
                            const nmod_mpoly_ctx_t ctx, mpoly_zipinfo_t zinfo,
                                                       flint_rand_t randstate);

typedef struct
{
    mpoly_ctx_t minfo;
    fq_nmod_ctx_t fqctx;
} fq_nmod_mpoly_ctx_struct;
typedef fq_nmod_mpoly_ctx_struct fq_nmod_mpoly_ctx_t[1];


/*
    fq_nmod_mpoly_t
    sparse multivariates with fq_nmod coefficients
*/
typedef struct
{
    fq_nmod_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    slong bits;     /* number of bits per exponent */
} fq_nmod_mpoly_struct;
typedef fq_nmod_mpoly_struct fq_nmod_mpoly_t[1];


FLINT_DLL void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                                                        mp_limb_t p, slong deg);

FLINT_DLL void fq_nmod_mpoly_ctx_clear(fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_ctx_change_modulus(fq_nmod_mpoly_ctx_t ctx,
                                                                    slong deg);

FLINT_DLL int fq_nmod_mpoly_is_canonical(const fq_nmod_mpoly_t poly,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_init(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_init2(fq_nmod_mpoly_t A, slong alloc,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_clear(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_zero(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_is_zero(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_is_one(const fq_nmod_mpoly_t poly,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A,
                            const char ** x_in, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_swap(fq_nmod_mpoly_t poly1, fq_nmod_mpoly_t poly2);

FLINT_DLL void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_fit_length(fq_nmod_struct ** coeff,
                              ulong ** exps, slong * alloc, slong len, slong N,
                                                    const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpoly_fit_bits(fq_nmod_mpoly_t A, slong bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_length(fq_nmod_mpoly_t A, slong newlen,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_length(fq_nmod_mpoly_t poly, slong newlen,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set(fq_nmod_mpoly_t poly1,
                   const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong _fq_nmod_mpoly_divides_monagan_pearce(
                  fq_nmod_struct ** coeff1,      ulong ** exp1, slong * alloc,
             const fq_nmod_struct * coeff2, const ulong * exp2, slong len2,
             const fq_nmod_struct * coeff3, const ulong * exp3, slong len3,
  mp_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);

FLINT_DLL int fq_nmod_mpoly_divides_monagan_pearce(fq_nmod_mpoly_t poly1,
                  const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong _fq_nmod_mpoly_mul_johnson(
                    fq_nmod_struct ** coeff1, ulong ** exp1, slong * alloc,
             const fq_nmod_struct * coeff2, const ulong * exp2, slong len2,
             const fq_nmod_struct * coeff3, const ulong * exp3, slong len3,
  mp_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpoly_mul_johnson(fq_nmod_mpoly_t poly1,
                    const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_sub(fq_nmod_mpoly_t poly1,
                    const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_neg(fq_nmod_mpoly_t poly1,
                   const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_ctx_t ctx);

/*
    fq_nmod_mpolyu_t
    sparse univariates with fq_nmod_mpoly_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
   fq_nmod_mpoly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   mp_bitcnt_t bits;    /* default bits to construct coeffs */
} fq_nmod_mpolyu_struct;
typedef fq_nmod_mpolyu_struct fq_nmod_mpolyu_t[1];

FLINT_DLL int fq_nmod_mpolyu_is_canonical(const fq_nmod_mpolyu_t poly,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_init(fq_nmod_mpolyu_t A, mp_bitcnt_t bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_clear(fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t uctx);

FLINT_DLL void fq_nmod_mpolyu_swap(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B);

FLINT_DLL void fq_nmod_mpolyu_zero(fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t uctx);

FLINT_DLL int fq_nmod_mpolyu_is_one(fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t uctx);

FLINT_DLL void fq_nmod_mpolyu_print_pretty(const fq_nmod_mpolyu_t poly,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_fit_length(fq_nmod_mpolyu_t A, slong length,
                                               const fq_nmod_mpoly_ctx_t uctx);

FLINT_DLL void fq_nmod_mpolyu_shift_right(fq_nmod_mpolyu_t A, ulong s);

FLINT_DLL void fq_nmod_mpolyu_shift_left(fq_nmod_mpolyu_t A, ulong s);

FLINT_DLL void fq_nmod_mpolyu_scalar_mul_fq_nmod(fq_nmod_mpolyu_t A,
                                   fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_set(fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B,
                                               const fq_nmod_mpoly_ctx_t uctx);

FLINT_DLL void fq_nmod_mpoly_cvtfrom_poly_notmain(fq_nmod_mpoly_t A,
                   fq_nmod_poly_t a, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtfrom_poly_notmain(fq_nmod_mpolyu_t A,
                   fq_nmod_poly_t a, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtto_poly(fq_nmod_poly_t a, fq_nmod_mpolyu_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtfrom_poly(fq_nmod_mpolyu_t A, fq_nmod_poly_t a,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_setform(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyu_divides(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_divexact_mpoly(fq_nmod_mpolyu_t A,
                                      fq_nmod_mpolyu_t B, fq_nmod_mpoly_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_mul_mpoly(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                             fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_set_fq_nmod_mpolyu(
                             nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx,
                          fq_nmod_mpolyu_t Ap, const fq_nmod_mpoly_ctx_t ctxp);

FLINT_DLL int nmod_mpolyun_CRT_fq_nmod_mpolyu(slong * lastdeg,
                             nmod_mpolyun_t H, const nmod_mpoly_ctx_t ctx,
            nmod_poly_t m, fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctxp);

FLINT_DLL void nmod_mpolyun_redto_fq_nmod_mpolyu(fq_nmod_mpolyu_t A, nmod_mpolyun_t B,
                  const fq_nmod_mpoly_ctx_t ffctx, const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpolyun_addinterp(slong * lastdeg,
                        nmod_mpolyun_t F, nmod_mpolyun_t T,
                        nmod_mpolyu_t A, nmod_poly_t modulus,  mp_limb_t alpha,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL int nmod_mpolyun_addinterp_fq_nmod_mpolyu(slong * lastdeg,
                                nmod_mpolyun_t F, nmod_mpolyun_t T,
                             nmod_poly_t m, const nmod_mpoly_ctx_t ctx,
                          fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ffctx);
/*
    fq_nmod_mpolyn_t
    multivariates with fq_nmod_poly_t coefficients
*/
typedef struct
{
   fq_nmod_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} fq_nmod_mpolyn_struct;
typedef fq_nmod_mpolyn_struct fq_nmod_mpolyn_t[1];

FLINT_DLL void fq_nmod_mpolyn_init(fq_nmod_mpolyn_t A, mp_bitcnt_t bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_clear(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_swap(fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B);

FLINT_DLL void fq_nmod_mpolyn_zero(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyn_is_zero(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_print_pretty(const fq_nmod_mpolyn_t A,
                             const char ** x_in, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_fit_length(fq_nmod_mpolyn_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_set_length(fq_nmod_mpolyn_t A, slong newlen,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_fit_bits(fq_nmod_mpolyn_t A, slong bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_set(fq_nmod_mpolyn_t A, const fq_nmod_mpolyn_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_mul_poly(fq_nmod_mpolyn_t A,
                         const fq_nmod_mpolyn_t B, const fq_nmod_poly_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

/*
    fq_nmod_mpolyun_t
    sparse univariates with fq_nmod_mpolyn_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
    fq_nmod_mpolyn_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    mp_bitcnt_t bits;   /* default bits to construct coeffs */
} fq_nmod_mpolyun_struct;
typedef fq_nmod_mpolyun_struct fq_nmod_mpolyun_t[1];

FLINT_DLL void fq_nmod_mpolyun_init(fq_nmod_mpolyun_t A, mp_bitcnt_t bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_clear(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_swap(fq_nmod_mpolyun_t A, fq_nmod_mpolyun_t B);

FLINT_DLL void fq_nmod_mpolyun_zero(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_print_pretty(const fq_nmod_mpolyun_t poly,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_fit_length(fq_nmod_mpolyun_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_shift_right(fq_nmod_mpolyun_t A, ulong s);

FLINT_DLL void fq_nmod_mpolyun_shift_left(fq_nmod_mpolyun_t A, ulong s);

FLINT_DLL slong fq_nmod_mpolyun_lastdeg(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_set(fq_nmod_mpolyun_t A,
                     const fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtto_mpolyun(fq_nmod_mpolyun_t A,
                   fq_nmod_mpolyu_t B, slong k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtfrom_mpolyun(fq_nmod_mpolyu_t A,
                fq_nmod_mpolyun_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_set_mpolyu(fq_nmod_mpolyun_t A,
                      const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_mul_poly(fq_nmod_mpolyun_t A,
                         const fq_nmod_mpolyun_t B, const fq_nmod_poly_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_content_last(fq_nmod_poly_t a,
                           fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_divexact_last(fq_nmod_mpolyun_t A,
                              fq_nmod_poly_t b, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_eval_last(fq_nmod_mpolyu_t B,
                                     fq_nmod_mpolyun_t A, fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyun_addinterp(slong * lastdeg,
            fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_mpolyu_t A,
       fq_nmod_poly_t modulus, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL nmod_gcds_ret_t fq_nmod_mpolyu_gcds_zippel(fq_nmod_mpolyu_t G,
                 fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, fq_nmod_mpolyu_t f,
                                     slong var, const fq_nmod_mpoly_ctx_t ctx,
                                     flint_rand_t randstate, slong * degbound);

FLINT_DLL int fq_nmod_mpolyu_gcdp_zippel(fq_nmod_mpolyu_t G,
                         fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, slong var,
                         const fq_nmod_mpoly_ctx_t ctx, mpoly_zipinfo_t zinfo,
                                                       flint_rand_t randstate);

int fq_nmod_mpolyu_gcdp_zippel_univar(fq_nmod_mpolyu_t G,
        fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx);


NMOD_MPOLY_INLINE fq_nmod_poly_struct *
fq_nmod_mpolyn_leadcoeff_ref(const fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

NMOD_MPOLY_INLINE fq_nmod_poly_struct *
fq_nmod_mpolyun_leadcoeff_ref(const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return fq_nmod_mpolyn_leadcoeff_ref(A->coeffs + 0, ctx);
}


NMOD_MPOLY_INLINE fq_nmod_struct *
fq_nmod_mpoly_leadcoeff_ref(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

NMOD_MPOLY_INLINE fq_nmod_struct *
fq_nmod_mpolyu_leadcoeff_ref(const fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return fq_nmod_mpoly_leadcoeff_ref(A->coeffs + 0, ctx);
}


/* Reduction *****************************************************************/

FLINT_DLL slong
_nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** polyq, 
       ulong ** polyr, ulong ** expr, slong * allocr, const ulong * poly2,
          const ulong * exp2, slong len2, nmod_mpoly_struct * const * poly3,
                        ulong * const * exp3, slong len, slong N, slong bits,
                       const nmod_mpoly_ctx_t ctx, const ulong * cmpmask);

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
                   slong bits, const mpoly_ctx_t mctx, const nmodf_ctx_t fctx);

FLINT_DLL int nmod_mpoly_fprint_pretty(FILE * file, 
         const nmod_mpoly_t poly, const char ** x, const nmod_mpoly_ctx_t ctx);

NMOD_MPOLY_INLINE
int _nmod_mpoly_print_pretty(const ulong * poly, 
                const ulong * exps, slong len, const char ** x,
                    slong bits, const mpoly_ctx_t mctx, const nmodf_ctx_t fctx)
{
   return _nmod_mpoly_fprint_pretty(stdout, poly, exps, len, x, bits, mctx, fctx);
}

NMOD_MPOLY_INLINE
int nmod_mpoly_print_pretty(const nmod_mpoly_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
{
   return nmod_mpoly_fprint_pretty(stdout, poly, x, ctx);
}

/* Random generation *********************************************************/

FLINT_DLL void nmod_mpoly_randtest_bounds(nmod_mpoly_t poly, flint_rand_t state,
                 slong length, ulong * exp_bounds, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_randtest_bound(nmod_mpoly_t poly, flint_rand_t state,
                    slong length, ulong exp_bound, const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpoly_randtest_bits(nmod_mpoly_t poly, flint_rand_t state,
               slong length, mp_bitcnt_t exp_bits, const nmod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/






/******************************************************************************

   Internal consistency checks

******************************************************************************/

FLINT_DLL void nmod_mpoly_assert_canonical(const nmod_mpoly_t poly,
                                                   const nmod_mpoly_ctx_t ctx);

/*
   test that the terms in poly are in the correct order
*/
NMOD_MPOLY_INLINE
void nmod_mpoly_test(const nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
   slong i;

    if (mpoly_monomials_overflow_test(poly->exps, poly->length, poly->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

    if (!mpoly_monomials_inorder_test(poly->exps, poly->length, poly->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents out of order");

    for (i = 0; i < poly->length; i++)
    {
        if (poly->coeffs[i] == 0)
            flint_throw(FLINT_ERROR, "Polynomial has a zero coefficient");
    }
}


/*
   test that r is a valid remainder upon division by g
   this means that no monomial of r is divisible by lm(g)
*/
NMOD_MPOLY_INLINE
void nmod_mpoly_remainder_strongtest(const nmod_mpoly_t r, const nmod_mpoly_t g,
                                                    const nmod_mpoly_ctx_t ctx)
{
   slong i, N, bits;
   ulong mask = 0;
   ulong * rexp, * gexp;

   bits = FLINT_MAX(r->bits, g->bits);
   N = mpoly_words_per_exp(bits, ctx->minfo);

   if (g->length == 0 )
      flint_throw(FLINT_ERROR, "Zero denominator in remainder test");

   if (r->length == 0 )
      return;

   rexp = (ulong *) flint_malloc(N*r->length*sizeof(ulong));
   gexp = (ulong *) flint_malloc(N*1        *sizeof(ulong));
   mpoly_repack_monomials(rexp, bits, r->exps, r->bits, r->length, ctx->minfo);
   mpoly_repack_monomials(gexp, bits, g->exps, g->bits, 1,         ctx->minfo);

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
