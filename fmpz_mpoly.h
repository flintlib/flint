/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MPOLY_H
#define FMPZ_MPOLY_H

#ifdef FMPZ_MPOLY_INLINES_C
#define FMPZ_MPOLY_INLINE FLINT_DLL
#else
#define FMPZ_MPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Context object ************************************************************/

typedef struct
{
    mpoly_ctx_t minfo;
} fmpz_mpoly_ctx_struct;

typedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1];

FLINT_DLL void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, 
                                            slong nvars, const ordering_t ord);

FLINT_DLL void fmpz_mpoly_ctx_init_rand(fmpz_mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars);


FLINT_DLL void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx);

/*  Type definitions *********************************************************/

typedef struct
{
   fmpz * coeffs; /* alloc fmpzs */
   ulong * exps;  
   slong alloc;
   slong length;
   slong bits;     /* number of bits per exponent */
} fmpz_mpoly_struct;

typedef fmpz_mpoly_struct fmpz_mpoly_t[1];

/* sparse univariates with multivariate coefficients */
typedef struct
{
   fmpz_mpoly_struct * coeffs; /* multivariate coefficients */
   ulong * exps;
   slong alloc;
   slong length;
   slong var; /* univariate variable number */
} fmpz_mpoly_univar_struct;

typedef fmpz_mpoly_univar_struct fmpz_mpoly_univar_t[1];





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
    fmpz * coeffs;
} fmpz_mpolyd_struct;

typedef fmpz_mpolyd_struct fmpz_mpolyd_t[1];

typedef struct
{
    slong nvars;
    slong * perm;
} fmpz_mpolyd_ctx_struct;

typedef fmpz_mpolyd_ctx_struct fmpz_mpolyd_ctx_t[1];





/* geobuckets ****************************************************************/
typedef struct fmpz_mpoly_geobucket
{
    fmpz_mpoly_struct polys[FLINT_BITS/2];
    slong length;
} fmpz_mpoly_geobucket_struct;

typedef fmpz_mpoly_geobucket_struct fmpz_mpoly_geobucket_t[1];

FLINT_DLL void fmpz_mpoly_geobucket_init(fmpz_mpoly_geobucket_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_clear(fmpz_mpoly_geobucket_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_empty(fmpz_mpoly_t p,
                         fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_print(fmpz_mpoly_geobucket_t B,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_fit_length(fmpz_mpoly_geobucket_t B,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_geobucket_fix(fmpz_mpoly_geobucket_t B, slong i,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_set(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_add(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_sub(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_set_fmpz(fmpz_mpoly_geobucket_t B,
                                         fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_gen(fmpz_mpoly_geobucket_t B, slong var,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_add_inplace(fmpz_mpoly_geobucket_t B1,
                        fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_sub_inplace(fmpz_mpoly_geobucket_t B1,
                        fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_neg_inplace(fmpz_mpoly_geobucket_t B1,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_mul_inplace(fmpz_mpoly_geobucket_t B1,
                        fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_pow_ui_inplace(fmpz_mpoly_geobucket_t B1,
                                          ulong k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_pow_fmpz_inplace(fmpz_mpoly_geobucket_t B1,
                                   const fmpz_t k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_geobucket_divides_inplace(fmpz_mpoly_geobucket_t B1,
                        fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

/*  Memory management ********************************************************/

FLINT_DLL void fmpz_mpoly_init(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_init2(fmpz_mpoly_t poly, slong alloc, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_realloc(fmpz ** poly, ulong ** exps,
                                            slong * alloc, slong len, slong N);

FLINT_DLL void fmpz_mpoly_realloc(fmpz_mpoly_t poly, slong alloc, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_fit_length(fmpz ** poly,
                             ulong ** exps, slong * alloc, slong len, slong N);

FLINT_DLL void fmpz_mpoly_fit_length(fmpz_mpoly_t poly, slong len, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_clear(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_set_length(fmpz_mpoly_t poly, slong newlen, 
                                                   const fmpz_mpoly_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; i++)
           _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = newlen;
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_truncate(fmpz_mpoly_t poly, slong newlen, 
                                                   const fmpz_mpoly_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        slong i;

        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);

        poly->length = newlen;
    }  
}

/*
   if poly->bits < bits, set poly->bits = bits and reallocate poly->exps
*/
FMPZ_MPOLY_INLINE
void fmpz_mpoly_fit_bits(fmpz_mpoly_t poly,
                                        slong bits, const fmpz_mpoly_ctx_t ctx)
{
   if (poly->bits < bits)
   {
      if (poly->alloc != 0)
      {
         slong N = mpoly_words_per_exp(bits, ctx->minfo);
         ulong * t = flint_malloc(N*poly->alloc*sizeof(ulong));
         mpoly_repack_monomials(t, bits, poly->exps,
                                         poly->bits, poly->length, ctx->minfo);
         flint_free(poly->exps);
         poly->exps = t;
      }

      poly->bits = bits;
   }
}

/*  Basic manipulation *******************************************************/

FMPZ_MPOLY_INLINE
int _fmpz_mpoly_fits_small(const fmpz * poly, slong len)
{
   slong i;
   for (i = 0; i < len; i++)
   {
      if (COEFF_IS_MPZ(poly[i]))
         return 0;
   }
   return 1;
}

FMPZ_MPOLY_INLINE
slong fmpz_mpoly_max_bits(const fmpz_mpoly_t poly)
{
    return _fmpz_vec_max_bits(poly->coeffs, poly->length);
}

FLINT_DLL void fmpz_mpoly_degrees_si(slong * degs, const fmpz_mpoly_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mpoly_degree_si(const fmpz_mpoly_t poly, slong var,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mpoly_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_degree_fmpz(fmpz_t degs, const fmpz_mpoly_t poly, slong var,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_gen(fmpz * poly, ulong * exps, slong i,
                               slong bits, slong n, int deg, int rev, slong N);

FLINT_DLL void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_is_gen(const fmpz_mpoly_t poly,
                                          slong k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_ui(fmpz_mpoly_t poly,
                                          ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_si(fmpz_mpoly_t poly,
                                          slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_fmpz(fmpz_mpoly_t poly,
                                   const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_equal_ui(const fmpz_mpoly_t poly,
                                          ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_equal_si(const fmpz_mpoly_t poly,
                                          slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_equal_fmpz(const fmpz_mpoly_t poly,
                                   const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_swap(fmpz_mpoly_t poly1, 
                                fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpoly_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_zero(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
   _fmpz_mpoly_set_length(poly, 0, ctx);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_one(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_ui(poly, UWORD(1), ctx);
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_is_zero(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
   return poly->length == 0;
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_is_one(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
   return fmpz_mpoly_equal_ui(poly, 1, ctx);
}


FLINT_DLL void _fmpz_mpoly_set_term_fmpz_fmpz(fmpz_mpoly_t poly,
                 const fmpz_t c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_term_fmpz_fmpz(fmpz_mpoly_t poly,
                const fmpz_t c, const fmpz ** exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_term_ui_fmpz(fmpz_mpoly_t poly,
                 const ulong c, const fmpz ** exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_term_si_fmpz(fmpz_mpoly_t poly,
                 const slong c, const fmpz ** exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_term_fmpz_ui(fmpz_mpoly_t poly,
                const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_term_ui_ui(fmpz_mpoly_t poly,
                 const ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_term_si_ui(fmpz_mpoly_t poly,
                 const slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_get_term_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t poly,
                                 const fmpz * exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_get_term_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t poly,
                                const fmpz ** exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL ulong fmpz_mpoly_get_term_ui_fmpz(           const fmpz_mpoly_t poly,
                                const fmpz ** exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL slong fmpz_mpoly_get_term_si_fmpz(           const fmpz_mpoly_t poly,
                                const fmpz ** exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_get_term_fmpz_ui(fmpz_t c, const fmpz_mpoly_t poly,
                                const ulong * exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL ulong fmpz_mpoly_get_term_ui_ui(           const fmpz_mpoly_t poly,
                                const ulong * exp, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL slong fmpz_mpoly_get_term_si_ui(           const fmpz_mpoly_t poly,
                                const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_coeff_fmpz(fmpz_t x, const fmpz_mpoly_t poly,
                                          slong n, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL ulong fmpz_mpoly_get_coeff_ui(           const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL slong fmpz_mpoly_get_coeff_si(           const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_fmpz(fmpz_mpoly_t poly, 
                          slong n, const fmpz_t x, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_coeff_ui(fmpz_mpoly_t poly,
                                 slong n, ulong x, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_coeff_si(fmpz_mpoly_t poly,
                                 slong n, slong x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_monomial_ui(ulong * exps, const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_get_monomial_fmpz(fmpz ** exps, const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_monomial_ui(fmpz_mpoly_t poly, 
                      slong n, const ulong * exps, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_set_monomial_fmpz(fmpz_mpoly_t poly, 
                      slong n,       fmpz ** exps, const fmpz_mpoly_ctx_t ctx);

#define fmpz_mpoly_get_coeff_ptr(poly, n, ctx) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

#define fmpz_mpoly_get_monomial_ptr(poly, n, ctx) \
    ((n) < (poly)->length ? (poly)->exps + \
                     (n)*(((ctx)->n - 1)/(FLINT_BITS/(poly)->bits) + 1) : NULL)


/* Set and negate ************************************************************/

FLINT_DLL void _fmpz_mpoly_set(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong N);

FLINT_DLL void fmpz_mpoly_set(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_neg(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong N);

FLINT_DLL void fmpz_mpoly_neg(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Comparison ****************************************************************/

FLINT_DLL int _fmpz_mpoly_equal(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong N);

FLINT_DLL int fmpz_mpoly_equal(const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Reverse *******************************************************************/

FLINT_DLL void _fmpz_mpoly_reverse(fmpz * poly1, ulong * exp1,
                   const fmpz * poly2, const ulong * exp2, slong len, slong N);

FLINT_DLL void fmpz_mpoly_reverse(fmpz_mpoly_t poly1,
                               fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);

/* Basic arithmetic **********************************************************/

FLINT_DLL void fmpz_mpoly_add_ui(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_add_si(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_add_fmpz(fmpz_mpoly_t poly1,
         const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_ui(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_si(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_fmpz(fmpz_mpoly_t poly1,
         const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_add(fmpz * poly1, ulong * exps1,
                 const fmpz * poly2, const ulong * exps2, slong len2,
                 const fmpz * poly3, const ulong * exps3, slong len3, slong N,
                                                        const ulong * cmpmask);

FLINT_DLL void fmpz_mpoly_add(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                         const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_sub(fmpz * poly1, ulong * exps1,
                 const fmpz * poly2, const ulong * exps2, slong len2,
                 const fmpz * poly3, const ulong * exps3, slong len3, slong N,
                                                        const ulong * cmpmask);

FLINT_DLL void fmpz_mpoly_sub(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                         const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

/* Scalar operations *********************************************************/

FLINT_DLL void _fmpz_mpoly_scalar_mul_fmpz(fmpz * poly1, ulong * exps1,
 const fmpz * poly2, const ulong * exps2, slong len2, slong N, const fmpz_t c);

FLINT_DLL void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t poly1,
         const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_mul_si(fmpz * poly1, ulong * exps1,
        const fmpz * poly2, const ulong * exps2, slong len2, slong N, slong c);

FLINT_DLL void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_mul_ui(fmpz * poly1, ulong * exps1,
        const fmpz * poly2, const ulong * exps2, slong len2, slong N, ulong c);

FLINT_DLL void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_divexact_fmpz(fmpz * poly1, ulong * exps1,
 const fmpz * poly2, const ulong * exps2, slong len2, slong N, const fmpz_t c);

FLINT_DLL void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t poly1,
         const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_divexact_si(fmpz * poly1, ulong * exps1,
        const fmpz * poly2, const ulong * exps2, slong len2, slong N, slong c);

FLINT_DLL void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_divexact_ui(fmpz * poly1, ulong * exps1,
        const fmpz * poly2, const ulong * exps2, slong len2, slong N, ulong c);

FLINT_DLL void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx);

/* Multiplication ************************************************************/

FLINT_DLL slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1,
        slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                const fmpz * poly3, const ulong * exp3, slong len3, mp_bitcnt_t bits, slong N,
                                                        const ulong * cmpmask);

FLINT_DLL void fmpz_mpoly_mul_johnson(fmpz_mpoly_t poly1,
                 const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_mul_heap_threaded(fmpz_mpoly_t poly1,
                 const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_mul_array(fmpz ** poly1, ulong ** exp1,
          slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, 
                         const fmpz * poly3, const ulong * exp3, slong len3, 
                                         slong * mults, slong num, slong bits);

FLINT_DLL int fmpz_mpoly_mul_array(fmpz_mpoly_t poly1, 
                 const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

FLINT_DLL slong _fmpz_mpoly_pow_fps(fmpz ** poly1, ulong ** exp1,
                slong * alloc, const fmpz * poly2, const ulong * exp2, 
        slong len2, ulong k, mp_bitcnt_t bits, slong N, const ulong * cmpmask);

FLINT_DLL void fmpz_mpoly_pow_fps(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                          ulong k, const fmpz_mpoly_ctx_t ctx);
FLINT_DLL void fmpz_mpoly_pow_fmpz(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                 const fmpz_t pow, const fmpz_mpoly_ctx_t ctx);
/* Calculus ******************************************************************/

FLINT_DLL void fmpz_mpoly_derivative(fmpz_mpoly_t poly1,
              const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_integral(fmpz_mpoly_t poly1, fmpz_t scale,
              const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx);

/* Divisibility **************************************************************/

FLINT_DLL slong _fmpz_mpoly_divides_array(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                        const fmpz * poly3, const ulong * exp3, slong len3,
                                         slong * mults, slong num, slong bits);

FLINT_DLL int fmpz_mpoly_divides_array(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1,
                      ulong ** exp1, slong * alloc, const fmpz * poly2,
                    const ulong * exp2, slong len2, const fmpz * poly3,
                          const ulong * exp3, slong len3, mp_bitcnt_t bits, slong N,
                                                         const ulong * cmpmask);

FLINT_DLL int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Division ******************************************************************/

FLINT_DLL slong _fmpz_mpoly_div_monagan_pearce(fmpz ** polyq,
           ulong ** expq, slong * allocq, const fmpz * poly2,
   const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, 
                       slong len3, slong bits, slong N, const ulong * cmpmask);

FLINT_DLL void fmpz_mpoly_div_monagan_pearce(fmpz_mpoly_t q,
                     const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_divrem_monagan_pearce(slong * lenr,
  fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr,
                  ulong ** expr, slong * allocr, const fmpz * poly2,
   const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, 
                       slong len3, slong bits, slong N, const ulong * cmpmask);

FLINT_DLL void fmpz_mpoly_divrem_monagan_pearce(fmpz_mpoly_t q, fmpz_mpoly_t r,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_divrem_array(slong * lenr,
       fmpz ** polyq, ulong ** expq, slong * allocq,
              fmpz ** polyr, ulong ** expr, slong * allocr, 
                const fmpz * poly2, const ulong * exp2, slong len2, 
        const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, 
                                                        slong num, slong bits);

FLINT_DLL int fmpz_mpoly_divrem_array(fmpz_mpoly_t q, fmpz_mpoly_t r,
                    const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_quasidivrem_heap(fmpz_t scale,
                        fmpz_mpoly_t q, fmpz_mpoly_t r,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_quasidiv_heap(fmpz_t scale, fmpz_mpoly_t q,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Evaluation ****************************************************************/

FLINT_DLL void fmpz_mpoly_evaluate_all_tree_fmpz(fmpz_t ev, fmpz_mpoly_t poly,
                                            fmpz ** val, fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t poly1, fmpz_mpoly_t poly2,
                                  slong var, fmpz_t val, fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_compose(fmpz_mpoly_t res, fmpz_mpoly_t poly1,
    fmpz_mpoly_struct ** polys2, fmpz_mpoly_ctx_t ctx1, fmpz_mpoly_ctx_t ctx2);

/* Univariates ***************************************************************/

FLINT_DLL void fmpz_mpoly_univar_init(fmpz_mpoly_univar_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_clear(fmpz_mpoly_univar_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_swap(fmpz_mpoly_univar_t poly1,
                        fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_fit_length(fmpz_mpoly_univar_t poly,
                                     slong length, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_print_pretty(const fmpz_mpoly_univar_t poly,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_assert_canonical(fmpz_mpoly_univar_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_from_univar(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_to_univar(fmpz_mpoly_univar_t poly1,
              const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_univar_equal(fmpz_mpoly_univar_t poly1,
                  const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_add(fmpz_mpoly_univar_t poly1,
            const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_univar_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_univar_mul(fmpz_mpoly_univar_t poly1,
            const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_univar_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_to_fmpz_poly(fmpz_poly_t poly1, slong * shift1,
              const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_from_fmpz_poly(fmpz_mpoly_t poly1,
        const fmpz_poly_t poly2, slong shift2, slong var,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_derivative(fmpz_mpoly_univar_t poly1,
                  const fmpz_mpoly_univar_t poly2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_univar_prem(fmpz_mpoly_univar_t polyA,
            const fmpz_mpoly_univar_t polyB, fmpz_mpoly_univar_t polyC,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_univar_pgcd(fmpz_mpoly_univar_t poly1,
            const fmpz_mpoly_univar_t polyP, const fmpz_mpoly_univar_t polyQ,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_univar_pgcd_ducos(fmpz_mpoly_univar_t poly1,
            const fmpz_mpoly_univar_t polyP, const fmpz_mpoly_univar_t polyQ,
                                                   const fmpz_mpoly_ctx_t ctx);

/* GCD ***********************************************************************/

FLINT_DLL void fmpz_mpoly_term_content(fmpz_mpoly_t poly1,
                         const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_prs(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                         const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_is_unit(const fmpz_mpoly_t a, const fmpz_mpoly_t b,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_resultant(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                        slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_discriminant(fmpz_mpoly_t poly1,
              const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_brown(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                         const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

/* Reduction *****************************************************************/

FLINT_DLL slong
_fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** polyq, 
       fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2,
          const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3,
                        ulong * const * exp3, slong len, slong N, slong bits,
                            const fmpz_mpoly_ctx_t ctx, const ulong * cmpmask);

FLINT_DLL void
fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** q, fmpz_mpoly_t r,
    const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void
fmpz_mpoly_quasidivrem_ideal_heap(fmpz_t scale,
                                 fmpz_mpoly_struct ** q, fmpz_mpoly_t r,
                const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3,
                                        slong len, const fmpz_mpoly_ctx_t ctx);


/* Input/output **************************************************************/

FLINT_DLL int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t poly, const char * str,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL char * _fmpz_mpoly_get_str_pretty(const fmpz * poly,
                          const ulong * exps, slong len, const char ** x, 
                                           slong bits, const mpoly_ctx_t mctx);

FLINT_DLL char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t poly,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_fprint_pretty(FILE * file, const fmpz * poly, 
                        const ulong * exps, slong len, const char ** x_in,
                                     mp_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL int fmpz_mpoly_fprint_pretty(FILE * file, 
         const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
int _fmpz_mpoly_print_pretty(const fmpz * poly, 
                       const ulong * exps, slong len, const char ** x,
                                            slong bits, const mpoly_ctx_t mctx)
{
    return _fmpz_mpoly_fprint_pretty(stdout, poly, exps, len, x, bits, mctx);
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_print_pretty(const fmpz_mpoly_t poly,
                                   const char ** x, const fmpz_mpoly_ctx_t ctx)
{
   return fmpz_mpoly_fprint_pretty(stdout, poly, x, ctx);
}

/* Random generation *********************************************************/

FLINT_DLL void fmpz_mpoly_randtest_bound(fmpz_mpoly_t poly, flint_rand_t state,
   slong length, mp_bitcnt_t coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_randtest_bits(fmpz_mpoly_t poly, flint_rand_t state,
   slong length, mp_bitcnt_t coeff_bits, mp_bitcnt_t exp_bits, const fmpz_mpoly_ctx_t ctx);


/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

/* Internal packing and conversion */

FLINT_DLL slong _fmpz_mpoly_from_ulong_array(fmpz ** poly1,
                         ulong ** exp1, slong * alloc, ulong * poly2,
                          const slong * mults, slong num, slong bits, slong k);

FLINT_DLL slong _fmpz_mpoly_from_ulong_array2(fmpz ** poly1,
                         ulong ** exp1, slong * alloc, ulong * poly2, 
                          const slong * mults, slong num, slong bits, slong k);

FLINT_DLL slong _fmpz_mpoly_from_ulong_array1(fmpz ** poly1,
                         ulong ** exp1, slong * alloc, ulong * poly2,
                          const slong * mults, slong num, slong bits, slong k);

FLINT_DLL slong _fmpz_mpoly_from_fmpz_array(fmpz ** poly1,
                         ulong ** exp1, slong * alloc, fmpz * poly2,
                          const slong * mults, slong num, slong bits, slong k);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong(ulong * poly1, 
                  const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong2(ulong * poly1, 
                  const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong1(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_fmpz(fmpz * poly1, 
                 const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong_1(ulong * poly1, 
                          slong d, const ulong exp2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_slong2_1(ulong * poly1, 
                           slong d, const ulong exp2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_submul_array1_fmpz_1(fmpz * poly1, 
                          const fmpz_t d, ulong exp2,
                           const fmpz * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_to_ulong_array2(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_ulong_array1(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_ulong_array(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_fmpz_array(fmpz * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_chunk_max_bits(slong * b1, slong * maxb1,
                          const fmpz * poly1, slong * i1, slong * n1, slong i);

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_sub_uiuiui_fmpz(ulong * c, const fmpz_t d)
{
   fmpz fc = *d;

   if (!COEFF_IS_MPZ(fc))
   {
        ulong f0, f1, f2;
        f0 = fc;
        f1 = f2 = FLINT_SIGN_EXT(f0);
        sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], f2, f1, f0);
   } else
   {
      slong size = fmpz_size(d);
      __mpz_struct * m = COEFF_TO_PTR(fc);
      if (fmpz_sgn(d) < 0)
         mpn_add(c, c, 3, m->_mp_d, size);
      else
         mpn_sub(c, c, 3, m->_mp_d, size);
   }
}

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_add_uiuiui_fmpz(ulong * c, const fmpz_t d)
{
    fmpz fc = *d;

    if (!COEFF_IS_MPZ(fc))
    {
        ulong f0, f1, f2;
        f0 = fc;
        f1 = f2 = FLINT_SIGN_EXT(f0);
        add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], f2, f1, f0);
   } else
   {
      slong size = fmpz_size(d);
      __mpz_struct * m = COEFF_TO_PTR(fc);
      if (fmpz_sgn(d) < 0)
         mpn_sub(c, c, 3, m->_mp_d, size);
      else
         mpn_add(c, c, 3, m->_mp_d, size);
   }
}

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_submul_uiuiui_fmpz(ulong * c, slong d1, slong d2)
{
    ulong p[2], p2;
    smul_ppmm(p[1], p[0], d1, d2);
    p2 = FLINT_SIGN_EXT(p[1]);
    sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], p2, p[1], p[0]);
}

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_addmul_uiuiui_fmpz(ulong * c, slong d1, slong d2)
{
    ulong p[2], p2;
    smul_ppmm(p[1], p[0], d1, d2);
    p2 = FLINT_SIGN_EXT(p[1]);
    add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], p2, p[1], p[0]);
}



/******************************************************************************

   Internal consistency checks

******************************************************************************/

FLINT_DLL void fmpz_mpoly_assert_canonical(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

/*
   test that r is a valid remainder upon division by g
   this means that if c*x^a is a term of r and x^a is divisible by the leading
   monomial of g, then |c| < |leading coefficient of g|
*/
FMPZ_MPOLY_INLINE
void fmpz_mpoly_remainder_test(const fmpz_mpoly_t r, const fmpz_mpoly_t g,
                                                    const fmpz_mpoly_ctx_t ctx)
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
      if (mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask)
         && fmpz_cmpabs(g->coeffs + 0, r->coeffs + i) <= 0)
      {
         flint_printf("fmpz_mpoly_remainder_test FAILED i = %wd\n", i);
         flint_printf("rem ");fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
         flint_printf("den ");fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
         flint_abort();
      }

   flint_free(rexp);
   flint_free(gexp);
}


/*
   test that r is a valid remainder upon division by g over Q
   this means that no term of r is divisible by lt(g)
*/
FMPZ_MPOLY_INLINE
void fmpz_mpoly_remainder_strongtest(const fmpz_mpoly_t r, const fmpz_mpoly_t g,
                                                    const fmpz_mpoly_ctx_t ctx)
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
         flint_printf("fmpz_mpoly_remainder_strongtest FAILED i = %wd\n", i);
         flint_printf("rem ");fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
         flint_printf("den ");fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
         flint_abort();
      }

   flint_free(rexp);
   flint_free(gexp);
}



#ifdef __cplusplus
}
#endif

#endif
