/*
    Copyright (C) 2016-2017 William Hart

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
#include "mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef enum {
   ORD_LEX, ORD_REVLEX, ORD_DEGLEX, ORD_DEGREVLEX
} ordering_t;

typedef struct
{
   slong n;        /* number of elements in exponent vector (including deg) */
   ordering_t ord; /* polynomial ordering */
} fmpz_mpoly_ctx_struct;

typedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1];

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
} fmpz_mpoly_univariate_struct;

typedef fmpz_mpoly_univariate_struct fmpz_mpoly_univariate_t[1];

/* Context object ************************************************************/

FLINT_DLL void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, 
                                            slong nvars, const ordering_t ord);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx)
{
   /* nothing to be done at the moment */
}

/*  Memory management ********************************************************/

FLINT_DLL void fmpz_mpoly_init(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_init2(fmpz_mpoly_t poly, slong alloc, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_realloc(fmpz_mpoly_t poly, slong alloc, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_fit_length(fmpz_mpoly_t poly, slong len, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_clear(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_normalise(fmpz_mpoly_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

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

        _fmpz_mpoly_normalise(poly, ctx);
    }  
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_fit_bits(fmpz_mpoly_t poly,
                                        slong bits, const fmpz_mpoly_ctx_t ctx)
{
   if (bits > poly->bits)
   {
      slong N = (bits*ctx->n - 1)/FLINT_BITS + 1;

      poly->exps = flint_realloc(poly->exps, N*poly->alloc*sizeof(ulong));
      poly->bits = bits;
   }   
}

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

/*  Basic manipulation *******************************************************/

FLINT_DLL void _fmpz_mpoly_max_degrees(ulong * max_degs, const ulong * exps,
                    slong len, slong bits, slong n, int deg, int rev, slong N);

FLINT_DLL void fmpz_mpoly_max_degrees(ulong * max_degs,
                          const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_gen(fmpz * poly, ulong * exps, slong i,
                               slong bits, slong n, int deg, int rev, slong N);

FLINT_DLL void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i,
                                                   const fmpz_mpoly_ctx_t ctx);

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

FLINT_DLL int fmpz_mpoly_is_gen(const fmpz_mpoly_t poly,
                                          slong k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_coeff_fmpz(fmpz_t x,
                 const fmpz_mpoly_t poly, slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL ulong fmpz_mpoly_get_coeff_ui(const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mpoly_get_coeff_si(const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_fmpz(fmpz_mpoly_t poly, 
                          slong n, const fmpz_t x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_ui(fmpz_mpoly_t poly,
                                 slong n, ulong x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_si(fmpz_mpoly_t poly,
                                 slong n, slong x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_monomial(ulong * exps, const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_monomial(fmpz_mpoly_t poly, 
                      slong n, const ulong * exps, const fmpz_mpoly_ctx_t ctx);

#define fmpz_mpoly_get_coeff_ptr(poly, n, ctx) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

#define fmpz_mpoly_get_monomial_ptr(poly, n, ctx) \
    ((n) < (poly)->length ? (poly)->exps + \
                     (n)*(((ctx)->n - 1)/(FLINT_BITS/(poly)->bits) + 1) : NULL)

FLINT_DLL void _fmpz_mpoly_renormalise(fmpz_mpoly_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Set and negate */

FLINT_DLL void _fmpz_mpoly_set(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong N);

FLINT_DLL void fmpz_mpoly_set(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_neg(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong N);

FLINT_DLL void fmpz_mpoly_neg(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Comparison */

FLINT_DLL int _fmpz_mpoly_equal(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong N);

FLINT_DLL int fmpz_mpoly_equal(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
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
                 const fmpz * poly3, const ulong * exps3, slong len3, slong N);

FLINT_DLL void fmpz_mpoly_add(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                         const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_sub(fmpz * poly1, ulong * exps1,
                 const fmpz * poly2, const ulong * exps2, slong len2,
                 const fmpz * poly3, const ulong * exps3, slong len3, slong N);

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
                  const fmpz * poly3, const ulong * exp3, slong len3, slong N);

FLINT_DLL void fmpz_mpoly_mul_johnson(fmpz_mpoly_t poly1,
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
                                                 slong len2, slong k, slong N);

FLINT_DLL void fmpz_mpoly_pow_fps(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                          slong k, const fmpz_mpoly_ctx_t ctx);

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
                          const ulong * exp3, slong len3, slong bits, slong N);

FLINT_DLL int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Division ******************************************************************/

FLINT_DLL slong _fmpz_mpoly_div_monagan_pearce(fmpz ** polyq,
           ulong ** expq, slong * allocq, const fmpz * poly2,
   const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, 
                                slong len3, ulong * maxn, slong bits, slong N);

FLINT_DLL void fmpz_mpoly_div_monagan_pearce(fmpz_mpoly_t q,
                     const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_divrem_monagan_pearce(slong * lenr,
  fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr,
                  ulong ** expr, slong * allocr, const fmpz * poly2,
   const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, 
                                slong len3, ulong * maxn, slong bits, slong N);

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

/* Reduction *****************************************************************/

FLINT_DLL slong _fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** polyq, 
       fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2,
          const ulong * exp2, slong len2, const fmpz_mpoly_struct ** poly3,
                        ulong * const * exp3, slong len, slong N, slong bits,
                                     ulong * maxn, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** q, fmpz_mpoly_t r,
         const fmpz_mpoly_t poly2, const fmpz_mpoly_struct ** poly3, slong len,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Input/output **************************************************************/

FLINT_DLL char * _fmpz_mpoly_get_str_pretty(const fmpz * poly,
                          const ulong * exps, slong len, const char ** x, 
                               slong bits, slong n, int deg, int rev, slong N);

FLINT_DLL char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t poly,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_fprint_pretty(FILE * file, const fmpz * poly, 
                           const ulong * exps, slong len, const char ** x,
                               slong bits, slong n, int deg, int rev, slong N);

FLINT_DLL int fmpz_mpoly_fprint_pretty(FILE * file, 
         const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx);

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

#ifdef __cplusplus
}
#endif

#endif
