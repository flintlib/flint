/*
    Copyright (C) 2016-2017 William Hart
    Copyright (C) 2017-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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
#include "nmod_mpoly.h"
#include "fmpz_mod.h"
#include "n_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Type definitions **********************************************************/

/*
    context object for fmpz_mpoly
*/
typedef struct
{
    mpoly_ctx_t minfo;
} fmpz_mpoly_ctx_struct;

typedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1];

/*
    fmpz_mpoly_t
    sparse multivariates with fmpz coeffs
*/
typedef struct
{
   fmpz * coeffs; /* alloc fmpzs */
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;     /* number of bits per exponent */
} fmpz_mpoly_struct;

typedef fmpz_mpoly_struct fmpz_mpoly_t[1];

FMPZ_MPOLY_INLINE
fmpz * fmpz_mpoly_term_coeff_ref(fmpz_mpoly_t A, slong i,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < A->length);
    return A->coeffs + i;
}

FMPZ_MPOLY_INLINE fmpz * fmpz_mpoly_leadcoeff(const fmpz_mpoly_t A)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

/* Internal type definitions *************************************************/

/*
    fmpz_mpoly_univar_t
    sparse univariates with multivariate coefficients
*/
typedef struct
{
   fmpz_mpoly_struct * coeffs; /* multivariate coefficients */
   fmpz * exps;
   slong alloc;
   slong length;
} fmpz_mpoly_univar_struct;

typedef fmpz_mpoly_univar_struct fmpz_mpoly_univar_t[1];

/*
    fmpz_mpolyd_t
    A dense mpoly is stored as a flat array of coeffcients.
    Suppose deg_bounds = {r0, r1, r2}. The coefficient of the monomial with
    exponents {e0, e1, e2} is stored at the coefficient of index
        e2 + r2*(e1 + r1*(e0 + r0*0))
*/
typedef struct
{
    slong nvars;
    slong degb_alloc;
    slong * deg_bounds;
    slong length;           /* usage is inconsistent currently */
    slong coeff_alloc;
    fmpz * coeffs;
} fmpz_mpolyd_struct;

typedef fmpz_mpolyd_struct fmpz_mpolyd_t[1];

/* Context object ************************************************************/

FLINT_DLL void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, 
                                            slong nvars, const ordering_t ord);

FLINT_DLL void fmpz_mpoly_ctx_init_rand(fmpz_mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars);


FLINT_DLL void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
slong fmpz_mpoly_ctx_nvars(const fmpz_mpoly_ctx_t ctx)
{
    return ctx->minfo->nvars;
}

FMPZ_MPOLY_INLINE
ordering_t fmpz_mpoly_ctx_ord(const fmpz_mpoly_ctx_t ctx)
{
    return ctx->minfo->ord;
}


/*  Memory management ********************************************************/

FLINT_DLL void fmpz_mpoly_init(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_init2(fmpz_mpoly_t A, slong alloc, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_init3(fmpz_mpoly_t A, slong alloc, flint_bitcnt_t bits,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_realloc(fmpz ** Acoeff, ulong ** Aexp,
                                           slong * Aalloc, slong len, slong N);

FLINT_DLL void fmpz_mpoly_realloc(fmpz_mpoly_t A, slong alloc, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_fit_length(fmpz ** Acoeff,
                            ulong ** Aexp, slong * Aalloc, slong len, slong N);

FLINT_DLL void fmpz_mpoly_fit_length(fmpz_mpoly_t A, slong len, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_fit_length_reset_bits(fmpz_mpoly_t A, slong len,
                              flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_clear(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_set_length(fmpz_mpoly_t A, slong newlen, 
                                                   const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    for (i = newlen; i < A->length; i++)
       _fmpz_demote(A->coeffs + i);

    A->length = newlen;
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_truncate(fmpz_mpoly_t A, slong newlen, 
                                                   const fmpz_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        slong i;

        for (i = newlen; i < A->length; i++)
            _fmpz_demote(A->coeffs + i);

        A->length = newlen;
    }  
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_fit_bits(fmpz_mpoly_t A,
                                  flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)
{
   if (A->bits < bits)
   {
      if (A->alloc != 0)
      {
         slong N = mpoly_words_per_exp(bits, ctx->minfo);
         ulong * t = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
         mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length, ctx->minfo);
         flint_free(A->exps);
         A->exps = t;
      }

      A->bits = bits;
   }
}


/* Input/output **************************************************************/

FLINT_DLL int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t A, const char * str,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL char * _fmpz_mpoly_get_str_pretty(const fmpz * poly,
                          const ulong * exps, slong len, const char ** x, 
                                           slong bits, const mpoly_ctx_t mctx);

FLINT_DLL char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t A,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_fprint_pretty(FILE * file, const fmpz * poly, 
                        const ulong * exps, slong len, const char ** x_in,
                                     flint_bitcnt_t bits, const mpoly_ctx_t mctx);

FLINT_DLL int fmpz_mpoly_fprint_pretty(FILE * file, 
            const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
int _fmpz_mpoly_print_pretty(const fmpz * poly, 
                       const ulong * exps, slong len, const char ** x,
                                            slong bits, const mpoly_ctx_t mctx)
{
    return _fmpz_mpoly_fprint_pretty(stdout, poly, exps, len, x, bits, mctx);
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_print_pretty(const fmpz_mpoly_t A,
                                   const char ** x, const fmpz_mpoly_ctx_t ctx)
{
   return fmpz_mpoly_fprint_pretty(stdout, A, x, ctx);
}


/*  Basic manipulation *******************************************************/

FLINT_DLL void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_is_gen(const fmpz_mpoly_t poly,
                                          slong k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_set(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong N);

FLINT_DLL void fmpz_mpoly_set(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_equal(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong N);

FLINT_DLL int fmpz_mpoly_equal(const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_swap(fmpz_mpoly_t A, 
                                fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpoly_struct t = *A;
   *A = *B;
   *B = t;
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
slong fmpz_mpoly_max_bits(const fmpz_mpoly_t A)
{
    return _fmpz_vec_max_bits(A->coeffs, A->length);
}

/* Constants *****************************************************************/

FLINT_DLL int fmpz_mpoly_is_fmpz(const fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_fmpz(fmpz_t c, const fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_fmpz(fmpz_mpoly_t A,
                                   const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_ui(fmpz_mpoly_t A,
                                          ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_si(fmpz_mpoly_t A,
                                          slong c, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_zero(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
   _fmpz_mpoly_set_length(A, 0, ctx);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_one(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_ui(A, UWORD(1), ctx);
}

FLINT_DLL int fmpz_mpoly_equal_fmpz(const fmpz_mpoly_t A,
                                   const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_equal_ui(const fmpz_mpoly_t A,
                                          ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_equal_si(const fmpz_mpoly_t A,
                                          slong c, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
int fmpz_mpoly_is_zero(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
   return A->length == 0;
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_is_one(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
   return fmpz_mpoly_equal_ui(A, UWORD(1), ctx);
}


/* Degrees *******************************************************************/

FMPZ_MPOLY_INLINE
int fmpz_mpoly_degrees_fit_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_degrees_si(slong * degs, const fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mpoly_t A, slong var,
                                                   const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

FMPZ_MPOLY_INLINE
slong fmpz_mpoly_degree_si(const fmpz_mpoly_t A, slong var,
                                                   const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_total_degree_fits_si(const fmpz_mpoly_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_total_degree_fmpz(fmpz_t td, const fmpz_mpoly_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MPOLY_INLINE
slong fmpz_mpoly_total_degree_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_used_vars(int * used, const fmpz_mpoly_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < ctx->minfo->nvars; i++)
        used[i] = 0;

    mpoly_used_vars_or(used, A->exps, A->length, A->bits, ctx->minfo);
}

/* Coefficients **************************************************************/

FLINT_DLL void fmpz_mpoly_get_coeff_fmpz_monomial(fmpz_t c,
                          const fmpz_mpoly_t A, const fmpz_mpoly_t M,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_fmpz_monomial(fmpz_mpoly_t A,
                                  const fmpz_t c, const fmpz_mpoly_t M,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t A,
                               fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL ulong fmpz_mpoly_get_coeff_ui_fmpz(           const fmpz_mpoly_t A,
                               fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mpoly_get_coeff_si_fmpz(           const fmpz_mpoly_t A,
                               fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mpoly_t A,
                                const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL ulong fmpz_mpoly_get_coeff_ui_ui(           const fmpz_mpoly_t A,
                                const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mpoly_get_coeff_si_ui(           const fmpz_mpoly_t A,
                                const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A,
                 const fmpz_t c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A,
               const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_ui_fmpz(fmpz_mpoly_t A,
                const ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_si_fmpz(fmpz_mpoly_t A,
                const slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_t A,
                const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_ui_ui(fmpz_mpoly_t A,
                 const ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_si_ui(fmpz_mpoly_t A,
                 const slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_coeff_vars_ui(fmpz_mpoly_t C,
             const fmpz_mpoly_t A, const slong * vars, const ulong * exps,
                                     slong length, const fmpz_mpoly_ctx_t ctx);

/* conversion ****************************************************************/

FLINT_DLL int fmpz_mpoly_is_fmpz_poly(const fmpz_mpoly_t A, slong var,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_get_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B,
                                        slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_t A, flint_bitcnt_t Abits,
      const fmpz * Bcoeffs, slong Blen, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_t A, const fmpz_poly_t B,
                                          slong v, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_set_fmpz_poly_one_var(fmpz_mpoly_t A,
                        flint_bitcnt_t Aminbits, fmpz * Acoeffs, slong Adeg,
                                                   const fmpz_mpoly_ctx_t ctx);

/* comparison ****************************************************************/

FLINT_DLL int fmpz_mpoly_cmp(const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

/* container operations ******************************************************/

FLINT_DLL int fmpz_mpoly_is_canonical(const fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
slong fmpz_mpoly_length(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    return A->length;
}

FLINT_DLL void fmpz_mpoly_resize(fmpz_mpoly_t A, slong new_length,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mpoly_t A,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL ulong fmpz_mpoly_get_term_coeff_ui(           const fmpz_mpoly_t A, 
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mpoly_get_term_coeff_si(           const fmpz_mpoly_t A, 
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_term_coeff_fmpz(fmpz_mpoly_t A, 
                          slong i, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_term_coeff_ui(fmpz_mpoly_t A,
                                 slong i, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_term_coeff_si(fmpz_mpoly_t A,
                                 slong i, slong c, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
int fmpz_mpoly_term_exp_fits_ui(const fmpz_mpoly_t A, slong i,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_ui(A->exps, A->bits, i, ctx->minfo);
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_term_exp_fits_si(const fmpz_mpoly_t A, slong i,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_si(A->exps, A->bits, i, ctx->minfo);
}

FLINT_DLL void fmpz_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mpoly_t A,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mpoly_t A,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_term_exp_si(slong * exp, const fmpz_mpoly_t A,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL ulong fmpz_mpoly_get_term_var_exp_ui(const fmpz_mpoly_t A, slong i,
                                        slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mpoly_get_term_var_exp_si(const fmpz_mpoly_t A, slong i,
                                        slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_term_exp_fmpz(fmpz_mpoly_t A,
                      slong i, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_term_exp_ui(fmpz_mpoly_t A,
                       slong i, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_term(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_term_monomial(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_push_term_fmpz_fmpz(fmpz_mpoly_t A,
               const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_push_term_ui_fmpz(fmpz_mpoly_t A,
                      ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_push_term_si_fmpz(fmpz_mpoly_t A,
                      slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_t A,
                const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_push_term_ui_ui(fmpz_mpoly_t A,
                       ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_push_term_si_ui(fmpz_mpoly_t A,
                       slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sort_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_combine_like_terms(fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_reverse(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_assert_canonical(const fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_radix_sort1(fmpz_mpoly_t A, slong left, slong right,
                              flint_bitcnt_t pos, ulong cmpmask, ulong totalmask);

FLINT_DLL void _fmpz_mpoly_radix_sort(fmpz_mpoly_t A, slong left, slong right,
                                    flint_bitcnt_t pos, slong N, ulong * cmpmask);

FLINT_DLL void _fmpz_mpoly_push_exp_ffmpz(fmpz_mpoly_t A,
                                 const fmpz * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_push_exp_pfmpz(fmpz_mpoly_t A,
                               fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_push_exp_ui(fmpz_mpoly_t A,
                                const ulong * exp, const fmpz_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

FLINT_DLL void fmpz_mpoly_randtest_bound(fmpz_mpoly_t A, flint_rand_t state,
                        slong length, flint_bitcnt_t coeff_bits, ulong exp_bound,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_randtest_bounds(fmpz_mpoly_t A, flint_rand_t state,
                     slong length, flint_bitcnt_t coeff_bits, ulong * exp_bounds,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_randtest_bits(fmpz_mpoly_t A, flint_rand_t state,
                   slong length, flint_bitcnt_t coeff_bits, flint_bitcnt_t exp_bits,
                                                   const fmpz_mpoly_ctx_t ctx);


/* Addition/Subtraction ******************************************************/

FLINT_DLL void fmpz_mpoly_add_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                   const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_add_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                          ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_add_si(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                          slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                   const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                          ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_si(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                          slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_add(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                             const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_add(fmpz * poly1, ulong * exps1,
                 const fmpz * poly2, const ulong * exps2, slong len2,
                 const fmpz * poly3, const ulong * exps3, slong len3, slong N,
                                                        const ulong * cmpmask);

FLINT_DLL void fmpz_mpoly_sub(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                             const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_sub(fmpz * poly1, ulong * exps1,
                 const fmpz * poly2, const ulong * exps2, slong len2,
                 const fmpz * poly3, const ulong * exps3, slong len3, slong N,
                                                        const ulong * cmpmask);


/* Scalar operations *********************************************************/

FLINT_DLL void fmpz_mpoly_neg(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t A,
             const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t A,
                    const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t A,
                    const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_scalar_fmma(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                         const fmpz_t c, const fmpz_mpoly_t D, const fmpz_t e,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t A,
             const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t A,
                    const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t A,
                    const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_scalar_divides_fmpz(fmpz_mpoly_t A,
             const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_scalar_divides_si(fmpz_mpoly_t A,
                    const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_scalar_divides_ui(fmpz_mpoly_t A,
                    const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx);


/* Differentiation/Integration ***********************************************/

FLINT_DLL void fmpz_mpoly_derivative(fmpz_mpoly_t A,
                  const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_integral(fmpz_mpoly_t A, fmpz_t scale,
                  const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx);


/* Evaluation ****************************************************************/

FLINT_DLL int _fmpz_pow_ui_is_not_feasible(flint_bitcnt_t bbits, ulong e);

FLINT_DLL int _fmpz_pow_fmpz_is_not_feasible(flint_bitcnt_t bbits, const fmpz_t e);

FLINT_DLL int fmpz_mpoly_evaluate_all_fmpz(fmpz_t ev, const fmpz_mpoly_t A,
                              fmpz * const * vals, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL mp_limb_t fmpz_mpoly_evaluate_all_nmod(const fmpz_mpoly_t A,
           const mp_limb_t * alphas, const fmpz_mpoly_ctx_t ctx, nmod_t fpctx);

FLINT_DLL void fmpz_mpoly_evaluate_all_fmpz_mod(fmpz_t ev,
                        const fmpz_mpoly_t A, const fmpz * alphas,
                       const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx);

FLINT_DLL int fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t A,
                           const fmpz_mpoly_t B, slong var, const fmpz_t val,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A,
                         const fmpz_mpoly_t B, fmpz_poly_struct * const * C,
                                                  const fmpz_mpoly_ctx_t ctxB);

FLINT_DLL void _fmpz_mpoly_compose_mat(fmpz_mpoly_t A,
                            const fmpz_mpoly_t B, const fmpz_mat_t M,
                    const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);

FLINT_DLL int fmpz_mpoly_compose_fmpz_mpoly_geobucket(fmpz_mpoly_t A,
                   const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C,
                    const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);

FLINT_DLL int fmpz_mpoly_compose_fmpz_mpoly_horner(fmpz_mpoly_t A,
                   const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C,
                    const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);

FLINT_DLL int fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_t A,
                   const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C,
                    const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);

FLINT_DLL void fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_t A,
                             const fmpz_mpoly_t B, const slong * c,
                    const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);


/* Multiplication ************************************************************/

FLINT_DLL void fmpz_mpoly_mul(fmpz_mpoly_t A,
       const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_mul_monomial(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                             const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_mul_johnson(fmpz_mpoly_t A,
       const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_mul_heap_threaded(fmpz_mpoly_t A,
       const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_mul_array(fmpz_mpoly_t A, 
       const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_mul_array_threaded(fmpz_mpoly_t A,
       const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_mul_dense(fmpz_mpoly_t A, 
       const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
                             flint_bitcnt_t bits, slong N, const ulong * cmpmask);

FLINT_DLL void _fmpz_mpoly_mul_johnson_maxfields(fmpz_mpoly_t A,
                                 const fmpz_mpoly_t B, fmpz * maxBfields,
                                 const fmpz_mpoly_t C, fmpz * maxCfields,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_mul_heap_threaded_pool_maxfields(fmpz_mpoly_t A,
           const fmpz_mpoly_t B, fmpz * maxBfields,
           const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);

FLINT_DLL int _fmpz_mpoly_mul_array_DEG(fmpz_mpoly_t A,
                                 const fmpz_mpoly_t B, fmpz * maxBfields,
                                 const fmpz_mpoly_t C, fmpz * maxCfields,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_mul_array_LEX(fmpz_mpoly_t A,
                                 const fmpz_mpoly_t B, fmpz * maxBfields,
                                 const fmpz_mpoly_t C, fmpz * maxCfields,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_mul_array_threaded_pool_DEG(fmpz_mpoly_t A,
           const fmpz_mpoly_t B, fmpz * maxBfields,
           const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);

FLINT_DLL int _fmpz_mpoly_mul_array_threaded_pool_LEX(fmpz_mpoly_t A,
           const fmpz_mpoly_t B, fmpz * maxBfields,
           const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);

FLINT_DLL int _fmpz_mpoly_mul_dense(fmpz_mpoly_t P,
                                 const fmpz_mpoly_t A, fmpz * maxAfields,
                                 const fmpz_mpoly_t B, fmpz * maxBfields,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Powering ******************************************************************/

FLINT_DLL int fmpz_mpoly_pow_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                  const fmpz_t k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_pow_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                          ulong k, const fmpz_mpoly_ctx_t ctx);

/* Division ******************************************************************/

FLINT_DLL int fmpz_mpoly_divides(fmpz_mpoly_t Q,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t Q,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_divides_heap_threaded(fmpz_mpoly_t Q,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_divides_heap_threaded_pool(fmpz_mpoly_t Q,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx,
                        const thread_pool_handle * handles, slong num_handles);

FLINT_DLL slong _fmpz_mpoly_divides_array(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                        const fmpz * poly3, const ulong * exp3, slong len3,
                                         slong * mults, slong num, slong bits);

FLINT_DLL int fmpz_mpoly_divides_array(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, 
                                                   const fmpz_mpoly_ctx_t ctx);


FLINT_DLL int mpoly_divides_select_exps(fmpz_mpoly_t S, fmpz_mpoly_ctx_t zctx,
                                slong nworkers, ulong * Aexp, slong Alen,
                                   ulong * Bexp, slong Blen, flint_bitcnt_t bits);

FLINT_DLL slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1,
                      ulong ** exp1, slong * alloc, const fmpz * poly2,
                    const ulong * exp2, slong len2, const fmpz * poly3,
                    const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N,
                                                        const ulong * cmpmask);

FLINT_DLL void fmpz_mpoly_divrem(fmpz_mpoly_t Q, fmpz_mpoly_t R,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_quasidivrem(fmpz_t scale, fmpz_mpoly_t Q,
     fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_div(fmpz_mpoly_t Q, const fmpz_mpoly_t A,
                             const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_quasidiv(fmpz_t scale, fmpz_mpoly_t Q,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** Q,
     fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B,
                                        slong len, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_quasidivrem_ideal(fmpz_t scale,
     fmpz_mpoly_struct ** Q, fmpz_mpoly_t R, const fmpz_mpoly_t A,
         fmpz_mpoly_struct * const * B, slong len, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_divexact(fmpz_mpoly_t Q, const fmpz_mpoly_t A,
                             const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_divides(Q, A, B, ctx))
        return;

    flint_throw(FLINT_ERROR, "fmpz_mpoly_divexact: nonexact division");
}

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

/* Square root ***************************************************************/

FLINT_DLL slong _fmpz_mpoly_sqrt_heap(fmpz ** polyq, ulong ** expq,
           slong * allocq, const fmpz * poly2, const ulong * exp2, slong len2,
                       flint_bitcnt_t bits, const mpoly_ctx_t mctx, int check);

FLINT_DLL int fmpz_mpoly_sqrt_heap(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, 
                                        const fmpz_mpoly_ctx_t ctx, int check);

FMPZ_MPOLY_INLINE
int fmpz_mpoly_sqrt(fmpz_mpoly_t q, const fmpz_mpoly_t poly2,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_sqrt_heap(q, poly2, ctx, 1);
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_is_square(const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    int res;
    fmpz_mpoly_t q;
    fmpz_mpoly_init(q, ctx);
    res = fmpz_mpoly_sqrt_heap(q, poly2, ctx, 1);
    fmpz_mpoly_clear(q, ctx);
    return res;
}

/* GCD ***********************************************************************/

FLINT_DLL void fmpz_mpoly_term_content(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_content_vars(fmpz_mpoly_t g, const fmpz_mpoly_t A,
                  slong * vars, slong vars_length, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd(fmpz_mpoly_t G,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_cofactors(fmpz_mpoly_t G,
                fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, const fmpz_mpoly_t A,
                             const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_deflation(fmpz * shift, fmpz * stride,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_deflate(fmpz_mpoly_t A, const fmpz_mpoly_t B,
          const fmpz * shift, const fmpz * stride, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_inflate(fmpz_mpoly_t A, const fmpz_mpoly_t B,
          const fmpz * shift, const fmpz * stride, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_hensel(fmpz_mpoly_t G,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_brown(fmpz_mpoly_t G,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_subresultant(fmpz_mpoly_t G,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_zippel(fmpz_mpoly_t G,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_gcd_zippel2(fmpz_mpoly_t G,
       const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

/* Univariates ***************************************************************/

FLINT_DLL void fmpz_mpoly_univar_init(fmpz_mpoly_univar_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_clear(fmpz_mpoly_univar_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_fit_length(fmpz_mpoly_univar_t A,
                                     slong length, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_print_pretty(const fmpz_mpoly_univar_t A,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_univar_assert_canonical(fmpz_mpoly_univar_t A,
                                                   const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_univar_zero(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)
{
    A->length = 0;
}

FLINT_DLL void fmpz_mpoly_univar_set_coeff_ui(fmpz_mpoly_univar_t A,
                    ulong e, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_to_univar(fmpz_mpoly_univar_t A,
                  const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_from_univar(fmpz_mpoly_t A, flint_bitcnt_t Abits,
           const fmpz_mpoly_univar_t B, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_from_univar(fmpz_mpoly_t A,
           const fmpz_mpoly_univar_t B, slong var, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_univar_swap(fmpz_mpoly_univar_t A, fmpz_mpoly_univar_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   fmpz_mpoly_univar_struct t = *A;
   *A = *B;
   *B = t;
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_univar_degree_fits_si(const fmpz_mpoly_univar_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    return A->length == 0 || fmpz_fits_si(A->exps + 0);
}

FMPZ_MPOLY_INLINE
slong fmpz_mpoly_univar_length(const fmpz_mpoly_univar_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    return A->length;
}

FMPZ_MPOLY_INLINE
slong fmpz_mpoly_univar_get_term_exp_si(fmpz_mpoly_univar_t A, slong i,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    return fmpz_get_si(A->exps + i);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_univar_get_term_coeff(fmpz_mpoly_t c,
              const fmpz_mpoly_univar_t A, slong i, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fmpz_mpoly_set(c, A->coeffs + i, ctx);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_univar_swap_term_coeff(fmpz_mpoly_t c,
                    fmpz_mpoly_univar_t A, slong i, const fmpz_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fmpz_mpoly_swap(c, A->coeffs + i, ctx);
}

FLINT_DLL int fmpz_mpoly_univar_pseudo_gcd(fmpz_mpoly_univar_t gx,
                const fmpz_mpoly_univar_t ax, const fmpz_mpoly_univar_t bx,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_univar_resultant(fmpz_mpoly_t d,
                const fmpz_mpoly_univar_t ax, const fmpz_mpoly_univar_t bx,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_univar_discriminant(fmpz_mpoly_t d,
                     const fmpz_mpoly_univar_t fx, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_resultant(fmpz_mpoly_t R, const fmpz_mpoly_t A,
                  const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_discriminant(fmpz_mpoly_t R, const fmpz_mpoly_t A,
                                        slong var, const fmpz_mpoly_ctx_t ctx);

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

FLINT_DLL void mpoly_void_ring_init_fmpz_mpoly_ctx(mpoly_void_ring_t R,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_pow_fps(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                          ulong k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpolyl_lead_coeff(fmpz_mpoly_t c, const fmpz_mpoly_t A,
                                   slong num_vars, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpolyl_content(fmpz_mpoly_t g, const fmpz_mpoly_t A,
                                   slong num_vars, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_to_fmpz_poly_deflate(fmpz_poly_t A,
                         const fmpz_mpoly_t B, slong var, const ulong * Bshift,
                            const ulong * Bstride, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_from_fmpz_poly_inflate(fmpz_mpoly_t A,
       flint_bitcnt_t Abits, const fmpz_poly_t B, slong var, const ulong * Ashift,
                            const ulong * Astride, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_repack_bits(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_repack_bits_inplace(fmpz_mpoly_t A, flint_bitcnt_t Abits,
                                                   const fmpz_mpoly_ctx_t ctx);

typedef struct _fmpz_mpoly_stripe_struct
{
    char * big_mem;
    slong big_mem_alloc;
    slong N;
    flint_bitcnt_t bits;
    const ulong * cmpmask;
    slong * startidx;
    slong * endidx;
    ulong * emin;
    ulong * emax;
    flint_bitcnt_t coeff_bits;
    int upperclosed;
    int flint_small;
} fmpz_mpoly_stripe_struct;

typedef fmpz_mpoly_stripe_struct fmpz_mpoly_stripe_t[1];


/* mpolyd ********************************************************************/

typedef struct
{
    slong nvars;
    slong * perm;
} fmpz_mpolyd_ctx_struct;

typedef fmpz_mpolyd_ctx_struct fmpz_mpolyd_ctx_t[1];


FLINT_DLL void fmpz_mpolyd_init(fmpz_mpolyd_t poly, slong nvars);

FLINT_DLL void fmpz_mpolyd_fit_length(fmpz_mpolyd_t poly, slong len);

FLINT_DLL void fmpz_mpolyd_clear(fmpz_mpolyd_t poly);

/*****************************************************************************/

typedef struct {
    fmpz * powers;
    slong length;
    slong alloc;
    fmpz_t tmp;
} fmpz_pow_cache_t[1];

void fmpz_pow_cache_init(fmpz_pow_cache_t T, const fmpz_t val);

void fmpz_pow_cache_clear(fmpz_pow_cache_t T);

int fmpz_pow_cache_mulpow_ui(fmpz_t a, const fmpz_t b, ulong k,
                                                          fmpz_pow_cache_t T);

int fmpz_pow_cache_mulpow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t k,
                                                          fmpz_pow_cache_t T);

/*****************************************************************************/

FLINT_DLL void fmpz_mpoly_to_mpoly_perm_deflate_threaded_pool(
                                fmpz_mpoly_t A, const fmpz_mpoly_ctx_t lctx,
                            const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx,
               const slong * perm, const ulong * shift, const ulong * stride,
                        const thread_pool_handle * handles, slong num_handles);

FLINT_DLL void fmpz_mpoly_from_mpoly_perm_inflate(
               fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx,
                          const fmpz_mpoly_t B,  const fmpz_mpoly_ctx_t lctx,
                const slong * perm, const ulong * shift, const ulong * stride);

FLINT_DLL void fmpz_mpoly_height(fmpz_t max,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_heights(fmpz_t max, fmpz_t sum,
                             const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

/* geobuckets ****************************************************************/

typedef struct fmpz_mpoly_geobucket
{
    fmpz_mpoly_struct polys[FLINT_BITS/2];
    fmpz_mpoly_struct temps[FLINT_BITS/2];
    slong length;
} fmpz_mpoly_geobucket_struct;

typedef fmpz_mpoly_geobucket_struct fmpz_mpoly_geobucket_t[1];

FLINT_DLL void fmpz_mpoly_geobucket_init(fmpz_mpoly_geobucket_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_clear(fmpz_mpoly_geobucket_t B,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_empty(fmpz_mpoly_t p,
                         fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_fit_length(fmpz_mpoly_geobucket_t B,
                                          slong i, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_set(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_add(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_geobucket_sub(fmpz_mpoly_geobucket_t B,
                                   fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

/* Helpers for array methods *************************************************/

FLINT_DLL void _fmpz_mpoly_mul_array_chunked_DEG(fmpz_mpoly_t P,
                             const fmpz_mpoly_t A, const fmpz_mpoly_t B, 
                                       ulong degb, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_mul_array_chunked_LEX(fmpz_mpoly_t P,
                             const fmpz_mpoly_t A, const fmpz_mpoly_t B, 
                              const ulong * mults, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_addmul_array1_slong1(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_addmul_array1_slong(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_addmul_array1_slong2(ulong * poly1, 
                 const slong * poly2, const ulong * exp2, slong len2,
                          const slong * poly3, const ulong * exp3, slong len3);

FLINT_DLL void _fmpz_mpoly_addmul_array1_fmpz(fmpz * poly1, 
                 const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3);

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

FLINT_DLL slong fmpz_mpoly_append_array_sm1_LEX(fmpz_mpoly_t P,
                        slong Plen, ulong * coeff_array,
                  const ulong * mults, slong num, slong array_size, slong top);
FLINT_DLL slong fmpz_mpoly_append_array_sm2_LEX(fmpz_mpoly_t P,
                        slong Plen, ulong * coeff_array,
                  const ulong * mults, slong num, slong array_size, slong top);
FLINT_DLL slong fmpz_mpoly_append_array_sm3_LEX(fmpz_mpoly_t P,
                         slong Plen, ulong * coeff_array,
                  const ulong * mults, slong num, slong array_size, slong top);
FLINT_DLL slong fmpz_mpoly_append_array_fmpz_LEX(fmpz_mpoly_t P,
                        slong Plen, fmpz * coeff_array,
                  const ulong * mults, slong num, slong array_size, slong top);

FLINT_DLL slong fmpz_mpoly_append_array_sm1_DEGLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_sm2_DEGLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_sm3_DEGLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_fmpz_DEGLEX(fmpz_mpoly_t P,
           slong Plen, fmpz * coeff_array, slong top, slong nvars, slong degb);

FLINT_DLL slong fmpz_mpoly_append_array_sm1_DEGREVLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_sm2_DEGREVLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_sm3_DEGREVLEX(fmpz_mpoly_t P,
          slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
FLINT_DLL slong fmpz_mpoly_append_array_fmpz_DEGREVLEX(fmpz_mpoly_t P,
           slong Plen, fmpz * coeff_array, slong top, slong nvars, slong degb);

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

FLINT_DLL void _fmpz_mpoly_to_ulong_array2(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_ulong_array1(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_ulong_array(ulong * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);

FLINT_DLL void _fmpz_mpoly_to_fmpz_array(fmpz * p, const fmpz * coeffs,
                                                const ulong * exps, slong len);


/* Misc arithmetic - has nothing to do with mpoly, should be moved out *******/

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


FLINT_DLL mpz_srcptr _fmpz_mpoly_get_mpz_signed_uiuiui(ulong * sm, fmpz x,
                                                                    mpz_ptr t);

FMPZ_MPOLY_INLINE
void flint_mpz_add_signed_uiuiui(mpz_ptr a, mpz_srcptr b,
                                                 ulong c2, ulong c1, ulong c0)
{
    ulong cs, d[3];
    mpz_t c;

    c->_mp_d = d;
    c->_mp_alloc = 3;

    cs = FLINT_SIGN_EXT(c2);

    sub_dddmmmsss(d[2], d[1], d[0], cs^c2, cs^c1, cs^c0, cs, cs, cs);

    c->_mp_size = d[2] != 0 ? 3 :
                  d[1] != 0 ? 2 :
                  d[0] != 0;
    if (cs != 0)
        c->_mp_size = -c->_mp_size;

    mpz_add(a, b, c);
}



/******************************************************************************

   Internal consistency checks

******************************************************************************/

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

    if (bits <= FLINT_BITS)
        mask = mpoly_overflow_mask_sp(bits);
    else
        mask = 0;

    for (i = 0; i < r->length; i++)
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides && fmpz_cmpabs(g->coeffs + 0, r->coeffs + i) <= 0)
        {
            flint_printf("fmpz_mpoly_remainder_test FAILED i = %wd\n", i);
            flint_printf("rem ");fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
            flint_printf("den ");fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
            flint_abort();
        }
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

    if (bits <= FLINT_BITS)
        mask = mpoly_overflow_mask_sp(bits);
    else
        mask = 0;

    for (i = 0; i < r->length; i++)
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides)
        {
            flint_printf("fmpz_mpoly_remainder_strongtest FAILED i = %wd\n", i);
            flint_printf("rem ");fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
            flint_printf("den ");fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
            flint_abort();
        }
    }

    flint_free(rexp);
    flint_free(gexp);
}

#ifdef __cplusplus
}
#endif

#endif

