/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_MPOLY_H
#define FQ_NMOD_MPOLY_H

#ifdef FQ_NMOD_MPOLY_INLINES_C
#define FQ_NMOD_MPOLY_INLINE FLINT_DLL
#else
#define FQ_NMOD_MPOLY_INLINE static __inline__
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

#include "nmod_mpoly.h"


#ifdef __cplusplus
 extern "C" {
#endif


/*  Type definitions *********************************************************/

/*
    context object for fq_nmod_mpoly
*/
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
    flint_bitcnt_t bits;     /* number of bits per exponent */
} fq_nmod_mpoly_struct;

typedef fq_nmod_mpoly_struct fq_nmod_mpoly_t[1];

FQ_NMOD_MPOLY_INLINE
fq_nmod_struct * fq_nmod_mpoly_term_coeff_ref(fq_nmod_mpoly_t A, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < A->length);
    return A->coeffs + i;
}


/* Internal type definitions *************************************************/

/*
    fq_nmod_mpoly_univar_t
    sparse univariates with multivariate coefficients
*/
typedef struct
{
   fq_nmod_mpoly_struct * coeffs; /* multivariate coefficients */
   fmpz * exps;
   slong alloc;
   slong length;
} fq_nmod_mpoly_univar_struct;

typedef fq_nmod_mpoly_univar_struct fq_nmod_mpoly_univar_t[1];

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
   flint_bitcnt_t bits;    /* default bits to construct coeffs */
} fq_nmod_mpolyu_struct;

typedef fq_nmod_mpolyu_struct fq_nmod_mpolyu_t[1];

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
    flint_bitcnt_t bits;   /* default bits to construct coeffs */
} fq_nmod_mpolyun_struct;

typedef fq_nmod_mpolyun_struct fq_nmod_mpolyun_t[1];

/*
    fq_nmod_mpolyd_t
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
    slong coeff_alloc;
    fq_nmod_struct * coeffs;
} fq_nmod_mpolyd_struct;

typedef fq_nmod_mpolyd_struct fq_nmod_mpolyd_t[1];

/*
    fq_nmod_mpoly_geobucket_t
    power of 4 increment
*/
typedef struct fq_nmod_mpoly_geobucket
{
    fq_nmod_mpoly_struct polys[FLINT_BITS/2];
    slong length;
} fq_nmod_mpoly_geobucket_struct;

typedef fq_nmod_mpoly_geobucket_struct fq_nmod_mpoly_geobucket_t[1];


/* Embeddings ****************************************************************/

/* see fq_nmod_mpoly/fq_nmod_embed.c for more info */

typedef struct bad_fq_nmod_embed
{
    const fq_nmod_ctx_struct * smctx; /* modulus is f */
    fq_nmod_poly_t phi_sm;      /* phi as an element of F_p[theta][x] */
    fq_nmod_poly_t h;
    const fq_nmod_ctx_struct * lgctx; /* modulus is g */
    fq_nmod_t theta_lg;         /* theta as an element of F_p[phi]/g(phi) */
    fq_nmod_t x_lg;             /* x as an element of F_p[phi]/g(phi) */
} bad_fq_nmod_embed_struct;

typedef bad_fq_nmod_embed_struct bad_fq_nmod_embed_t[1];


FLINT_DLL void bad_fq_nmod_embed_clear(bad_fq_nmod_embed_t emb);

FLINT_DLL void bad_fq_nmod_embed_array_clear(bad_fq_nmod_embed_struct * emb, slong m);

FLINT_DLL void bad_fq_nmod_embed_array_init(
    bad_fq_nmod_embed_struct * emb,
    const fq_nmod_ctx_t bigctx, /* F_p[phi]/g(phi) */
    const fq_nmod_ctx_t smallctx);

FLINT_DLL void bad_fq_nmod_embed_sm_to_lg(
    fq_nmod_t out,            /* element of lgctx */
    const fq_nmod_poly_t in,  /* poly over smctx */
    const bad_fq_nmod_embed_t emb);

FLINT_DLL void bad_fq_nmod_embed_lg_to_sm(
    fq_nmod_poly_t out,  /* poly over smctx */
    const fq_nmod_t in,  /* element of lgctx */
    const bad_fq_nmod_embed_t emb);


typedef struct bad_fq_nmod_mpoly_embed_chooser
{
    bad_fq_nmod_embed_struct * embed;
    slong m; /* degree of the extension F_q / F_p */
    slong n; /* degree of the extension F_q^n / F_q */
    slong k; /* index of current in embed */
    mp_limb_t p;
} bad_fq_nmod_mpoly_embed_chooser_struct;

typedef bad_fq_nmod_mpoly_embed_chooser_struct bad_fq_nmod_mpoly_embed_chooser_t[1];

FLINT_DLL bad_fq_nmod_embed_struct * bad_fq_nmod_mpoly_embed_chooser_init(
    bad_fq_nmod_mpoly_embed_chooser_t embc,
    fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate);

FLINT_DLL void bad_fq_nmod_mpoly_embed_chooser_clear(
    bad_fq_nmod_mpoly_embed_chooser_t embc,
    fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate);

FLINT_DLL bad_fq_nmod_embed_struct * bad_fq_nmod_mpoly_embed_chooser_next(
    bad_fq_nmod_mpoly_embed_chooser_t embc,
    fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t randstate);


/* Context object ************************************************************/

FLINT_DLL void fq_nmod_mpoly_ctx_init_deg(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                                 const ordering_t ord, mp_limb_t p, slong deg);

FLINT_DLL void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                              const ordering_t ord, const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpoly_ctx_init_rand(fq_nmod_mpoly_ctx_t ctx,
                                       flint_rand_t state, slong max_nvars,
                                          flint_bitcnt_t p_bits, slong deg_bound);

FLINT_DLL void fq_nmod_mpoly_ctx_clear(fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_ctx_nvars(const fq_nmod_mpoly_ctx_t ctx)
{
    return ctx->minfo->nvars;
}

FQ_NMOD_MPOLY_INLINE
ordering_t fq_nmod_mpoly_ctx_ord(const fq_nmod_mpoly_ctx_t ctx)
{
    return ctx->minfo->ord;
}


/*  Memory management ********************************************************/

FLINT_DLL void fq_nmod_mpoly_init(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_init2(fq_nmod_mpoly_t A, slong alloc,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_init3(fq_nmod_mpoly_t A, slong alloc,
                              flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_realloc(fq_nmod_mpoly_t A,
                                   slong alloc, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_fit_length(fq_nmod_struct ** coeff,
                              ulong ** exps, slong * alloc, slong len, slong N,
                                                    const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpoly_clear(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);


FQ_NMOD_MPOLY_INLINE
void _fq_nmod_mpoly_set_length(fq_nmod_mpoly_t A, slong newlen, 
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(newlen <= A->alloc);
    A->length = newlen;
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_truncate(fq_nmod_mpoly_t A, slong newlen, 
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    if (A->length > newlen)
    {
        A->length = newlen;
    }
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_fit_bits(fq_nmod_mpoly_t A, slong bits,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    if (A->bits < bits)
    {
        if (A->alloc != 0)
        {
            slong N = mpoly_words_per_exp(bits, ctx->minfo);
            ulong * t = (ulong *) flint_malloc(N*A->alloc*sizeof(ulong));
            mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length,
                                                                   ctx->minfo);
            flint_free(A->exps);
            A->exps = t;
        }

        A->bits = bits;
    }
}


/* Input/output **************************************************************/

FLINT_DLL int fq_nmod_mpoly_set_str_pretty(fq_nmod_mpoly_t A, const char * str,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL char * fq_nmod_mpoly_get_str_pretty(const fq_nmod_mpoly_t A,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_fprint_pretty(FILE * file, 
      const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A,
                                const char ** x, const fq_nmod_mpoly_ctx_t ctx)
{
   return fq_nmod_mpoly_fprint_pretty(stdout, A, x, ctx);
}


/*  Basic manipulation *******************************************************/

FLINT_DLL void fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, slong var,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_is_gen(const fq_nmod_mpoly_t A,
                                     slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_equal(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_swap(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
   fq_nmod_mpoly_struct t = *A;
   *A = *B;
   *B = t;
}


/* Constants *****************************************************************/

FLINT_DLL int fq_nmod_mpoly_is_fq_nmod(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_fq_nmod(fq_nmod_mpoly_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_ui(fq_nmod_mpoly_t A, ulong c,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_fq_nmod_gen(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_equal_fq_nmod(const fq_nmod_mpoly_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_zero(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
   _fq_nmod_mpoly_set_length(A, 0, ctx);
}

FLINT_DLL void fq_nmod_mpoly_one(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_is_zero(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
   return A->length == 0;
}

FLINT_DLL int fq_nmod_mpoly_is_one(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* Degrees *******************************************************************/

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_degrees_fit_si(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
               : mpoly_degrees_fit_si(A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_degrees_fmpz(fmpz ** degs, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_degrees_pfmpz(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_degrees_si(slong * degs, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_degrees_si(degs, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_degree_fmpz(fmpz_t deg, const fq_nmod_mpoly_t A, slong var,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_degree_fmpz(deg, A->exps, A->length, A->bits, var, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_degree_si(const fq_nmod_mpoly_t A, slong var,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return mpoly_degree_si(A->exps, A->length, A->bits, var, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_total_degree_fits_si(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_fits_si(A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_total_degree_fmpz(fmpz_t td, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    mpoly_total_degree_fmpz(td, A->exps, A->length, A->bits, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_total_degree_si(const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return mpoly_total_degree_si(A->exps, A->length, A->bits, ctx->minfo);
}


/* Coefficients **************************************************************/

FLINT_DLL void fq_nmod_mpoly_get_coeff_fq_nmod_monomial(fq_nmod_t c,
                          const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t M,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_coeff_fq_nmod_monomial(fq_nmod_mpoly_t A,
                                  const fq_nmod_t c, const fq_nmod_mpoly_t M,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(fq_nmod_t c,
                             const fq_nmod_mpoly_t A, fmpz * const * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_coeff_fq_nmod_ui(fq_nmod_t c,
                                const fq_nmod_mpoly_t A, const ulong * exp,
                                               const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A,
                                       const fq_nmod_t c, const fmpz * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A,
                                     const fq_nmod_t c, fmpz * const * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_coeff_fq_nmod_ui(fq_nmod_mpoly_t A,
                                      const fq_nmod_t c, const ulong * exp,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_coeff_vars_ui(fq_nmod_mpoly_t C,
           const fq_nmod_mpoly_t A, slong * vars, ulong * exps, slong length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE fq_nmod_struct * fq_nmod_mpoly_leadcoeff(
                        const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}


/* comparison ****************************************************************/

FLINT_DLL int fq_nmod_mpoly_cmp(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* container operations ******************************************************/

FLINT_DLL int fq_nmod_mpoly_is_canonical(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_length(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length;
}

FLINT_DLL void fq_nmod_mpoly_resize(fq_nmod_mpoly_t A, slong new_length,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term_coeff_fq_nmod(fq_nmod_t c,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_term_coeff_fq_nmod(fq_nmod_mpoly_t A, 
                    slong i, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_term_exp_fits_ui(const fq_nmod_mpoly_t A, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_ui(A->exps, A->bits, i, ctx->minfo);
}

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_term_exp_fits_si(const fq_nmod_mpoly_t A, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_si(A->exps, A->bits, i, ctx->minfo);
}

FLINT_DLL void fq_nmod_mpoly_get_term_exp_fmpz(fmpz ** exp,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term_exp_ui(ulong * exp,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term_exp_si(slong * exp,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL ulong fq_nmod_mpoly_get_term_var_exp_ui(const fq_nmod_mpoly_t A,
                            slong i, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong fq_nmod_mpoly_get_term_var_exp_si(const fq_nmod_mpoly_t A,
                            slong i, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_term_exp_fmpz(fq_nmod_mpoly_t A,
                   slong i, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_set_term_exp_ui(fq_nmod_mpoly_t A,
                    slong i, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A,
                                       slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_get_term_monomial(fq_nmod_mpoly_t M,
              const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_push_term_fq_nmod_fmpz(fq_nmod_mpoly_t A,
         const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_push_term_fq_nmod_ui(fq_nmod_mpoly_t A,
          const fq_nmod_t c, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_sort_terms(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_combine_like_terms(fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_reverse(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_assert_canonical(const fq_nmod_mpoly_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_radix_sort1(fq_nmod_mpoly_t A, slong left,
                 slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask);

FLINT_DLL void _fq_nmod_mpoly_radix_sort(fq_nmod_mpoly_t A, slong left,
                       slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask);

FLINT_DLL void _fq_nmod_mpoly_push_exp_ffmpz(fq_nmod_mpoly_t A,
                              const fmpz * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_push_exp_pfmpz(fq_nmod_mpoly_t A,
                            fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_push_exp_ui(fq_nmod_mpoly_t A,
                             const ulong * exp, const fq_nmod_mpoly_ctx_t ctx);


/* Random generation *********************************************************/

FLINT_DLL void fq_nmod_mpoly_randtest_bound(fq_nmod_mpoly_t A, flint_rand_t state,
                 slong length, ulong exp_bound, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_randtest_bounds(fq_nmod_mpoly_t A, flint_rand_t state,
              slong length, ulong * exp_bounds, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_randtest_bits(fq_nmod_mpoly_t A, flint_rand_t state,
            slong length, flint_bitcnt_t exp_bits, const fq_nmod_mpoly_ctx_t ctx);


/* Addition/Subtraction ******************************************************/

FLINT_DLL slong _fq_nmod_mpoly_add(
                         fq_nmod_struct * coeff1,       ulong * exp1,
                         fq_nmod_struct * coeff2, const ulong * exp2, slong len2,
                         fq_nmod_struct * coeff3, const ulong * exp3, slong len3,
                   slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpoly_add_fq_nmod(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_sub_fq_nmod(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_add(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_sub(fq_nmod_mpoly_t A,
                            const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C,
                                                const fq_nmod_mpoly_ctx_t ctx);


/* Scalar operations *********************************************************/

FLINT_DLL void fq_nmod_mpoly_neg(fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_scalar_mul_fq_nmod(fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_make_monic(fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);


/* Differentiation **********************************************************/

FLINT_DLL void fq_nmod_mpoly_derivative(fq_nmod_mpoly_t A,
            const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);


/* Multiplication ************************************************************/

FLINT_DLL void fq_nmod_mpoly_mul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                       const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_mul_johnson(fq_nmod_mpoly_t poly1,
                    const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong _fq_nmod_mpoly_mul_johnson(
                    fq_nmod_struct ** coeff1, ulong ** exp1, slong * alloc,
             const fq_nmod_struct * coeff2, const ulong * exp2, slong len2,
             const fq_nmod_struct * coeff3, const ulong * exp3, slong len3,
  flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);


/* Powering ******************************************************************/

FLINT_DLL int fq_nmod_mpoly_pow_fmpz(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_pow_ui(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                       ulong k, const fq_nmod_mpoly_ctx_t ctx);


/* Division ******************************************************************/

FLINT_DLL int fq_nmod_mpoly_divides(fq_nmod_mpoly_t Q,
                         const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_div(fq_nmod_mpoly_t Q,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_divrem(fq_nmod_mpoly_t Q, fq_nmod_mpoly_t R,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_divrem_ideal(fq_nmod_mpoly_struct ** Q,
                                  fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A,
                                  fq_nmod_mpoly_struct * const * B, slong len,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_divides_monagan_pearce(fq_nmod_mpoly_t poly1,
                  const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_div_monagan_pearce(fq_nmod_mpoly_t q,
                      const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_divrem_monagan_pearce(fq_nmod_mpoly_t q, fq_nmod_mpoly_t r,
                      const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_divrem_ideal_monagan_pearce(
                        fq_nmod_mpoly_struct ** q, fq_nmod_mpoly_t r,
            const fq_nmod_mpoly_t poly2, fq_nmod_mpoly_struct * const * poly3,
                                      slong len, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong _fq_nmod_mpoly_divides_monagan_pearce(
                  fq_nmod_struct ** coeff1,      ulong ** exp1, slong * alloc,
             const fq_nmod_struct * coeff2, const ulong * exp2, slong len2,
             const fq_nmod_struct * coeff3, const ulong * exp3, slong len3,
  flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx);


/* GCD ***********************************************************************/

FLINT_DLL int fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int _fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_gcd_cofactors(fq_nmod_mpoly_t G,
          fq_nmod_mpoly_t Abar, fq_nmod_mpoly_t Bbar, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int _fq_nmod_mpoly_gcd_cofactors(
                                fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
                                fq_nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
                                fq_nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
                             const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_gcd_brown(fq_nmod_mpoly_t G,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_gcd_zippel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_deflation(fmpz * shift, fmpz * stride,
                       const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_deflate(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
          const fmpz * shift, const fmpz * stride, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_inflate(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
       const fmpz * shift, const fmpz * stride, const fq_nmod_mpoly_ctx_t ctx);


/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

FLINT_DLL void fq_nmod_mpoly_pow_rmul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                       ulong k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_to_fq_nmod_poly_deflate(fq_nmod_poly_t A,
                   const fq_nmod_mpoly_t B, slong var, const ulong * Bshift,
                         const ulong * Bstride, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void _fq_nmod_mpoly_from_fq_nmod_poly_inflate(fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits, const fq_nmod_poly_t B, slong var, const ulong * Ashift,
                         const ulong * Astride, const fq_nmod_mpoly_ctx_t ctx);


FLINT_DLL int fq_nmod_mpoly_repack_bits(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                          flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_repack_bits_inplace(fq_nmod_mpoly_t A,
                         flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx);


FLINT_DLL void fq_nmod_mpoly_ctx_change_modulus(fq_nmod_mpoly_ctx_t ctx,
                                                                    slong deg);


/* Univariates ***************************************************************/

FLINT_DLL void fq_nmod_mpoly_univar_init(fq_nmod_mpoly_univar_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_univar_clear(fq_nmod_mpoly_univar_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_univar_fit_length(fq_nmod_mpoly_univar_t A,
                                  slong length, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_univar_print_pretty(const fq_nmod_mpoly_univar_t A,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_univar_assert_canonical(fq_nmod_mpoly_univar_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_to_univar(fq_nmod_mpoly_univar_t A,
            const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_from_univar_bits(fq_nmod_mpoly_t A, flint_bitcnt_t Abits,
     const fq_nmod_mpoly_univar_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_from_univar(fq_nmod_mpoly_t A,
     const fq_nmod_mpoly_univar_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_univar_swap(fq_nmod_mpoly_univar_t A,
                       fq_nmod_mpoly_univar_t B, const fq_nmod_mpoly_ctx_t ctx)
{
   fq_nmod_mpoly_univar_struct t = *A;
   *A = *B;
   *B = t;
}

FQ_NMOD_MPOLY_INLINE
int fq_nmod_mpoly_univar_degree_fits_si(const fq_nmod_mpoly_univar_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length == 0 || fmpz_fits_si(A->exps + 0);
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_univar_length(const fq_nmod_mpoly_univar_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length;
}

FQ_NMOD_MPOLY_INLINE
slong fq_nmod_mpoly_univar_get_term_exp_si(fq_nmod_mpoly_univar_t A, slong i,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    return fmpz_get_si(A->exps + i);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_univar_get_term_coeff(fq_nmod_mpoly_t c,
        const fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fq_nmod_mpoly_set(c, A->coeffs + i, ctx);
}

FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_univar_swap_term_coeff(fq_nmod_mpoly_t c,
              fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT((ulong)i < (ulong)A->length);
    fq_nmod_mpoly_swap(c, A->coeffs + i, ctx);
}


/* mpolyd ********************************************************************/

typedef struct
{
    slong nvars;
    slong * perm;
    fq_nmod_ctx_t fqctx;
} fq_nmod_mpolyd_ctx_struct;
typedef fq_nmod_mpolyd_ctx_struct fq_nmod_mpolyd_ctx_t[1];

FLINT_DLL void fq_nmod_mpolyd_ctx_init(fq_nmod_mpolyd_ctx_t dctx, slong nvars,
                                                       mp_limb_t p, slong deg);

FLINT_DLL void fq_nmod_mpolyd_ctx_init2(fq_nmod_mpolyd_ctx_t dctx, slong nvars,
                                                    const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpolyd_ctx_clear(fq_nmod_mpolyd_ctx_t dctx);

FLINT_DLL int fq_nmod_mpolyd_ctx_set_for_gcd(fq_nmod_mpolyd_ctx_t dctx,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyd_init(fq_nmod_mpolyd_t poly, slong nvars,
                                                    const fq_nmod_ctx_t fqctx);

FQ_NMOD_MPOLY_INLINE void fq_nmod_mpolyd_swap(fq_nmod_mpolyd_t poly1,
                                                        fq_nmod_mpolyd_t poly2)
{
   fq_nmod_mpolyd_struct t = *poly1;
   *poly1 = *poly2;
   *poly2 = t;
}

FLINT_DLL void fq_nmod_mpolyd_fit_length(fq_nmod_mpolyd_t poly, slong len,
                                                    const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpolyd_set_nvars(fq_nmod_mpolyd_t poly, slong nvars);

FLINT_DLL void fq_nmod_mpolyd_zero(fq_nmod_mpolyd_t poly,
                                                    const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpolyd_clear(fq_nmod_mpolyd_t poly,
                                                    const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpolyd_print(fq_nmod_mpolyd_t poly,
                                                    const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_mpoly_convert_to_fq_nmod_mpolyd(
                         fq_nmod_mpolyd_t A, const fq_nmod_mpolyd_ctx_t dctx,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx);

/* mpolyu ********************************************************************/

FLINT_DLL int fq_nmod_mpolyu_is_canonical(const fq_nmod_mpolyu_t poly,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_init(fq_nmod_mpolyu_t A, flint_bitcnt_t bits,
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

FLINT_DLL void fq_nmod_mpolyu_one(fq_nmod_mpolyu_t A,
                                               const fq_nmod_mpoly_ctx_t uctx);

FLINT_DLL void fq_nmod_mpolyu_shift_right(fq_nmod_mpolyu_t A, ulong s);

FLINT_DLL void fq_nmod_mpolyu_shift_left(fq_nmod_mpolyu_t A, ulong s);

FLINT_DLL int fq_nmod_mpolyu_content_mpoly(fq_nmod_mpoly_t g,
                      const fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_scalar_mul_fq_nmod(fq_nmod_mpolyu_t A,
                                   fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_set(fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B,
                                               const fq_nmod_mpoly_ctx_t uctx);

FLINT_DLL void fq_nmod_mpolyu_cvtto_poly(fq_nmod_poly_t a, fq_nmod_mpolyu_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtfrom_poly(fq_nmod_mpolyu_t A, fq_nmod_poly_t a,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_cvtfrom_poly_notmain(fq_nmod_mpoly_t A,
                   fq_nmod_poly_t a, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtfrom_poly_notmain(fq_nmod_mpolyu_t A,
                   fq_nmod_poly_t a, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_setform(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_to_mpolyu_perm_deflate(
                         fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx,
                 const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx,
                const slong * perm, const ulong * shift, const ulong * stride);

FLINT_DLL void fq_nmod_mpoly_from_mpolyu_perm_inflate(
        fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx,
          const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t uctx,
                const slong * perm, const ulong * shift, const ulong * stride);

FLINT_DLL int fq_nmod_mpolyuu_divides(fq_nmod_mpolyu_t Q,
          const fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B, slong nmainvars,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_divexact_mpoly(fq_nmod_mpolyu_t A,
                                      fq_nmod_mpolyu_t B, fq_nmod_mpoly_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_divexact_mpoly_inplace(fq_nmod_mpolyu_t A,
                             fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_mul_mpoly(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                             fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_mul_mpoly_inplace(fq_nmod_mpolyu_t A,
                             fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyu_gcdm_zippel(fq_nmod_mpolyu_t G,
             fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar,  fq_nmod_mpolyu_t A,
            fq_nmod_mpolyu_t B, fq_nmod_mpoly_ctx_t ctx, mpoly_zipinfo_t zinfo,
                                                       flint_rand_t randstate);

FQ_NMOD_MPOLY_INLINE fq_nmod_struct * fq_nmod_mpolyu_leadcoeff(
                       const fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return fq_nmod_mpoly_leadcoeff(A->coeffs + 0, ctx);
}


/* mpolyn ********************************************************************/

FLINT_DLL void fq_nmod_mpolyn_init(fq_nmod_mpolyn_t A, flint_bitcnt_t bits,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_clear(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_swap(fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B);

FLINT_DLL int fq_nmod_mpolyn_is_canonical(const fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_zero(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_one(fq_nmod_mpolyn_t A,
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

FQ_NMOD_MPOLY_INLINE fq_nmod_struct * fq_nmod_mpolyn_leadcoeff(
                             fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_poly_struct * leadpoly;
    FLINT_ASSERT(A->length > 0);
    leadpoly = A->coeffs + 0;
    FLINT_ASSERT(leadpoly->length > 0);
    return leadpoly->coeffs + leadpoly->length - 1;
}

FQ_NMOD_MPOLY_INLINE fq_nmod_poly_struct * fq_nmod_mpolyn_leadcoeff_poly(
                       const fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return A->coeffs + 0;
}

/* mpolyun *******************************************************************/

FLINT_DLL void fq_nmod_mpoly_to_mpolyn_perm_deflate(fq_nmod_mpolyn_t A,
                const fq_nmod_mpoly_ctx_t nctx,
                const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx,
                const slong * perm, const ulong * shift, const ulong * stride);

FLINT_DLL void fq_nmod_mpoly_from_mpolyn_perm_inflate(fq_nmod_mpoly_t A,
                flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx,
                const fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t nctx,
                const slong * perm, const ulong * shift, const ulong * stride);

FLINT_DLL void fq_nmod_mpolyun_init(fq_nmod_mpolyun_t A,
                           flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_clear(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyun_is_canonical(const fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_print_pretty(const fq_nmod_mpolyun_t poly,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_swap(fq_nmod_mpolyun_t A, fq_nmod_mpolyun_t B);

FLINT_DLL void fq_nmod_mpolyun_zero(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_fit_length(fq_nmod_mpolyun_t A,
                                  slong length, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_one(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyn_is_nonzero_fq_nmod(const fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyun_is_nonzero_fq_nmod(const fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_scalar_mul_fq_nmod(fq_nmod_mpolyn_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_scalar_mul_fq_nmod(fq_nmod_mpolyun_t A,
                             const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_shift_right(fq_nmod_mpolyun_t A, ulong s);

FLINT_DLL void fq_nmod_mpolyun_shift_left(fq_nmod_mpolyun_t A, ulong s);

FLINT_DLL void fq_nmod_mpolyn_mul_poly(fq_nmod_mpolyn_t A,
                        const fq_nmod_mpolyn_t B, const fq_nmod_poly_t c,
                              const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t t);

FLINT_DLL void fq_nmod_mpolyun_mul_poly(fq_nmod_mpolyun_t A,
                         const fq_nmod_mpolyun_t B, const fq_nmod_poly_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_divexact_poly(fq_nmod_mpolyn_t A,
                      const fq_nmod_mpolyn_t B, const fq_nmod_poly_t c,
            const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t q, fq_nmod_poly_t r);

FLINT_DLL void fq_nmod_mpolyun_divexact_poly(fq_nmod_mpolyun_t A,
                        const fq_nmod_mpolyun_t B, const fq_nmod_poly_t c,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_content_poly(fq_nmod_poly_t g,
                            fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_content_poly(fq_nmod_poly_t g,
                           fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong fq_nmod_mpolyn_lastdeg(fq_nmod_mpolyn_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL slong fq_nmod_mpolyun_lastdeg(fq_nmod_mpolyun_t A,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_set(
    fq_nmod_mpolyun_t A,
    const fq_nmod_mpolyun_t B,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_cvtto_mpolyn(
    fq_nmod_mpolyn_t A,
    fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtto_mpolyun(
    fq_nmod_mpolyun_t A,
    fq_nmod_mpolyu_t B,
    slong k,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_cvtfrom_mpolyn(
    fq_nmod_mpoly_t A,
    fq_nmod_mpolyn_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyu_cvtfrom_mpolyun(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyun_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx);

FQ_NMOD_MPOLY_INLINE fq_nmod_poly_struct * fq_nmod_mpolyun_leadcoeff_poly(
                      const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return fq_nmod_mpolyn_leadcoeff_poly(A->coeffs + 0, ctx);
}


/*****************************************************************************/

FLINT_DLL int fq_nmod_next(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx);

FLINT_DLL nmod_gcds_ret_t fq_nmod_mpolyu_gcds_zippel(fq_nmod_mpolyu_t G,
                fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, fq_nmod_mpolyu_t f,
                                    slong var, const fq_nmod_mpoly_ctx_t ctx,
                                     flint_rand_t randstate, slong * degbound);

FLINT_DLL int fq_nmod_mpolyu_gcdp_zippel_univar(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t Abar,
    fq_nmod_mpolyu_t Bbar,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyu_gcdp_zippel_univar_no_cofactors(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyu_gcdp_zippel(fq_nmod_mpolyu_t G,
              fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar, fq_nmod_mpolyu_t A,
                  fq_nmod_mpolyu_t B, slong var, const fq_nmod_mpoly_ctx_t ctx,
                                mpoly_zipinfo_t zinfo, flint_rand_t randstate);

FLINT_DLL int fq_nmod_mpolyn_gcd_brown_smprime(fq_nmod_mpolyn_t G,
         fq_nmod_mpolyn_t Abar, fq_nmod_mpolyn_t Bbar, fq_nmod_mpolyn_t A,
                 fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyn_gcd_brown_lgprime(fq_nmod_mpolyn_t G,
         fq_nmod_mpolyn_t Abar, fq_nmod_mpolyn_t Bbar, fq_nmod_mpolyn_t A,
                 fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

/* gcd_helper_eval_interp ****************************************************/

FLINT_DLL void nmod_mpolyn_interp_reduce_lg_poly(fq_nmod_poly_t E,
                            const fq_nmod_ctx_t fqctx, const nmod_mpolyn_t A,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_interp_lift_lg_poly(slong * lastdeg_,
                         nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx,
                            const fq_nmod_poly_t B, const fq_nmod_ctx_t fqctx);

FLINT_DLL int nmod_mpolyn_interp_crt_lg_poly(slong * lastdeg_,
     nmod_mpolyn_t F, nmod_mpolyn_t T, nmod_poly_t modulus,
     const nmod_mpoly_ctx_t ctx, fq_nmod_poly_t A,  const fq_nmod_ctx_t fqctx);

FLINT_DLL void nmod_mpolyn_interp_reduce_lg_mpolyn(fq_nmod_mpolyn_t E,
                      fq_nmod_mpoly_ctx_t ectx, nmod_mpolyn_t A, slong var,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyn_interp_lift_lg_mpolyn(slong * lastdeg,
                nmod_mpolyn_t A, slong var, const nmod_mpoly_ctx_t ctx,
                           fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ectx);

FLINT_DLL int nmod_mpolyn_interp_crt_lg_mpolyn(slong * lastdeg_,
            nmod_mpolyn_t F, nmod_mpolyn_t T, nmod_poly_t modulus, slong var,
            const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t A,
                                               const fq_nmod_mpoly_ctx_t ectx);

FLINT_DLL void nmod_mpolyun_interp_reduce_lg_mpolyu(fq_nmod_mpolyu_t A,
                          nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ffctx,
                                                   const nmod_mpoly_ctx_t ctx);

FLINT_DLL void nmod_mpolyun_interp_lift_lg_mpolyu(nmod_mpolyun_t A,
                            const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t Ap,
                                               const fq_nmod_mpoly_ctx_t ctxp);

FLINT_DLL int nmod_mpolyun_interp_crt_lg_mpolyu(slong * lastdeg,
                          nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_poly_t m,
                          const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t A,
                                              const fq_nmod_mpoly_ctx_t ffctx);

FLINT_DLL int nmod_mpolyun_interp_mcrt_lg_mpolyu(slong * lastdeg,
               nmod_mpolyun_t H, const nmod_mpoly_ctx_t ctx, nmod_poly_t m,
                           fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctxp);

FLINT_DLL void fq_nmod_mpolyn_interp_reduce_sm_poly(fq_nmod_poly_t E,
                            const fq_nmod_mpolyn_t A, const fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_interp_lift_sm_poly(fq_nmod_mpolyn_t A,
                        const fq_nmod_poly_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyn_interp_crt_sm_poly(slong * lastdeg_,
            fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, const fq_nmod_poly_t A,
            const fq_nmod_poly_t modulus, const fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_interp_reduce_sm_mpolyn(fq_nmod_mpolyn_t E,
                             fq_nmod_mpolyn_t A, slong var, fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_interp_lift_sm_mpolyn(fq_nmod_mpolyn_t A,
                 fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyn_interp_crt_sm_mpolyn(slong * lastdeg_,
        fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_mpolyn_t A,  slong var,
        fq_nmod_poly_t modulus, const fq_nmod_t alpha,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyn_interp_reduce_lg_poly(fq_nmod_poly_t E,
                        const fq_nmod_mpoly_ctx_t ectx, fq_nmod_mpolyn_t A,
                    const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb);

FLINT_DLL void fq_nmod_mpolyn_interp_lift_lg_poly(slong * lastdeg_,
      fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t B,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

FLINT_DLL int fq_nmod_mpolyn_interp_crt_lg_poly(slong * lastdeg_,
            fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_poly_t modulus,
            const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t A,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

FLINT_DLL void fq_nmod_mpolyn_interp_reduce_lg_mpolyn(fq_nmod_mpolyn_t E,
               const fq_nmod_mpoly_ctx_t ectx, fq_nmod_mpolyn_t A, slong var,
                    const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb);

FLINT_DLL void fq_nmod_mpolyn_interp_lift_lg_mpolyn(slong * lastdeg_,
                         fq_nmod_mpolyn_t A, slong var,
                         const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t B,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

FLINT_DLL int fq_nmod_mpolyn_interp_crt_lg_mpolyn(slong * lastdeg_,
    fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_poly_t modulus, slong var,
    const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t A,
                  const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

FLINT_DLL void fq_nmod_mpolyun_interp_reduce_sm_mpolyu(fq_nmod_mpolyu_t B,
          fq_nmod_mpolyun_t A, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_interp_lift_sm_mpolyu(fq_nmod_mpolyun_t A,
                      const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpolyun_interp_crt_sm_mpolyu(slong * lastdeg,
            fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_mpolyu_t A,
       fq_nmod_poly_t modulus, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpolyun_interp_reduce_lg_mpolyu(fq_nmod_mpolyu_t A,
                        fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ectx,
                    const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb);

FLINT_DLL void fq_nmod_mpolyun_interp_lift_lg_mpolyu(fq_nmod_mpolyun_t A,
                          const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t B,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

FLINT_DLL int fq_nmod_mpolyun_interp_crt_lg_mpolyu(slong * lastdeg,
                fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_poly_t m,
                const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t A,
                   const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb);

FLINT_DLL int fq_nmod_mpolyun_interp_mcrt_lg_mpolyu(slong * lastdeg,
                        fq_nmod_mpolyun_t H, const fq_nmod_mpoly_ctx_t ctx,
                        fq_nmod_poly_t m, fq_nmod_mpolyu_t A,
                         const fq_nmod_mpoly_ctx_t ectx, bad_fq_nmod_embed_t emb);

/* geobuckets ****************************************************************/

FLINT_DLL void fq_nmod_mpoly_geobucket_init(fq_nmod_mpoly_geobucket_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_clear(fq_nmod_mpoly_geobucket_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_empty(fq_nmod_mpoly_t p,
                   fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_print(fq_nmod_mpoly_geobucket_t B,
                               const char ** x, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_fit_length(fq_nmod_mpoly_geobucket_t B,
                                       slong i, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_set(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_add(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_sub(fq_nmod_mpoly_geobucket_t B,
                             fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_set_ui(fq_nmod_mpoly_geobucket_t B,
                                       ulong c, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_set_fq_nmod_gen(fq_nmod_mpoly_geobucket_t B,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_gen(fq_nmod_mpoly_geobucket_t B, slong var,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_add_inplace(fq_nmod_mpoly_geobucket_t B1,
                  fq_nmod_mpoly_geobucket_t B2, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_sub_inplace(fq_nmod_mpoly_geobucket_t B1,
                  fq_nmod_mpoly_geobucket_t B2, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_neg_inplace(fq_nmod_mpoly_geobucket_t B1,
                                                const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_mul_inplace(fq_nmod_mpoly_geobucket_t B1,
                  fq_nmod_mpoly_geobucket_t B2, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_pow_ui_inplace(fq_nmod_mpoly_geobucket_t B1,
                                       ulong k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL void fq_nmod_mpoly_geobucket_pow_fmpz_inplace(fq_nmod_mpoly_geobucket_t B1,
                                const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx);

FLINT_DLL int fq_nmod_mpoly_geobucket_divides_inplace(fq_nmod_mpoly_geobucket_t B1,
                  fq_nmod_mpoly_geobucket_t B2, const fq_nmod_mpoly_ctx_t ctx);


/******************************************************************************

   Internal consistency checks

******************************************************************************/

/*
   test that r is a valid remainder upon division by g
   this means that no monomial of r is divisible by lm(g)
*/
FQ_NMOD_MPOLY_INLINE
void fq_nmod_mpoly_remainder_strongtest(const fq_nmod_mpoly_t r,
                        const fq_nmod_mpoly_t g, const fq_nmod_mpoly_ctx_t ctx)
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
    {
        int divides;

        if (bits <= FLINT_BITS)
            divides = mpoly_monomial_divides_test(rexp + i*N, gexp + 0*N, N, mask);
        else
            divides = mpoly_monomial_divides_mp_test(rexp + i*N, gexp + 0*N, N, bits);

        if (divides)
        {
            flint_printf("fq_nmod_mpoly_remainder_strongtest FAILED i = %wd\n", i);
            flint_printf("rem ");fq_nmod_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
            flint_printf("den ");fq_nmod_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
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

