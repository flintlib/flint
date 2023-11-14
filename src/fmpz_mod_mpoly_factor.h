/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MPOLY_FACTOR_H
#define FMPZ_MOD_MPOLY_FACTOR_H

#ifdef FMPZ_MOD_MPOLY_FACTOR_INLINES_C
#define FMPZ_MOD_MPOLY_FACTOR_INLINE
#else
#define FMPZ_MOD_MPOLY_FACTOR_INLINE static inline
#endif

#include "thread_pool.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_mpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_init(fmpz_mod_mpoly_factor_t f,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_init_set_ui(f->constant, 1);
    f->poly  = NULL;
    f->exp   = NULL;
    f->num   = 0;
    f->alloc = 0;
}

void fmpz_mod_mpoly_factor_init2(fmpz_mod_mpoly_factor_t f,
                                  slong alloc, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_factor_realloc(fmpz_mod_mpoly_factor_t f,
                                  slong alloc, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_factor_fit_length(fmpz_mod_mpoly_factor_t f,
                                    slong len, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_factor_clear(fmpz_mod_mpoly_factor_t f,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
slong fmpz_mod_mpoly_factor_length(const fmpz_mod_mpoly_factor_t f,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    return f->num;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_get_constant_fmpz(fmpz_t c,
               const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_set(c, f->constant);
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_get_base(fmpz_mod_mpoly_t p,
     const fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fmpz_mod_mpoly_set(p, f->poly + i, ctx);
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_swap_base(fmpz_mod_mpoly_t p,
            fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fmpz_mod_mpoly_swap(p, f->poly + i, ctx);
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
slong fmpz_mod_mpoly_factor_get_exp_si(fmpz_mod_mpoly_factor_t f,
                                       slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    return fmpz_get_si(f->exp + i);
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_swap(fmpz_mod_mpoly_factor_t f,
                     fmpz_mod_mpoly_factor_t g, const fmpz_mod_mpoly_ctx_t ctx)
{
   fmpz_mod_mpoly_factor_struct t = *f;
   *f = *g;
   *g = t;
}

void fmpz_mod_mpoly_factor_set(fmpz_mod_mpoly_factor_t f,
              const fmpz_mod_mpoly_factor_t g, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_factor_print_pretty(const fmpz_mod_mpoly_factor_t f,
                           const char ** vars, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_factor_content(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_factor_squarefree(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_factor(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_factor_sort(fmpz_mod_mpoly_factor_t f,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_factor_cmp(const fmpz_mod_mpoly_factor_t A,
              const fmpz_mod_mpoly_factor_t B, const fmpz_mod_mpoly_ctx_t ctx);

/******************************************************************************

   Internal functions (guaranteed to change without notice)

******************************************************************************/

/* fmpz_mod_poly extras ******************************************************/

FMPZ_MOD_MPOLY_FACTOR_INLINE
slong _fmpz_mod_poly_degree(const fmpz_mod_poly_t a)
{
    return a->length - 1;
}

void fmpz_mod_poly_scalar_addmul_fmpz_mod(fmpz_mod_poly_t A,
           const fmpz_mod_poly_t B, const fmpz_mod_poly_t C, const fmpz_t d0,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_addmul_linear(fmpz_mod_poly_t A,
                    const fmpz_mod_poly_t B, const fmpz_mod_poly_t C,
                   const fmpz_t d1, const fmpz_t d0, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(
         fmpz_mod_poly_t A, slong k, const fmpz_t c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_eval_pow(fmpz_t eval, const fmpz_mod_poly_t P,
                           fmpz_mod_poly_t alphapow, const fmpz_mod_ctx_t ctx);

void fmpz_mod_poly_eval2_pow(fmpz_t evalp, fmpz_t evalm,
                            const fmpz_mod_poly_t P, fmpz_mod_poly_t alphapow,
                                                     const fmpz_mod_ctx_t ctx);

/* fmpz_mod_mat extras *******************************************************/

int fmpz_mod_mat_is_reduced(const fmpz_mod_mat_t N);

void fmpz_mod_mat_init_nullspace_tr(fmpz_mod_mat_t X,
                                 fmpz_mod_mat_t tmp, const fmpz_mod_ctx_t ctx);

/*****************************************************************************/

typedef struct
{
    fmpz_mod_poly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_mod_bpoly_struct;

typedef fmpz_mod_bpoly_struct fmpz_mod_bpoly_t[1];


typedef struct
{
    fmpz_mod_bpoly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_mod_tpoly_struct;

typedef fmpz_mod_tpoly_struct fmpz_mod_tpoly_t[1];


typedef struct
{
    ulong * exps;
    fmpz * coeffs;
    slong length;
    slong alloc;
} fmpz_mod_polyu_struct;

typedef fmpz_mod_polyu_struct fmpz_mod_polyu_t[1];


typedef struct
{
    fmpz_mod_poly_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
} fmpz_mod_polyun_struct;

typedef fmpz_mod_polyun_struct fmpz_mod_polyun_t[1];

/*
    fmpz_mod_mpolyu_t
    sparse univariates with fmpz_mod_mpoly_t coefficients
        with uniform bits and LEX ordering
*/
typedef struct
{
   fmpz_mod_mpoly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;    /* default bits to construct coeffs */
} fmpz_mod_mpolyu_struct;
typedef fmpz_mod_mpolyu_struct fmpz_mod_mpolyu_t[1];

/*
    fmpz_mod_mpolyn_t
    sparse multivariates with fmpz_mod_poly_t coefficients
        with LEX ordering
*/
typedef struct
{
   fmpz_mod_poly_struct * coeffs;
   ulong * exps;
   slong alloc;
   slong length;
   slong bits;
} fmpz_mod_mpolyn_struct;

typedef fmpz_mod_mpolyn_struct fmpz_mod_mpolyn_t[1];

/*****************************************************************************/

int fmpz_mod_mpoly_factor_separable(fmpz_mod_mpoly_factor_t f,
            const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, int sep);

int fmpz_mod_mpoly_factor_expand(fmpz_mod_mpoly_t A,
              const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
int fmpz_mod_mpoly_factor_matches(const fmpz_mod_mpoly_t a,
               const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
{
    int matches;
    fmpz_mod_mpoly_t t;
    fmpz_mod_mpoly_init(t, ctx);
    fmpz_mod_mpoly_factor_expand(t, f, ctx);
    matches = fmpz_mod_mpoly_equal(t, a, ctx);
    fmpz_mod_mpoly_clear(t, ctx);
    return matches;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_append_fmpz_swap(fmpz_mod_mpoly_factor_t f,
            fmpz_mod_mpoly_t A, const fmpz_t e, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i = f->num;
    fmpz_mod_mpoly_factor_fit_length(f, i + 1, ctx);
    fmpz_mod_mpoly_swap(f->poly + i, A, ctx);
    fmpz_set(f->exp + i, e);
    f->num = i + 1;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpoly_factor_one(fmpz_mod_mpoly_factor_t f,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_one(f->constant);
	f->num = 0;
}

void _fmpz_mod_mpoly_get_lead0(
    fmpz_mod_mpoly_t c,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_set_lead0(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_t c,
    const fmpz_mod_mpoly_ctx_t ctx);

/*****************************************************************************/

void _fmpz_mod_mpoly_factor_set_nmod_mpoly_factor(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_ctx_t ctx,
    const nmod_mpoly_factor_t nf,
    const nmod_mpoly_ctx_t nctx);


/* stack *********************************************************************/

typedef struct
{
    fmpz_mod_poly_struct ** array;
    slong alloc;
    slong top;
} fmpz_mod_poly_stack_struct;

typedef fmpz_mod_poly_stack_struct fmpz_mod_poly_stack_t[1];


typedef struct
{
    fmpz_mod_bpoly_struct ** array;
    slong alloc;
    slong top;
} fmpz_mod_bpoly_stack_struct;

typedef fmpz_mod_bpoly_stack_struct fmpz_mod_bpoly_stack_t[1];

typedef struct
{
    fmpz_mod_polyun_struct ** array;
    slong alloc;
    slong top;
} fmpz_mod_polyun_stack_struct;

typedef fmpz_mod_polyun_stack_struct fmpz_mod_polyun_stack_t[1];

typedef struct
{
    fmpz_mod_mpolyn_struct ** array;
    slong alloc;
    slong top;
    flint_bitcnt_t bits;
} fmpz_mod_mpolyn_stack_struct;

typedef fmpz_mod_mpolyn_stack_struct fmpz_mod_mpolyn_stack_t[1];

typedef struct {
    fmpz_mod_poly_stack_t poly_stack;
    fmpz_mod_bpoly_stack_t bpoly_stack;
} fmpz_mod_poly_bpoly_stack_struct;

typedef fmpz_mod_poly_bpoly_stack_struct fmpz_mod_poly_bpoly_stack_t[1];

typedef struct {
    fmpz_mod_poly_stack_t poly_stack;
    fmpz_mod_polyun_stack_t polyun_stack;
} fmpz_mod_poly_polyun_stack_struct;

typedef fmpz_mod_poly_polyun_stack_struct fmpz_mod_poly_polyun_stack_t[1];

typedef struct {
    fmpz_mod_poly_stack_t poly_stack;
    fmpz_mod_polyun_stack_t polyun_stack;
    fmpz_mod_mpolyn_stack_t mpolyn_stack;
} fmpz_mod_poly_polyun_mpolyn_stack_struct;

typedef fmpz_mod_poly_polyun_mpolyn_stack_struct
                                       fmpz_mod_poly_polyun_mpolyn_stack_t[1];

void fmpz_mod_poly_stack_init(fmpz_mod_poly_stack_t S);

void fmpz_mod_poly_stack_clear(fmpz_mod_poly_stack_t S);

fmpz_mod_poly_struct ** fmpz_mod_poly_stack_fit_request(
                                             fmpz_mod_poly_stack_t S, slong k);

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_poly_struct ** fmpz_mod_poly_stack_request(
                                              fmpz_mod_poly_stack_t S, slong k)
{
    fmpz_mod_poly_struct ** poly_top;
    poly_top = fmpz_mod_poly_stack_fit_request(S, k);
    S->top += k;
    return poly_top;
}

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_poly_struct * fmpz_mod_poly_stack_take_top(fmpz_mod_poly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    fmpz_mod_poly_struct ** poly_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    poly_top = S->array + S->top;
    S->top += 1;
    return poly_top[0];
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_poly_stack_give_back(fmpz_mod_poly_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_poly_stack_size(const fmpz_mod_poly_stack_t S)
{
    return S->top;
}

void fmpz_mod_bpoly_stack_init(fmpz_mod_bpoly_stack_t S);

void fmpz_mod_bpoly_stack_clear(fmpz_mod_bpoly_stack_t S);

fmpz_mod_bpoly_struct ** fmpz_mod_bpoly_stack_fit_request(
                                            fmpz_mod_bpoly_stack_t S, slong k);

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_bpoly_struct ** fmpz_mod_bpoly_stack_request(
                                             fmpz_mod_bpoly_stack_t S, slong k)
{
    fmpz_mod_bpoly_struct ** bpoly_top;
    bpoly_top = fmpz_mod_bpoly_stack_fit_request(S, k);
    S->top += k;
    return bpoly_top;
}

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_bpoly_struct * fmpz_mod_bpoly_stack_take_top(fmpz_mod_bpoly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    fmpz_mod_bpoly_struct ** bpoly_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    bpoly_top = S->array + S->top;
    S->top += 1;
    return bpoly_top[0];
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_bpoly_stack_give_back(fmpz_mod_bpoly_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_bpoly_stack_size(const fmpz_mod_bpoly_stack_t S)
{
    return S->top;
}


void fmpz_mod_polyun_stack_init(fmpz_mod_polyun_stack_t S);

void fmpz_mod_polyun_stack_clear(fmpz_mod_polyun_stack_t S);

fmpz_mod_polyun_struct ** fmpz_mod_polyun_stack_fit_request(
                                           fmpz_mod_polyun_stack_t S, slong k);

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_polyun_struct ** fmpz_mod_polyun_stack_request(
                                            fmpz_mod_polyun_stack_t S, slong k)
{
    fmpz_mod_polyun_struct ** polyun_top;
    polyun_top = fmpz_mod_polyun_stack_fit_request(S, k);
    S->top += k;
    return polyun_top;
}

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_polyun_struct * fmpz_mod_polyun_stack_take_top(fmpz_mod_polyun_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    fmpz_mod_polyun_struct ** polyun_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    polyun_top = S->array + S->top;
    S->top += 1;
    return polyun_top[0];
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_polyun_stack_give_back(fmpz_mod_polyun_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_polyun_stack_size(const fmpz_mod_polyun_stack_t S)
{
    return S->top;
}



void fmpz_mod_mpolyn_stack_init(fmpz_mod_mpolyn_stack_t S,
                          flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_stack_clear(fmpz_mod_mpolyn_stack_t S,
                                               const fmpz_mod_mpoly_ctx_t ctx);

fmpz_mod_mpolyn_struct ** fmpz_mod_mpolyn_stack_fit_request(
           fmpz_mod_mpolyn_stack_t S, slong k, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_mpolyn_struct ** fmpz_mod_mpolyn_stack_request(
            fmpz_mod_mpolyn_stack_t S, slong k, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpolyn_struct ** mpolyn_top;
    mpolyn_top = fmpz_mod_mpolyn_stack_fit_request(S, k, ctx);
    S->top += k;
    return mpolyn_top;
}

FMPZ_MOD_MPOLY_INLINE
fmpz_mod_mpolyn_struct * fmpz_mod_mpolyn_stack_take_top(fmpz_mod_mpolyn_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    fmpz_mod_mpolyn_struct ** mpolyn_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    mpolyn_top = S->array + S->top;
    S->top += 1;
    return mpolyn_top[0];
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_mpolyn_stack_give_back(fmpz_mod_mpolyn_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

FMPZ_MOD_MPOLY_INLINE
slong fmpz_mod_mpolyn_stack_size(const fmpz_mod_mpolyn_stack_t S)
{
    return S->top;
}

/* poly_vec ******************************************************************/

slong _fmpz_mod_poly_vec_max_degree(const fmpz_mod_poly_struct * A,
                                         slong Alen, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_vec_content(fmpz_mod_poly_t g,
         const fmpz_mod_poly_struct * A, slong Alen, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_vec_remove_content(fmpz_mod_poly_t g,
               fmpz_mod_poly_struct * A, slong Alen, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_vec_mul_poly(fmpz_mod_poly_struct * A,
                slong Alen, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_vec_divexact_poly(fmpz_mod_poly_struct * A,
                slong Alen, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_poly_vec_mul_fmpz_mod(fmpz_mod_poly_struct * A,
                         slong Alen, const fmpz_t g, const fmpz_mod_ctx_t ctx);

/* polyun ********************************************************************/

FMPZ_MOD_MPOLY_FACTOR_INLINE
ulong fmpz_mod_polyu1n_bidegree(const fmpz_mod_polyun_t A)
{
    ulong x, y;
    FLINT_ASSERT(A->length > 0);
    x = A->exps[0];
    y = A->coeffs[0].length - 1;
    return (x << (FLINT_BITS/2)) + y;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
const fmpz * fmpz_mod_polyun_leadcoeff(const fmpz_mod_polyun_t A)
{
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0].length > 0);
    return A->coeffs[0].coeffs + A->coeffs[0].length - 1;
}


FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_polyun_swap(fmpz_mod_polyun_t A, fmpz_mod_polyun_t B)
{
    fmpz_mod_polyun_struct t = *A;
    *A = *B;
    *B = t;
}

int fmpz_mod_polyun_is_canonical(const fmpz_mod_polyun_t A,
                                                    const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_polyun_init(fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx)
{
    A->length = 0;
    A->alloc = 0;
    A->coeffs = NULL;
    A->exps = NULL;
}

void fmpz_mod_polyun_clear(fmpz_mod_polyun_t A,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_polyun_realloc(fmpz_mod_polyun_t A,
                                          slong len, const fmpz_mod_ctx_t ctx);

void fmpz_mod_polyu2n_print_pretty(const fmpz_mod_polyun_t A,
                const char * var0, const char * var1, const char * varlast,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_polyun_equal(fmpz_mod_polyun_t A,
                          const fmpz_mod_polyun_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_polyun_set(fmpz_mod_polyun_t A,
                          const fmpz_mod_polyun_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_polyu3n_print_pretty( const fmpz_mod_polyun_t A,
                    const char * var0, const char * var1, const char * var2,
                               const char * varlast, const fmpz_mod_ctx_t ctx);

void fmpz_mod_polyu1n_print_pretty(const fmpz_mod_polyun_t A,
            const char * var0, const char * varlast, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_polyun_fit_length(fmpz_mod_polyun_t A, slong len,
                                                      const fmpz_mod_ctx_t ctx)
{
    if (len > A->alloc)
        fmpz_mod_polyun_realloc(A, len, ctx);
}

void fmpz_mod_polyun_one(fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mpoly_get_polyu1n(fmpz_mod_polyun_t A,
                            const fmpz_mod_mpoly_t B, slong varx, slong vary,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_polyu1n(fmpz_mod_mpoly_t B,
                            const fmpz_mod_polyun_t A, slong varx, slong vary,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* mpolyn ********************************************************************/

void fmpz_mod_mpolyn_init(fmpz_mod_mpolyn_t A, flint_bitcnt_t bits,
                                               const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_INLINE void fmpz_mod_mpolyn_swap(fmpz_mod_mpolyn_t A,
                           fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx)
{
   fmpz_mod_mpolyn_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mod_mpolyn_fit_length(fmpz_mod_mpolyn_t A,
                                 slong length, const fmpz_mod_mpoly_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
const fmpz * fmpz_mod_mpolyn_leadcoeff(const fmpz_mod_mpolyn_t A)
{
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0].length > 0);
    return A->coeffs[0].coeffs + A->coeffs[0].length - 1;
}

int fmpz_mod_mpolyn_is_canonical(const fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

slong fmpz_mod_mpolyn_lastdeg(const fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_clear(fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_one(fmpz_mod_mpolyn_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(fmpz_mod_mpolyn_t A,
                               const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyn_equal(const fmpz_mod_mpolyn_t A,
                   const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_print_pretty(const fmpz_mod_mpolyn_t poly,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_cvtfrom_mpolyn(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_cvtto_mpolyn(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_to_mpolyn_perm_deflate(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t nctx,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

void fmpz_mod_mpoly_from_mpolyn_perm_inflate(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_mpoly_ctx_t ctx,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t nctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride);

void fmpz_mod_mpolyn_set(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyn_is_nonzero_fmpz(
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyn_divides(
    fmpz_mod_mpolyn_t Q,
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx);

/* interp ********************************************************************/

void fmpz_mod_polyu1n_interp_reduce_2sm_poly(
    fmpz_mod_poly_t E,
    fmpz_mod_poly_t F,
    const fmpz_mod_polyun_t A,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_polyu1n_interp_lift_2sm_poly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx);

int fmpz_mod_polyu1n_interp_crt_2sm_poly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    fmpz_mod_polyun_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_mpolyn_interp_reduce_sm_poly(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpolyn_t A,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_interp_lift_sm_poly(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyn_interp_crt_sm_poly(
    slong * lastdeg_,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(
    fmpz_mod_mpolyn_t E,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t A,
    slong var,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(
    slong * lastdeg,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyn_interp_lift_sm_mpoly(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyn_interp_crt_sm_mpoly(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpoly_t A,
    fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx);


/* polyu *********************************************************************/

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_polyu_swap(
    fmpz_mod_polyu_t A,
    fmpz_mod_polyu_t B)
{
    fmpz_mod_polyu_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mod_polyu_init(fmpz_mod_polyu_t A);

void fmpz_mod_polyu_clear(fmpz_mod_polyu_t A);

void fmpz_mod_polyu_realloc(fmpz_mod_polyu_t A, slong len);

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_polyu_fit_length(
    fmpz_mod_polyu_t a,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    if (len > a->alloc)
        fmpz_mod_polyu_realloc(a, len);
}

void fmpz_mod_polyu3_degrees(
    slong * deg0,
    slong * deg1,
    slong * deg2,
    const fmpz_mod_polyu_t A);

void fmpz_mod_polyu3_print_pretty(
    const fmpz_mod_polyu_t A,
    const char * var0,
    const char * var1,
    const char * var2,
    const fmpz_mod_ctx_t ctx);

/* mpolyu ********************************************************************/

int fmpz_mod_mpolyu_is_canonical(const fmpz_mod_mpolyu_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyu3_print_pretty(const fmpz_mod_mpolyu_t A,
                    const char * var0, const char * var1, const char * var2,
                           const char ** vars, const fmpz_mod_mpoly_ctx_t ctx);


/* bpoly *********************************************************************/

int fmpz_mod_bpoly_is_canonical(const fmpz_mod_bpoly_t A,
                                                     const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_bpoly_init(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_mod_bpoly_clear(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_bpoly_swap(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}


FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_bpoly_get_coeff(fmpz_t c, const fmpz_mod_bpoly_t A,
                                  slong e0, slong e1, const fmpz_mod_ctx_t ctx)
{
    if (e0 >= A->length)
        fmpz_zero(c);
    else
        fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + e0, e1, ctx);
}


FMPZ_MOD_MPOLY_FACTOR_INLINE
slong fmpz_mod_bpoly_degree0(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    return A->length - 1;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_bpoly_normalise(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    while (A->length > 0 && fmpz_mod_poly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;
}

int fmpz_mod_bpoly_equal(
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_set(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx);


void fmpz_mod_bpoly_set_poly_gen1(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx);


void fmpz_mod_bpoly_set_poly_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_one(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MPOLY_FACTOR_INLINE
int fmpz_mod_bpoly_is_one(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    return A->length == 1 && fmpz_mod_poly_is_one(A->coeffs + 0, ctx);
}

slong fmpz_mod_bpoly_degree1(const fmpz_mod_bpoly_t A,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_print_pretty(const fmpz_mod_bpoly_t A,
               const char * xvar, const char * yvar, const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_fit_length(fmpz_mod_bpoly_t A, slong len,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_set_coeff(fmpz_mod_bpoly_t A, slong xi, slong yi,
                                     const fmpz_t c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_zero(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_reverse_vars(fmpz_mod_bpoly_t A,
                           const fmpz_mod_bpoly_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_taylor_shift_gen1(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_sub(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_add(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_make_primitive(
    fmpz_mod_poly_t g,
    fmpz_mod_bpoly_t A,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_mul(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_mul_series(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    slong order,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_divrem_series(
    fmpz_mod_bpoly_t Q,
    fmpz_mod_bpoly_t R,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    slong order,
    const fmpz_mod_ctx_t ctx);

int fmpz_mod_bpoly_divides(
    fmpz_mod_bpoly_t Q,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_taylor_shift_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_derivative_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_bpoly_make_monic_series(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    slong order,
    const fmpz_mod_ctx_t ctx);


FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_tpoly_init(fmpz_mod_tpoly_t A, const fmpz_mod_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FMPZ_MOD_MPOLY_INLINE
void fmpz_mod_tpoly_swap(fmpz_mod_tpoly_t A, fmpz_mod_tpoly_t B,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_tpoly_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mod_tpoly_fit_length(fmpz_mod_tpoly_t A, slong len,
                                                     const fmpz_mod_ctx_t ctx);

void fmpz_mod_tpoly_clear(fmpz_mod_tpoly_t A, const fmpz_mod_ctx_t ctx);


void fmpz_mod_mpoly_get_fmpz_mod_bpoly(fmpz_mod_bpoly_t A,
                            const fmpz_mod_mpoly_t B, slong var0, slong var1,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_set_fmpz_mod_bpoly(fmpz_mod_mpoly_t A,
        flint_bitcnt_t Abits, const fmpz_mod_bpoly_t B, slong var0, slong var1,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_bpoly_factor_smprime(fmpz_mod_poly_t c,
                    fmpz_mod_tpoly_t F, fmpz_mod_bpoly_t B, int allow_shift,
                                                     const fmpz_mod_ctx_t ctx);

/*****************************************************************************/

int _fmpz_mod_zip_vand_solve(
    fmpz * coeffs,             /* in Fp: size mlength */
    const fmpz * monomials,    /* in Fp: size mlength */
    slong mlength,
    const fmpz * evals,        /* in Fp: size elength */
    slong elength,
    const fmpz * master,       /* in Fp: size mlength + 1 */
    fmpz * scratch,            /* in Fp: size mlength */
    const fmpz_mod_ctx_t ctx);

void _fmpz_mod_zip_eval_step(
    fmpz_t ev,
    fmpz * cur,            /* in Fp */
    const fmpz * inc,      /* in Fp */
    const fmpz * coeffs,   /* in Fp */
    slong length,
    const fmpz_mod_ctx_t ctx);


/*****************************************************************************/

typedef struct
{
    fmpz_mod_mpoly_struct * coeffs;
    slong alloc;
    slong length;
} fmpz_mod_mpolyv_struct;

typedef fmpz_mod_mpolyv_struct fmpz_mod_mpolyv_t[1];

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpolyv_init(fmpz_mod_mpolyv_t A, const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FMPZ_MOD_MPOLY_FACTOR_INLINE
void fmpz_mod_mpolyv_swap(fmpz_mod_mpolyv_t A, fmpz_mod_mpolyv_t B,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
   fmpz_mod_mpolyv_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mod_mpolyv_clear(fmpz_mod_mpolyv_t A,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyv_print_pretty(const fmpz_mod_mpolyv_t poly,
                              const char ** x, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyv_fit_length(fmpz_mod_mpolyv_t A, slong length,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpolyv_set_coeff(fmpz_mod_mpolyv_t A, slong i,
                           fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_to_mpolyv(fmpz_mod_mpolyv_t A,
                    const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t xalpha,
                                               const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_from_mpolyv(fmpz_mod_mpoly_t A,
                 flint_bitcnt_t Abits, const fmpz_mod_mpolyv_t B,
                const fmpz_mod_mpoly_t xalpha, const fmpz_mod_mpoly_ctx_t ctx);

int _fmpz_mod_mpoly_vec_content_mpoly(fmpz_mod_mpoly_t g,
                                const fmpz_mod_mpoly_struct * A, slong Alen,
                                              const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_vec_divexact_mpoly(fmpz_mod_mpoly_struct * A,
         slong Alen, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_vec_mul_mpoly(fmpz_mod_mpoly_struct * A,
         slong Alen, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx);

/*****************************************************************************/

int _fmpz_mod_mpoly_factor_separable(fmpz_mod_mpoly_factor_t f,
            const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, int sep);

int fmpz_mod_mpoly_factor_lcc_wang(fmpz_mod_mpoly_struct * lc_divs,
             const fmpz_mod_mpoly_factor_t lcAfac, const fmpz_mod_poly_t Auc,
             const fmpz_mod_bpoly_struct * Auf, slong r,
           const fmpz_mod_poly_struct * alpha, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_factor_irred_smprime_zassenhaus(
                            fmpz_mod_mpolyv_t fac, const fmpz_mod_mpoly_t A,
                           const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state);

int fmpz_mod_mpoly_factor_irred_smprime_wang(fmpz_mod_mpolyv_t fac,
                const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_factor_t lcAfac,
                const fmpz_mod_mpoly_t lcA, const fmpz_mod_mpoly_ctx_t ctx,
                                                           flint_rand_t state);

int fmpz_mod_mpoly_factor_irred_smprime_zippel(fmpz_mod_mpolyv_t fac,
                const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_factor_t lcAfac,
                const fmpz_mod_mpoly_t lcA, const fmpz_mod_mpoly_ctx_t ctx,
                                                           flint_rand_t state);

/*****************************************************************************/

void fmpz_mod_mpoly_compression_do(fmpz_mod_mpoly_t L,
                 const fmpz_mod_mpoly_ctx_t Lctx, fmpz * Acoeffs, slong Alen,
                                                        mpoly_compression_t M);

void fmpz_mod_mpoly_compression_undo(fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t Actx, fmpz_mod_mpoly_t L,
                       const fmpz_mod_mpoly_ctx_t Lctx, mpoly_compression_t M);

/*****************************************************************************/

typedef struct {
    flint_bitcnt_t bits;
    slong w;
    slong r;
    fmpz_mod_poly_struct * inv_prod_dbetas;
    fmpz_mod_mpoly_struct * inv_prod_dbetas_mvar;
    fmpz_mod_poly_struct * dbetas;
    fmpz_mod_mpoly_struct * dbetas_mvar;
    fmpz_mod_mpoly_struct * prod_mbetas;
    fmpz_mod_mpolyv_struct * prod_mbetas_coeffs;
    fmpz_mod_mpoly_struct * mbetas;
    fmpz_mod_mpoly_struct * deltas;
    fmpz_mod_mpoly_struct * xalpha;
    fmpz_mod_mpoly_struct * q;
    fmpz_mod_mpoly_geobucket_struct * G;
    fmpz_mod_mpoly_struct * qt;
    fmpz_mod_mpoly_struct * newt;
    fmpz_mod_mpolyv_struct * delta_coeffs;
    fmpz_mod_mpoly_t T;
    fmpz_mod_mpoly_t Q;
    fmpz_mod_mpoly_t R;
} fmpz_mod_mpoly_pfrac_struct;

typedef fmpz_mod_mpoly_pfrac_struct fmpz_mod_mpoly_pfrac_t[1];


int fmpz_mod_mpoly_pfrac_init(fmpz_mod_mpoly_pfrac_t I,
    flint_bitcnt_t bits, slong l, slong r, const fmpz_mod_mpoly_struct * betas,
                          const fmpz * alpha, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_pfrac_clear(fmpz_mod_mpoly_pfrac_t I,
                                              const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_pfrac(slong r, fmpz_mod_mpoly_t t, const slong * deg,
                     fmpz_mod_mpoly_pfrac_t I, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_hlift(slong m, fmpz_mod_mpoly_struct * f, slong r,
            const fmpz * alpha, const fmpz_mod_mpoly_t A, const slong * degs,
                                               const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_bpoly_pfrac(slong r, fmpz_mod_bpoly_struct * C,
          slong * C_deg1_bound, fmpz_mod_bpoly_t A, fmpz_mod_bpoly_struct * B,
                                                     const fmpz_mod_ctx_t ctx);

int fmpz_mod_bpoly_hlift2(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B0,
                  fmpz_mod_bpoly_t B1, const fmpz_t alpha, slong degree_inner,
                     const fmpz_mod_ctx_t ctx, fmpz_mod_poly_bpoly_stack_t St);

int fmpz_mod_bpoly_hlift(slong r, fmpz_mod_bpoly_t A,
           fmpz_mod_bpoly_struct * B, const fmpz_t alpha, slong degree_inner,
                     const fmpz_mod_ctx_t ctx, fmpz_mod_poly_bpoly_stack_t St);

int fmpz_mod_polyu3_hlift(slong r, fmpz_mod_polyun_struct * BB,
              fmpz_mod_polyu_t A, fmpz_mod_polyu_struct * B, const fmpz_t beta,
                                 slong degree_inner, const fmpz_mod_ctx_t ctx);

int fmpz_mod_mpoly_hlift_zippel(slong m, fmpz_mod_mpoly_struct * B,
     slong r, const fmpz * alpha, const fmpz_mod_mpoly_t A, const slong * degs,
                           const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state);

int fmpz_mod_mpoly_factor_algo(fmpz_mod_mpoly_factor_t f,
                    const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx,
                                                            unsigned int algo);

int fmpz_mod_mpoly_factor_zassenhaus(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_factor_wang(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpoly_factor_zippel(fmpz_mod_mpoly_factor_t f,
                     const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx);

int _fmpz_mod_mpoly_evaluate_rest_fmpz_mod_poly(fmpz_mod_poly_struct * E,
                slong * starts, slong * ends, slong * stops, ulong * es,
                const fmpz * Acoeffs, const ulong * Aexps, slong Alen, slong var,
                const fmpz_mod_poly_struct * alphas,
                const slong * offsets, const slong * shifts,
                   slong N, ulong mask, slong nvars, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_mpoly_eval_rest_to_fmpz_mod_bpoly(fmpz_mod_bpoly_t E,
             const fmpz_mod_mpoly_t A, const fmpz_mod_poly_struct * alphabetas,
                                                const fmpz_mod_mpoly_ctx_t ctx);

void _fmpz_mod_mpoly_set_fmpz_mod_bpoly_var1_zero(fmpz_mod_mpoly_t A,
                     flint_bitcnt_t Abits, const fmpz_mod_bpoly_t B, slong var,
                                               const fmpz_mod_mpoly_ctx_t ctx);

/* gcd ***********************************************************************/

int _fmpz_mod_mpoly_gcd_algo(fmpz_mod_mpoly_t G,
                            fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar,
                            const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                            const fmpz_mod_mpoly_ctx_t ctx, unsigned int algo);

int fmpz_mod_polyu1n_gcd_brown_smprime(
    fmpz_mod_polyun_t G,
    fmpz_mod_polyun_t Abar,
    fmpz_mod_polyun_t Bbar,
    fmpz_mod_polyun_t A,
    fmpz_mod_polyun_t B,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_stack_t St_poly,
    fmpz_mod_polyun_stack_t St_polyun);

int fmpz_mod_mpolyn_gcd_brown_smprime(
    fmpz_mod_mpolyn_t G,
    fmpz_mod_mpolyn_t Abar,
    fmpz_mod_mpolyn_t Bbar,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I,
    fmpz_mod_poly_polyun_mpolyn_stack_t St);

int fmpz_mod_mpolyl_gcdp_zippel(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    fmpz_mod_mpoly_t A,
    fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx,
    flint_rand_t state);

int fmpz_mod_mpolyl_gcd_zippel2_smprime(
    fmpz_mod_mpoly_t rG, const slong * rGdegs,
    fmpz_mod_mpoly_t rAbar,
    fmpz_mod_mpoly_t rBbar,
    const fmpz_mod_mpoly_t A, const slong * Adegs,
    const fmpz_mod_mpoly_t B, const slong * Bdegs,
    const fmpz_mod_mpoly_t gamma, const slong * gammadegs,
    const fmpz_mod_mpoly_ctx_t ctx);

int fmpz_mod_mpolyl_gcd_hensel_smprime(
    fmpz_mod_mpoly_t G, slong Gdeg,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx);

/* pow_cache */

void fmpz_mod_pow_cache_start(const fmpz_t b,
                                  fmpz_mod_poly_t c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_pow_cache_mulpow_ui(fmpz_t a, const fmpz_t b,
                         ulong e, fmpz_mod_poly_t c, const fmpz_mod_ctx_t ctx);

/* zip helpers */

void mpoly_monomial_evals_fmpz_mod(
    fmpz_mod_poly_t EH,
    const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
    fmpz_mod_poly_struct * alpha_caches,
    slong start,
    slong stop,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fpctx);

void mpoly1_monomial_evals_fmpz_mod(
    fmpz_mod_polyun_t EH,
    const ulong * Aexps, flint_bitcnt_t Abits,
    const ulong * Amarks, slong Amarkslen,
    fmpz_mod_poly_struct * alpha_caches,
    slong m,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fpctx);

void mpoly2_monomial_evals_fmpz_mod(
    fmpz_mod_polyun_t EH,
    const ulong * Aexps, flint_bitcnt_t Abits,
    ulong * Amarks, slong Amarkslen,
    fmpz_mod_poly_struct * alpha_caches,
    slong m,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fpctx);

void fmpz_mod_mpoly_mock_eval_coeff(
    fmpz_mod_polyun_t mock,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_polyun_t Aeh_inc,
    const fmpz_mod_mpoly_ctx_t ctx);

slong fmpz_mod_polyun_product_roots(
    fmpz_mod_polyun_t M,
    const fmpz_mod_polyun_t H,
    const fmpz_mod_ctx_t ctx);

void fmpz_mod_polyun_zip_start(
    fmpz_mod_polyun_t Z,
    fmpz_mod_polyun_t H,
    slong req_images,
    const fmpz_mod_ctx_t fctx);

void fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(
    fmpz_mod_polyun_t E,
    fmpz_mod_polyun_t Acur,
    const fmpz_mod_polyun_t Ainc,
    const fmpz_mod_polyun_t Acoeff,
    const fmpz_mod_ctx_t ctx);

int fmpz_mod_polyun_zip_solve(
    fmpz_mod_mpoly_t A,
    fmpz_mod_polyun_t Z,
    fmpz_mod_polyun_t H,
    fmpz_mod_polyun_t M,
    const fmpz_mod_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

