/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_H
#define FQ_NMOD_H

#ifdef FQ_NMOD_INLINES_C
#define FQ_NMOD_INLINE FLINT_DLL
#define FQ_TEMPLATES_INLINE FLINT_DLL
#else
#define FQ_NMOD_INLINE static __inline__
#define FQ_TEMPLATES_INLINE static __inline__
#endif

#include "fmpz_mini.h"
#include "nmod.h"
#include "nmod_poly_mini.h"

/* Data types and context ****************************************************/
#ifdef __cplusplus
extern "C" {
#endif

FLINT_DLL void fq_nmod_ctx_init(fq_nmod_ctx_t ctx,
                      const fmpz_t p, slong d, const char *var);

FLINT_DLL int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx,
                             const fmpz_t p, slong d, const char *var);

FLINT_DLL void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx,
                             const fmpz_t p, slong d, const char *var);

FLINT_DLL void fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx,
                              const nmod_poly_t modulus,
                              const char *var);

FLINT_DLL void fq_nmod_ctx_randtest(fq_nmod_ctx_t ctx, flint_rand_t state);

FLINT_DLL void fq_nmod_ctx_randtest_reducible(fq_nmod_ctx_t ctx, flint_rand_t state);

FLINT_DLL void fq_nmod_ctx_clear(fq_nmod_ctx_t ctx);

FQ_NMOD_INLINE const nmod_poly_struct* fq_nmod_ctx_modulus(const fq_nmod_ctx_t ctx)
{
    return ctx->modulus;
}

FQ_NMOD_INLINE slong fq_nmod_ctx_degree(const fq_nmod_ctx_t ctx)
{
    return ctx->modulus->length - 1;
}

#define fq_nmod_ctx_prime(ctx)  (&((ctx)->p))

FLINT_DLL void fq_nmod_ctx_order(fmpz_t f, const fq_nmod_ctx_t ctx);

/* Memory managment  *********************************************************/

FQ_NMOD_INLINE void fq_nmod_init(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_poly_init2_preinv(rop, ctx->mod.n, ctx->mod.ninv, fq_nmod_ctx_degree(ctx));
}

FQ_NMOD_INLINE void fq_nmod_init2(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_poly_init2_preinv(rop, ctx->mod.n, ctx->mod.ninv, fq_nmod_ctx_degree(ctx));
}

FQ_NMOD_INLINE void fq_nmod_clear(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_poly_clear(rop);
}

FQ_NMOD_INLINE 
void _fq_nmod_sparse_reduce(ulong *R, slong lenR, const fq_nmod_ctx_t ctx)
{
    slong i, k;
    const slong d = ctx->j[ctx->len - 1];

    NMOD_VEC_NORM(R, lenR);

    for (i = lenR - 1; i >= d; i--)
    {
        for (k = ctx->len - 2; k >= 0; k--)
        {
            /* TODO clean this mess up */
            R[ctx->j[k] + i - d] = n_submod(R[ctx->j[k] + i - d],
                                            n_mulmod2_preinv(R[i], ctx->a[k], ctx->mod.n, ctx->mod.ninv),
                                            ctx->mod.n);
        }
        R[i] = UWORD(0);
    }
}

FLINT_DLL void _fq_nmod_dense_reduce(ulong* R, slong lenR, const fq_nmod_ctx_t ctx);

FQ_NMOD_INLINE void _fq_nmod_reduce(ulong* R, slong lenR, const fq_nmod_ctx_t ctx)
{
    if (ctx->sparse_modulus)
        _fq_nmod_sparse_reduce(R, lenR, ctx);
    else
        _fq_nmod_dense_reduce(R, lenR, ctx);    
}

FQ_NMOD_INLINE void fq_nmod_reduce(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(rop->length <= 2*(ctx->modulus->length - 1));
    _fq_nmod_reduce(rop->coeffs, rop->length, ctx);
    rop->length = FLINT_MIN(rop->length, ctx->modulus->length - 1);
    _nmod_poly_normalise(rop);
}

/* Basic arithmetic **********************************************************/

FLINT_DLL void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1,
                                 const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_sub(fq_nmod_t rop, const fq_nmod_t op1,
                                 const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_sub_one(fq_nmod_t rop,
                                 const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_neg(fq_nmod_t rop,
                                 const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mul(fq_nmod_t rop,
            const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mul_fmpz(fq_nmod_t rop,
                  const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mul_si(fq_nmod_t rop,
                         const fq_nmod_t op, slong x, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mul_ui(fq_nmod_t rop,
                         const fq_nmod_t op, ulong x, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_sqr(fq_nmod_t rop,
                                  const fq_nmod_t op, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_inv(fq_nmod_t rop,
                                 const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_pow(ulong *rop, const ulong *op,
                           slong len, const fmpz_t e, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_pow(fq_nmod_t rop, const fq_nmod_t op1,
                                      const fmpz_t e, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_pow_ui(fq_nmod_t rop,
                  const fq_nmod_t op1, const ulong e, const fq_nmod_ctx_t ctx);

/* Roots ********************************************************************/

FLINT_DLL int fq_nmod_sqrt(fq_nmod_t rop, const fq_nmod_t op,
                                                      const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_pth_root(fq_nmod_t rop,
                                 const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

FLINT_DLL int fq_nmod_is_square(const fq_nmod_t op, const fq_nmod_ctx_t ctx);

/* Randomisation *************************************************************/

FLINT_DLL void fq_nmod_randtest(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_randtest_dense(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_randtest_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_rand(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_rand_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);


/* Comparison ****************************************************************/

FQ_NMOD_INLINE int fq_nmod_equal(const fq_nmod_t op1, const fq_nmod_t op2,
                                    const fq_nmod_ctx_t ctx)
{
    return nmod_poly_equal(op1, op2);
}

FQ_NMOD_INLINE int fq_nmod_is_zero(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    return nmod_poly_is_zero(op);
}

FQ_NMOD_INLINE int fq_nmod_is_one(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    return nmod_poly_is_one(op);
}

FLINT_DLL int fq_nmod_cmp(const fq_nmod_t a, const fq_nmod_t b,
                                                      const fq_nmod_ctx_t ctx);


/* Assignments and conversions ***********************************************/

FQ_NMOD_INLINE void fq_nmod_set(fq_nmod_t rop, const fq_nmod_t op,
                                   const fq_nmod_ctx_t ctx)
{
    nmod_poly_set(rop, op);
}

FLINT_DLL void fq_nmod_set_fmpz(fq_nmod_t rop, const fmpz_t x, const fq_nmod_ctx_t ctx);

FQ_NMOD_INLINE void fq_nmod_set_si(fq_nmod_t rop, const slong x, const fq_nmod_ctx_t ctx)
{
    ulong rx = x < 0 ? -x : x;
    rx =  n_mod2_preinv(rx, ctx->mod.n, ctx->mod.ninv);
    if (x < 0)
        rx = ctx->mod.n - rx;

    nmod_poly_zero(rop);
    nmod_poly_set_coeff_ui(rop, 0, rx);
}

FQ_NMOD_INLINE void fq_nmod_set_ui(fq_nmod_t rop, const ulong x, const fq_nmod_ctx_t ctx)
{
    nmod_poly_zero(rop);
    nmod_poly_set_coeff_ui(rop, 0, n_mod2_preinv(x, ctx->mod.n, ctx->mod.ninv));
}

FQ_NMOD_INLINE void fq_nmod_swap(fq_nmod_t op1, fq_nmod_t op2,
                                    const fq_nmod_ctx_t ctx)
{
    nmod_poly_swap(op1, op2);
}

FQ_NMOD_INLINE void fq_nmod_zero(fq_nmod_t rop,  const fq_nmod_ctx_t ctx)
{
    nmod_poly_zero(rop);
}

FQ_NMOD_INLINE void fq_nmod_one(fq_nmod_t rop,  const fq_nmod_ctx_t ctx)
{
    nmod_poly_one(rop);
}

FQ_NMOD_INLINE void fq_nmod_gen(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    if (ctx->modulus->length == 2)
    {
        nmod_poly_set_coeff_ui(rop, 0,
              nmod_neg(nmod_div(ctx->modulus->coeffs[0],
              ctx->modulus->coeffs[1], ctx->mod), ctx->mod));
    }
    else
    {
        nmod_poly_zero(rop);
        nmod_poly_set_coeff_ui(rop, 0, 0);
        nmod_poly_set_coeff_ui(rop, 1, 1);
    }
}

FLINT_DLL int fq_nmod_get_fmpz(fmpz_t a, const fq_nmod_t b,
                                                      const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_get_nmod_poly(nmod_poly_t a, const fq_nmod_t b,
                                                      const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_set_nmod_poly(fq_nmod_t a, const nmod_poly_t b,
                                                      const fq_nmod_ctx_t ctx);

/* Output ********************************************************************/

#if defined (H_STDIO)               \
  || defined (_H_STDIO)             \
  || defined (_STDIO_H)             \
  || defined (_STDIO_H_)            \
  || defined (__STDIO_H)            \
  || defined (__STDIO_H__)          \
  || defined (_STDIO_INCLUDED)      \
  || defined (__dj_include_stdio_h_)\
  || defined (_FILE_DEFINED)        \
  || defined (__STDIO__)            \
  || defined (_MSL_STDIO_H)         \
  || defined (_STDIO_H_INCLUDED)    \
  || defined (_ISO_STDIO_ISO_H)     \
  || defined (__STDIO_LOADED)       \
  || defined (_STDIO)

#include "nmod_poly.h"

FLINT_DLL int fq_nmod_ctx_fprint(FILE * file, const fq_nmod_ctx_t ctx);

FQ_NMOD_INLINE void fq_nmod_ctx_print(const fq_nmod_ctx_t ctx)
{
    fq_nmod_ctx_fprint(stdout, ctx);
}

FQ_NMOD_INLINE 
int fq_nmod_fprint(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    return nmod_poly_fprint(file, op);
}

FQ_NMOD_INLINE 
void fq_nmod_print(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    nmod_poly_print(op);
}

FQ_NMOD_INLINE 
int fq_nmod_fprint_pretty(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    return nmod_poly_fprint_pretty(file, op, ctx->var);
}

FQ_NMOD_INLINE 
void fq_nmod_print_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    nmod_poly_print_pretty(op, ctx->var);
}

#endif

FLINT_DLL char * fq_nmod_get_str(const fq_nmod_t op, const fq_nmod_ctx_t ctx);

FLINT_DLL char * fq_nmod_get_str_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx);

/* Special functions *********************************************************/

FLINT_DLL void _fq_nmod_trace(fmpz_t rop, const ulong *op, slong len, 
                    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_trace(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_frobenius(ulong *rop, const ulong *op, slong len, slong e, 
                        const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_frobenius(fq_nmod_t rop, const fq_nmod_t op, slong e, const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_norm(fmpz_t rop, const ulong *op, slong len, 
                   const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_norm(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

/* Bit packing ******************************************************/

FLINT_DLL void fq_nmod_bit_pack(fmpz_t f, const fq_nmod_t op, flint_bitcnt_t bit_size,
                 const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bit_unpack(fq_nmod_t rop, const fmpz_t f, flint_bitcnt_t bit_size,
                   const fq_nmod_ctx_t ctx);

/* Miscellaneous *************************************************************/

FLINT_DLL int fq_nmod_next(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx);

FLINT_DLL void fq_nmod_next_not_zero(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx);

/* Inlines *******************************************************************/

FLINT_DLL void __fq_nmod_ctx_prime(fmpz_t p, fq_nmod_ctx_t ctx);

#ifdef T
#undef T
#endif

#define T fq_nmod
#define CAP_T FQ_NMOD
#define B nmod
#include "fq_templates.h"
#undef B
#undef CAP_T
#undef T

#ifdef __cplusplus
}
#endif

#endif
