/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
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

#include "nmod_poly.h"
#include "ulong_extras.h"

/* Data types and context ****************************************************/
#ifdef __cplusplus
extern "C" {
#endif

typedef nmod_polydr_t fq_nmod_t;
typedef nmod_polydr_struct fq_nmod_struct;

typedef struct
{
    fmpz p;
    nmod_ctx_t fpctx;

    int sparse_modulus;

    mp_limb_t *a;
    slong *j;
    slong len;

    nmod_polydr_t modulus;
    nmod_polydr_t inv;

    char *var;
}
fq_nmod_ctx_struct;

typedef fq_nmod_ctx_struct fq_nmod_ctx_t[1];

FLINT_DLL void fq_nmod_ctx_init(fq_nmod_ctx_t ctx,
                      const fmpz_t p, slong d, const char *var);

FLINT_DLL int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx,
                             const fmpz_t p, slong d, const char *var);

FLINT_DLL void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx,
                             const fmpz_t p, slong d, const char *var);

FLINT_DLL void fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx,
                              const nmod_poly_t modulus,
                              const char *var);

FLINT_DLL void fq_nmod_ctx_init_modulusdr(fq_nmod_ctx_t ctx,
                              const nmod_polydr_t modulus, mp_limb_t p,
                              const char *var);

FLINT_DLL void fq_nmod_ctx_randtest(fq_nmod_ctx_t ctx, flint_rand_t state);

FLINT_DLL void fq_nmod_ctx_randtest_reducible(fq_nmod_ctx_t ctx, flint_rand_t state);

FLINT_DLL void fq_nmod_ctx_clear(fq_nmod_ctx_t ctx);

FQ_NMOD_INLINE const nmod_polydr_struct* fq_nmod_ctx_modulus(const fq_nmod_ctx_t ctx)
{
    return ctx->modulus;
}

FQ_NMOD_INLINE slong fq_nmod_ctx_degree(const fq_nmod_ctx_t ctx)
{
    return ctx->modulus->length - 1;
}

#define fq_nmod_ctx_prime(ctx)  (&((ctx)->p))

FQ_NMOD_INLINE void fq_nmod_ctx_order(fmpz_t f, const fq_nmod_ctx_t ctx)
{
    fmpz_set(f, fq_nmod_ctx_prime(ctx));
    fmpz_pow_ui(f, f, fq_nmod_ctx_degree(ctx));
}

/* TODO */
FQ_NMOD_INLINE int fq_nmod_ctx_fprint(FILE * file, const fq_nmod_ctx_t ctx)
{
    int r;
    slong i, k;

    r = flint_fprintf(file, "p = ");
    if (r <= 0)
        return r;

    r = fmpz_fprint(file, fq_nmod_ctx_prime(ctx));
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "\nd = %wd\nf(X) = ", ctx->j[ctx->len - 1]);
    if (r <= 0)
        return r;

    r = flint_fprintf(file, "%wu", ctx->a[0]);
    if (r <= 0)
        return r;

    for (k = 1; k < ctx->len; k++)
    {
        i = ctx->j[k];
        r = flint_fprintf(file, " + ");
        if (r <= 0)
            return r;

        if (ctx->a[k] == UWORD(1))
        {
            if (i == 1)
                r = flint_fprintf(file, "X");
            else
                r = flint_fprintf(file, "X^%wd", i);
            if (r <= 0)
                return r;
        }
        else
        {
            r = flint_fprintf(file, "%wu", ctx->a[k]);
            if (r <= 0)
                return r;

            if (i == 1)
                r = flint_fprintf(file, "*X");
            else
                r = flint_fprintf(file, "*X^%wd", i);
            if (r <= 0)
                return r;
        }
    }
    r = flint_fprintf(file, "\n");
    return r;
}

FQ_NMOD_INLINE void fq_nmod_ctx_print(const fq_nmod_ctx_t ctx)
{
    fq_nmod_ctx_fprint(stdout, ctx);
}

/* Memory managment  *********************************************************/

FQ_NMOD_INLINE void fq_nmod_init(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_polydr_init2(rop, fq_nmod_ctx_degree(ctx), ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_init2(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_polydr_init2(rop, fq_nmod_ctx_degree(ctx), ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_clear(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_polydr_clear(rop, ctx->fpctx);
}

FQ_NMOD_INLINE 
void _fq_nmod_sparse_reduce(mp_limb_t *R, slong lenR, const fq_nmod_ctx_t ctx)
{
    const slong d = ctx->j[ctx->len - 1];

    NMOD_VEC_NORM(R, lenR);

    if (lenR > d)
    {
        slong i, k;

        for (i = lenR - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
            {
                
                R[ctx->j[k] + i - d] = n_submod(R[ctx->j[k] + i - d],
                                                n_mulmod2_preinv(R[i], ctx->a[k], ctx->fpctx->mod.n, ctx->fpctx->mod.ninv),
                                                ctx->fpctx->mod.n);
            }
            R[i] = UWORD(0);
        }
    }
}

FQ_NMOD_INLINE void _fq_nmod_dense_reduce(mp_limb_t* R, slong lenR, const fq_nmod_ctx_t ctx)
{
    mp_limb_t  *q, *r;

    if (lenR < ctx->modulus->length)
    {
        _nmod_vec_reduce(R, R, lenR, ctx->fpctx->mod);
        return;
    }
    
    q = _nmod_vec_init(lenR - ctx->modulus->length + 1);
    r = _nmod_vec_init(ctx->modulus->length - 1);

    _nmod_poly_divrem_newton_n_preinv(q, r, R, lenR, 
                                      ctx->modulus->coeffs, ctx->modulus->length,
                                      ctx->inv->coeffs, ctx->inv->length,
                                      ctx->fpctx->mod);

    _nmod_vec_set(R, r, ctx->modulus->length - 1);
    _nmod_vec_clear(q);
    _nmod_vec_clear(r);

}

FQ_NMOD_INLINE void _fq_nmod_reduce(mp_limb_t* R, slong lenR, const fq_nmod_ctx_t ctx)
{
    if (ctx->sparse_modulus)
        _fq_nmod_sparse_reduce(R, lenR, ctx);
    else
        _fq_nmod_dense_reduce(R, lenR, ctx);    
}

FQ_NMOD_INLINE void fq_nmod_reduce(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    _fq_nmod_reduce(rop->coeffs, rop->length, ctx);
    rop->length = FLINT_MIN(rop->length, ctx->modulus->length - 1);
    _nmod_polydr_normalise(rop);
}

/* Basic arithmetic **********************************************************/

FLINT_DLL void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_sub(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_sub_one(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_neg(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mul(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mul_fmpz(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mul_si(fq_nmod_t rop, const fq_nmod_t op, slong x, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_mul_ui(fq_nmod_t rop, const fq_nmod_t op, ulong x, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_sqr(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_inv(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_pow(mp_limb_t *rop, const mp_limb_t *op, slong len, const fmpz_t e, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_pow(fq_nmod_t rop, const fq_nmod_t op1, const fmpz_t e, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_pow_ui(fq_nmod_t rop, const fq_nmod_t op1, const ulong e, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_pth_root(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

/* Randomisation *************************************************************/

FLINT_DLL void fq_nmod_randtest(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_randtest_dense(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_randtest_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

/* Comparison ****************************************************************/

FQ_NMOD_INLINE int fq_nmod_equal(const fq_nmod_t op1, const fq_nmod_t op2,
                                    const fq_nmod_ctx_t ctx)
{
    return nmod_polydr_equal(op1, op2, ctx->fpctx);
}

FQ_NMOD_INLINE int fq_nmod_is_zero(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    return nmod_polydr_is_zero(op, ctx->fpctx);
}

FQ_NMOD_INLINE int fq_nmod_is_one(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    return nmod_polydr_is_one(op, ctx->fpctx);
}

/* Assignments and conversions ***********************************************/

FQ_NMOD_INLINE void fq_nmod_set(fq_nmod_t rop, const fq_nmod_t op,
                                   const fq_nmod_ctx_t ctx)
{
    nmod_polydr_set(rop, op, ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_set_fmpz(fq_nmod_t rop, const fmpz_t x, const fq_nmod_ctx_t ctx)
{
    nmod_polydr_set_nmod(rop, fmpz_fdiv_ui(x, ctx->fpctx->mod.n), ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_set_si(fq_nmod_t rop, const slong x, const fq_nmod_ctx_t ctx)
{
    mp_limb_t rx = x < 0 ? -x : x;
    rx =  n_mod2_preinv(rx, ctx->fpctx->mod.n, ctx->fpctx->mod.ninv);
    if (x < 0)
        rx = ctx->fpctx->mod.n - rx;

    nmod_polydr_set_nmod(rop, rx, ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_set_ui(fq_nmod_t rop, const ulong x, const fq_nmod_ctx_t ctx)
{
    nmod_polydr_set_nmod(rop, n_mod2_preinv(x, ctx->fpctx->mod.n, ctx->fpctx->mod.ninv), ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_swap(fq_nmod_t op1, fq_nmod_t op2,
                                    const fq_nmod_ctx_t ctx)
{
    nmod_polydr_swap(op1, op2, ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_zero(fq_nmod_t rop,  const fq_nmod_ctx_t ctx)
{
    nmod_polydr_zero(rop, ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_one(fq_nmod_t rop,  const fq_nmod_ctx_t ctx)
{
    nmod_polydr_one(rop, ctx->fpctx);
}

FQ_NMOD_INLINE void fq_nmod_gen(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_polydr_zero(rop, ctx->fpctx);
    if (ctx->modulus->length == 2)
    {
        mp_limb_t t;
        t = nmod_div(ctx->modulus->coeffs[0], ctx->modulus->coeffs[1], ctx->fpctx->mod);
        t = nmod_neg(t, ctx->fpctx->mod);
        nmod_polydr_set_coeff_ui(rop, 0, t, ctx->fpctx);
    }
    else
    {
        nmod_polydr_set_coeff_ui(rop, 0, 0, ctx->fpctx);
        nmod_polydr_set_coeff_ui(rop, 1, 1, ctx->fpctx);
    }
}

/* Output ********************************************************************/

FQ_NMOD_INLINE 
int fq_nmod_fprint(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    return nmod_polydr_fprint(file, op, ctx->fpctx);
}

FQ_NMOD_INLINE 
void fq_nmod_print(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    nmod_polydr_print(op, ctx->fpctx);
}

FQ_NMOD_INLINE 
int fq_nmod_fprint_pretty(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    return nmod_polydr_fprint_pretty(file, op, ctx->var, ctx->fpctx);
}

FQ_NMOD_INLINE 
void fq_nmod_print_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    nmod_polydr_print_pretty(op, ctx->var, ctx->fpctx);
}

FLINT_DLL char * fq_nmod_get_str(const fq_nmod_t op, const fq_nmod_ctx_t ctx);

FLINT_DLL char * fq_nmod_get_str_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx);

/* Special functions *********************************************************/

FLINT_DLL void _fq_nmod_trace(fmpz_t rop, const mp_limb_t *op, slong len, 
                    const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_trace(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_frobenius(mp_limb_t *rop, const mp_limb_t *op, slong len, slong e, 
                        const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_frobenius(fq_nmod_t rop, const fq_nmod_t op, slong e, const fq_nmod_ctx_t ctx);

FLINT_DLL void _fq_nmod_norm(fmpz_t rop, const mp_limb_t *op, slong len, 
                   const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_norm(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

/* Bit packing ******************************************************/

FLINT_DLL void fq_nmod_bit_pack(fmpz_t f, const fq_nmod_t op, mp_bitcnt_t bit_size,
                 const fq_nmod_ctx_t ctx);

FLINT_DLL void fq_nmod_bit_unpack(fq_nmod_t rop, const fmpz_t f, mp_bitcnt_t bit_size,
                   const fq_nmod_ctx_t ctx);

#ifdef T
#undef T
#endif

#define T fq_nmod
#define CAP_T FQ_NMOD
#include "fq_templates.h"
#undef CAP_T
#undef T

#ifdef __cplusplus
}
#endif

#endif

