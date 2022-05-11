/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_MINI_H
#define FQ_MINI_H

#ifdef FQ_INLINES_C
#define FQ_INLINE FLINT_DLL
#define FQ_TEMPLATES_INLINE FLINT_DLL
#else
#define FQ_INLINE static __inline__
#define FQ_TEMPLATES_INLINE static __inline__
#endif

#include "fmpz_mini.h"
#include "fmpz_poly_mini.h"

/* Define here to avoid inclusion of fmpz_mod.h */
#define FMPZ_MOD_CTX_MODULUS(ctx) (ctx)->n

#ifdef __cplusplus
extern "C" {
#endif

/* memory management **********************************************************/

FLINT_DLL void _fq_sparse_reduce(fmpz *R, slong lenR, const fq_ctx_t ctx);
FLINT_DLL void _fq_dense_reduce(fmpz* R, slong lenR, const fq_ctx_t ctx);
FQ_INLINE void _fq_reduce(fmpz* R, slong lenR, const fq_ctx_t ctx)
{
    if (ctx->sparse_modulus)
        _fq_sparse_reduce(R, lenR, ctx);
    else
        _fq_dense_reduce(R, lenR, ctx);    
}

FQ_INLINE void fq_reduce(fq_t rop, const fq_ctx_t ctx)
{
    _fq_reduce(rop->coeffs, rop->length, ctx);
    rop->length = FLINT_MIN(rop->length, ctx->modulus->length - 1);
    _fmpz_poly_normalise(rop);
}

/* parameters *****************************************************************/

FQ_INLINE const fmpz_mod_poly_struct* fq_ctx_modulus(const fq_ctx_t ctx)
{
    return ctx->modulus;
}

FQ_INLINE slong fq_ctx_degree(const fq_ctx_t ctx)
{
    return ctx->modulus->length - 1;
}

FQ_INLINE const fmpz * fq_ctx_prime(const fq_ctx_t ctx)
{
    return FMPZ_MOD_CTX_MODULUS(ctx->ctxp);
}

FQ_INLINE void fq_ctx_order(fmpz_t f, const fq_ctx_t ctx)
{
    fmpz_pow_ui(f, fq_ctx_prime(ctx), fq_ctx_degree(ctx));
}

/* assignment *****************************************************************/

FQ_INLINE void fq_set(fq_t rop, const fq_t op, const fq_ctx_t ctx)
{
    fmpz_poly_set(rop, op);
}

FQ_INLINE void fq_set_fmpz(fq_t rop, const fmpz_t x, const fq_ctx_t ctx)
{
    fmpz_poly_set_fmpz(rop, x);
    fq_reduce(rop, ctx);
}

FQ_INLINE void fq_zero(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_zero(rop);
}

FQ_INLINE void fq_swap(fq_t op1, fq_t op2, const fq_ctx_t ctx)
{
    fmpz_poly_swap(op1, op2);
}

/* comparison *****************************************************************/

FQ_INLINE int fq_equal(const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
{
    return fmpz_poly_equal(op1, op2);
}

FQ_INLINE int fq_is_zero(const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_is_zero(op);
}

FQ_INLINE int fq_is_one(const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_is_one(op);
}

/* arithmetics ****************************************************************/

FLINT_DLL void fq_inv(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
