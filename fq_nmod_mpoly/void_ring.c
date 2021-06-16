/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

static void fq_nmod_mpoly_void_init(void * a, const void * ctx)
{
    fq_nmod_mpoly_init(a, ctx);
}

static void fq_nmod_mpoly_void_clear(void * a, const void * ctx)
{
    fq_nmod_mpoly_clear(a, ctx);
}

static int fq_nmod_mpoly_void_is_zero(const void * a, const void * ctx)
{
    return fq_nmod_mpoly_is_zero(a, ctx);
}

static void fq_nmod_mpoly_void_zero(void * a, const void * ctx)
{
    fq_nmod_mpoly_zero(a, ctx);
}

static void fq_nmod_mpoly_void_one(void * a, const void * ctx)
{
    fq_nmod_mpoly_one(a, ctx);
}

static void fq_nmod_mpoly_void_set(void * a, const void * b, const void * ctx)
{
    fq_nmod_mpoly_set(a, b, ctx);
}

static void fq_nmod_mpoly_void_set_fmpz(void * a, const fmpz_t b, const void * ctx)
{
    fq_nmod_mpoly_set_fmpz(a, b, ctx);
}

static void fq_nmod_mpoly_void_swap(void * a, void * b, const void * ctx)
{
    fq_nmod_mpoly_swap(a, b, ctx);
}

static void fq_nmod_mpoly_void_neg(void * a, const void * b, const void * ctx)
{
    fq_nmod_mpoly_neg(a, b, ctx);
}

static void fq_nmod_mpoly_void_add(void * a, const void * b, const void * c,
                                                              const void * ctx)
{
    fq_nmod_mpoly_add(a, b, c, ctx);
}

static void fq_nmod_mpoly_void_sub(void * a, const void * b, const void * c,
                                                              const void * ctx)
{
    fq_nmod_mpoly_sub(a, b, c, ctx);
}

static void fq_nmod_mpoly_void_mul(void * a, const void * b, const void * c,
                                                              const void * ctx)
{
    fq_nmod_mpoly_mul(a, b, c, ctx);
}

static void fq_nmod_mpoly_void_mul_fmpz(void * a, const void * b,
                                             const fmpz_t c, const void * ctx_)
{
    const fq_nmod_mpoly_ctx_struct * ctx = ctx_;
    fq_nmod_t C;
    fq_nmod_init(C, ctx->fqctx);
    fq_nmod_set_fmpz(C, c, ctx->fqctx);
    fq_nmod_mpoly_scalar_mul_fq_nmod(a, b, C, ctx);
    fq_nmod_clear(C, ctx->fqctx);
}

static void fq_nmod_mpoly_void_divexact(void * a, const void * b,
                                             const void * c, const void * ctx)
{
    if (!fq_nmod_mpoly_divides(a, b, c, ctx))
        flint_throw(FLINT_ERROR, "fq_nmod_mpoly_void_divexact: nonexact");
}

static int fq_nmod_mpoly_void_divides(void * a, const void * b,
                                              const void * c, const void * ctx)
{
    return fq_nmod_mpoly_divides(a, b, c, ctx);
}

static int fq_nmod_mpoly_void_pow_fmpz(void * a, const void * b,
                                              const fmpz_t c, const void * ctx)
{
    return fq_nmod_mpoly_pow_fmpz(a, b, c, ctx);
}

static slong fq_nmod_mpoly_void_length(const void * a, const void * ctx)
{
    return fq_nmod_mpoly_length(a, ctx);
}

void mpoly_void_ring_init_fq_nmod_mpoly_ctx(
    mpoly_void_ring_t R,
    const fq_nmod_mpoly_ctx_t ctx)
{
    R->elem_size = sizeof(fq_nmod_mpoly_struct);
    R->ctx = ctx;
    R->init = fq_nmod_mpoly_void_init;
    R->clear = fq_nmod_mpoly_void_clear;
    R->is_zero = fq_nmod_mpoly_void_is_zero;
    R->zero = fq_nmod_mpoly_void_zero;
    R->one = fq_nmod_mpoly_void_one;
    R->set = fq_nmod_mpoly_void_set;
    R->set_fmpz = fq_nmod_mpoly_void_set_fmpz;
    R->swap = fq_nmod_mpoly_void_swap;
    R->neg = fq_nmod_mpoly_void_neg;
    R->add = fq_nmod_mpoly_void_add;
    R->sub = fq_nmod_mpoly_void_sub;
    R->mul = fq_nmod_mpoly_void_mul;
    R->mul_fmpz = fq_nmod_mpoly_void_mul_fmpz;
    R->divexact = fq_nmod_mpoly_void_divexact;
    R->divides = fq_nmod_mpoly_void_divides;
    R->pow_fmpz = fq_nmod_mpoly_void_pow_fmpz;
    R->length = fq_nmod_mpoly_void_length;
}

