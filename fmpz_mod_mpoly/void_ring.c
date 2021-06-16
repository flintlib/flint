/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

static void fmpz_mod_mpoly_void_init(void * a, const void * ctx)
{
    fmpz_mod_mpoly_init(a, ctx);
}

static void fmpz_mod_mpoly_void_clear(void * a, const void * ctx)
{
    fmpz_mod_mpoly_clear(a, ctx);
}

static int fmpz_mod_mpoly_void_is_zero(const void * a, const void * ctx)
{
    return fmpz_mod_mpoly_is_zero(a, ctx);
}

static void fmpz_mod_mpoly_void_zero(void * a, const void * ctx)
{
    fmpz_mod_mpoly_zero(a, ctx);
}

static void fmpz_mod_mpoly_void_one(void * a, const void * ctx)
{
    fmpz_mod_mpoly_one(a, ctx);
}

static void fmpz_mod_mpoly_void_set(void * a, const void * b, const void * ctx)
{
    fmpz_mod_mpoly_set(a, b, ctx);
}

static void fmpz_mod_mpoly_void_set_fmpz(void * a, const fmpz_t b, const void * ctx)
{
    fmpz_mod_mpoly_set_fmpz(a, b, ctx);
}

static void fmpz_mod_mpoly_void_swap(void * a, void * b, const void * ctx)
{
    fmpz_mod_mpoly_swap(a, b, ctx);
}

static void fmpz_mod_mpoly_void_neg(void * a, const void * b, const void * ctx)
{
    fmpz_mod_mpoly_neg(a, b, ctx);
}

static void fmpz_mod_mpoly_void_add(void * a, const void * b, const void * c,
                                                              const void * ctx)
{
    fmpz_mod_mpoly_add(a, b, c, ctx);
}

static void fmpz_mod_mpoly_void_sub(void * a, const void * b, const void * c,
                                                              const void * ctx)
{
    fmpz_mod_mpoly_sub(a, b, c, ctx);
}

static void fmpz_mod_mpoly_void_mul(void * a, const void * b, const void * c,
                                                              const void * ctx)
{
    fmpz_mod_mpoly_mul(a, b, c, ctx);
}

static void fmpz_mod_mpoly_void_mul_fmpz(void * a, const void * b,
                                              const fmpz_t c, const void * ctx)
{
    fmpz_mod_mpoly_scalar_mul_fmpz(a, b, c, ctx);
}

static void fmpz_mod_mpoly_void_divexact(void * a, const void * b,
                                             const void * c, const void * ctx)
{
    if (!fmpz_mod_mpoly_divides(a, b, c, ctx))
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_void_divexact: nonexact");
}

static int fmpz_mod_mpoly_void_divides(void * a, const void * b,
                                              const void * c, const void * ctx)
{
    return fmpz_mod_mpoly_divides(a, b, c, ctx);
}

static int fmpz_mod_mpoly_void_pow_fmpz(void * a, const void * b,
                                              const fmpz_t c, const void * ctx)
{
    return fmpz_mod_mpoly_pow_fmpz(a, b, c, ctx);
}

static slong fmpz_mod_mpoly_void_length(const void * a, const void * ctx)
{
    return fmpz_mod_mpoly_length(a, ctx);
}

void mpoly_void_ring_init_fmpz_mod_mpoly_ctx(
    mpoly_void_ring_t R,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    R->elem_size = sizeof(fmpz_mod_mpoly_struct);
    R->ctx = ctx;
    R->init = fmpz_mod_mpoly_void_init;
    R->clear = fmpz_mod_mpoly_void_clear;
    R->is_zero = fmpz_mod_mpoly_void_is_zero;
    R->zero = fmpz_mod_mpoly_void_zero;
    R->one = fmpz_mod_mpoly_void_one;
    R->set = fmpz_mod_mpoly_void_set;
    R->set_fmpz = fmpz_mod_mpoly_void_set_fmpz;
    R->swap = fmpz_mod_mpoly_void_swap;
    R->neg = fmpz_mod_mpoly_void_neg;
    R->add = fmpz_mod_mpoly_void_add;
    R->sub = fmpz_mod_mpoly_void_sub;
    R->mul = fmpz_mod_mpoly_void_mul;
    R->mul_fmpz = fmpz_mod_mpoly_void_mul_fmpz;
    R->divexact = fmpz_mod_mpoly_void_divexact;
    R->divides = fmpz_mod_mpoly_void_divides;
    R->pow_fmpz = fmpz_mod_mpoly_void_pow_fmpz;
    R->length = fmpz_mod_mpoly_void_length;
}

