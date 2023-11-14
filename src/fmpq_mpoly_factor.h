/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_MPOLY_FACTOR_H
#define FMPQ_MPOLY_FACTOR_H

#ifdef FMPQ_MPOLY_FACTOR_INLINES_C
#define FMPQ_MPOLY_FACTOR_INLINE
#else
#define FMPQ_MPOLY_FACTOR_INLINE static inline
#endif

#include "fmpq_mpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

void fmpq_mpoly_factor_init(fmpq_mpoly_factor_t f,
                                                   const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_factor_realloc(fmpq_mpoly_factor_t f,
                                      slong alloc, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_factor_fit_length(fmpq_mpoly_factor_t f,
                                        slong len, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_factor_clear(fmpq_mpoly_factor_t f,
                                                   const fmpq_mpoly_ctx_t ctx);

FMPQ_MPOLY_FACTOR_INLINE
slong fmpq_mpoly_factor_length(const fmpq_mpoly_factor_t f,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    return f->num;
}

FMPQ_MPOLY_FACTOR_INLINE
void fmpq_mpoly_factor_get_constant_fmpq(fmpq_t c,
                      const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_set(c, f->constant);
}

FMPQ_MPOLY_FACTOR_INLINE
void fmpq_mpoly_factor_get_base(fmpq_mpoly_t p, const fmpq_mpoly_factor_t f,
                                           slong i, const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fmpq_mpoly_set(p, f->poly + i, ctx);
}

FMPQ_MPOLY_FACTOR_INLINE
void fmpq_mpoly_factor_swap_base(fmpq_mpoly_t p, fmpq_mpoly_factor_t f,
                                           slong i, const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    fmpq_mpoly_swap(p, f->poly + i, ctx);
}

FMPQ_MPOLY_FACTOR_INLINE
slong fmpq_mpoly_factor_get_exp_si(fmpq_mpoly_factor_t f,
                                           slong i, const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(i < (ulong) f->num);
    return fmpz_get_si(f->exp + i);
}


void fmpq_mpoly_factor_sort(fmpq_mpoly_factor_t f,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_factor_make_monic(fmpq_mpoly_factor_t f,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_factor_make_integral(fmpq_mpoly_factor_t f,
                                                   const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f,
                             const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_factor(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx);

void _fmpq_mpoly_factor_swap_fmpz_mpoly_factor(fmpq_mpoly_factor_t f,
            fmpz_mpoly_factor_t g, const fmpq_t c, const fmpq_mpoly_ctx_t ctx);

int fmpq_mpoly_factor_expand(fmpq_mpoly_t A,
                      const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

