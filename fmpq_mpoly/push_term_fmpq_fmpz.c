/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_push_term_fmpq_fmpz(fmpq_mpoly_t A,
                const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t C;
    fmpz_init_set(fmpq_numref(C), fmpq_numref(c));
    fmpz_init_set(fmpq_denref(C), fmpq_denref(c));
    _fmpq_mpoly_push_rescale(A, C, ctx);
    _fmpz_mpoly_push_exp_pfmpz(A->zpoly, exp, ctx->zctx);
    fmpz_swap(A->zpoly->coeffs + A->zpoly->length - 1, fmpq_numref(C));
    fmpq_clear(C);
}

void fmpq_mpoly_push_term_fmpz_fmpz(fmpq_mpoly_t A,
                const fmpz_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t C;
    fmpz_init_set(fmpq_numref(C), c);
    fmpz_init_set_ui(fmpq_denref(C), UWORD(1));
    _fmpq_mpoly_push_rescale(A, C, ctx);
    _fmpz_mpoly_push_exp_pfmpz(A->zpoly, exp, ctx->zctx);
    fmpz_swap(A->zpoly->coeffs + A->zpoly->length - 1, fmpq_numref(C));
    fmpq_clear(C);
}

void fmpq_mpoly_push_term_ui_fmpz(fmpq_mpoly_t A,
                       ulong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t C;
    fmpz_init_set_ui(fmpq_numref(C), c);
    fmpz_init_set_ui(fmpq_denref(C), UWORD(1));
    _fmpq_mpoly_push_rescale(A, C, ctx);
    _fmpz_mpoly_push_exp_pfmpz(A->zpoly, exp, ctx->zctx);
    fmpz_swap(A->zpoly->coeffs + A->zpoly->length - 1, fmpq_numref(C));
    fmpq_clear(C);
}

void fmpq_mpoly_push_term_si_fmpz(fmpq_mpoly_t A,
                       slong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t C;
    fmpz_init_set_si(fmpq_numref(C), c);
    fmpz_init_set_ui(fmpq_denref(C), UWORD(1));
    _fmpq_mpoly_push_rescale(A, C, ctx);
    _fmpz_mpoly_push_exp_pfmpz(A->zpoly, exp, ctx->zctx);
    fmpz_swap(A->zpoly->coeffs + A->zpoly->length - 1, fmpq_numref(C));
    fmpq_clear(C);
}
