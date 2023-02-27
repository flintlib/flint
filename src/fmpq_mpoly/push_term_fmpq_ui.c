/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void _fmpq_mpoly_push_rescale(fmpq_mpoly_t A,
                                          fmpq_t C, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_struct * Z = A->zpoly;

    if (!fmpz_is_one(fmpq_numref(A->content)))
    {
        _fmpz_vec_scalar_mul_fmpz(Z->coeffs, Z->coeffs, Z->length,
                                                      fmpq_numref(A->content));
        fmpz_one(fmpq_numref(A->content));
    }

    fmpq_mul_fmpz(C, C, fmpq_denref(A->content));
    if (!fmpz_is_one(fmpq_denref(C)))
    {
        _fmpz_vec_scalar_mul_fmpz(Z->coeffs, Z->coeffs, Z->length,
                                                               fmpq_denref(C));
        fmpz_mul(fmpq_denref(A->content), fmpq_denref(A->content),
                                                               fmpq_denref(C));
    }
}

void fmpq_mpoly_push_term_fmpq_ui(fmpq_mpoly_t A,
                 const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t C;
    fmpz_init_set(fmpq_numref(C), fmpq_numref(c));
    fmpz_init_set(fmpq_denref(C), fmpq_denref(c));
    _fmpq_mpoly_push_rescale(A, C, ctx);
    _fmpz_mpoly_push_exp_ui(A->zpoly, exp, ctx->zctx);
    fmpz_swap(A->zpoly->coeffs + A->zpoly->length - 1, fmpq_numref(C));
    fmpq_clear(C);
}


void fmpq_mpoly_push_term_fmpz_ui(fmpq_mpoly_t A,
                 const fmpz_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t C;
    fmpz_init_set(fmpq_numref(C), c);
    fmpz_init_set_ui(fmpq_denref(C), UWORD(1));
    _fmpq_mpoly_push_rescale(A, C, ctx);
    _fmpz_mpoly_push_exp_ui(A->zpoly, exp, ctx->zctx);
    fmpz_swap(A->zpoly->coeffs + A->zpoly->length - 1, fmpq_numref(C));
    fmpq_clear(C);
}

void fmpq_mpoly_push_term_ui_ui(fmpq_mpoly_t A,
                        ulong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t C;
    fmpz_init_set_ui(fmpq_numref(C), c);
    fmpz_init_set_ui(fmpq_denref(C), UWORD(1));
    _fmpq_mpoly_push_rescale(A, C, ctx);
    _fmpz_mpoly_push_exp_ui(A->zpoly, exp, ctx->zctx);
    fmpz_swap(A->zpoly->coeffs + A->zpoly->length - 1, fmpq_numref(C));
    fmpq_clear(C);
}

void fmpq_mpoly_push_term_si_ui(fmpq_mpoly_t A,
                        slong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t C;
    fmpz_init_set_si(fmpq_numref(C), c);
    fmpz_init_set_ui(fmpq_denref(C), UWORD(1));
    _fmpq_mpoly_push_rescale(A, C, ctx);
    _fmpz_mpoly_push_exp_ui(A->zpoly, exp, ctx->zctx);
    fmpz_swap(A->zpoly->coeffs + A->zpoly->length - 1, fmpq_numref(C));
    fmpq_clear(C);
}
