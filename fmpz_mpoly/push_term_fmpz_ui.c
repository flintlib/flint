/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    emplaceterm clears c
*/
void _fmpz_mpoly_emplacebackterm_fmpz_ui(fmpz_mpoly_t A,
                       fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    mp_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_ui(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(A, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fmpz_mpoly_fit_length(A, old_length + 1, ctx);
    fmpz_swap(A->coeffs + old_length, c);
    A->length = old_length + 1;
    fmpz_clear(c);
    mpoly_set_monomial_ui(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}


void fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_t A,
                 const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init_set(C, c);
    _fmpz_mpoly_emplacebackterm_fmpz_ui(A, C, exp, ctx);    
}

void fmpz_mpoly_push_term_ui_ui(fmpz_mpoly_t A,
                        ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init_set_ui(C, c);
    _fmpz_mpoly_emplacebackterm_fmpz_ui(A, C, exp, ctx);    
}

void fmpz_mpoly_push_term_si_ui(fmpz_mpoly_t A,
                        slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init(C);
    fmpz_set_si(C, c);
    _fmpz_mpoly_emplacebackterm_fmpz_ui(A, C, exp, ctx);    
}

