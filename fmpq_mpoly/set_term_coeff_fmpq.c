/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_set_term_coeff_fmpq(fmpq_mpoly_t A,
                           slong i, const fmpq_t x, const fmpq_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) fmpq_mpoly_length(A, ctx))
    {
        flint_throw(FLINT_ERROR, "index out of range in fmpq_mpoly_set_term_coeff_fmpq");
    }

    if (fmpq_is_zero(x))
    {
        /* we can easily get a zero coeff by zeroing the coeff of the zpoly */
        fmpz_mpoly_set_term_coeff_fmpz(A->zpoly, i, fmpq_numref(x), ctx->zctx);

    }
    else if (fmpq_is_zero(A->content))
    {
        /*
            if the A->content is zero,
                set the nth coeff to 1
                set all other coeffs to 0
                set the content to x
        */
        fmpz_t t;
        fmpz_init_set_ui(t, UWORD(1));
        fmpq_set(A->content, x);
        _fmpz_vec_zero(A->zpoly->coeffs, A->zpoly->length);
        fmpz_mpoly_set_term_coeff_fmpz(A->zpoly, i, t, ctx->zctx);
        fmpz_clear(t);
    }
    else
    {
        /*
            if A->content != 0 and x != 0,
                compute nun/den = x / A->content
                scale zpoly by den
                divide A->content by den
                set the nth coeff of zpoly to num
        */
        fmpq_t t;
        fmpq_init(t);
        fmpq_div(t, x, A->content);
        if (!fmpz_is_one(fmpq_denref(t)))
        {
            fmpq_div_fmpz(A->content, A->content, fmpq_denref(t));
            _fmpz_vec_scalar_mul_fmpz(A->zpoly->coeffs, A->zpoly->coeffs,
                                             A->zpoly->length, fmpq_denref(t));            
        }
        fmpz_mpoly_set_term_coeff_fmpz(A->zpoly, i, fmpq_numref(t), ctx->zctx);
        fmpq_clear(t);
    }
}
