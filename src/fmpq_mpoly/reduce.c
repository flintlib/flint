/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

/*
    We are supposed to produce a canonical A assuming that A->zpoly is itself
    canonical. The code should not break if A->zpoly is not canonical.
*/
void fmpq_mpoly_reduce(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t g;

    if (A->zpoly->length < 1)
    {
        fmpq_zero(A->content);
        return;
    }
    else if (fmpq_is_zero(A->content))
    {
        fmpz_mpoly_zero(A->zpoly, ctx->zctx);
        return;
    }

    /* A is nonzero at this point (assuming A->zpoly is canonical) */

    fmpz_init(g);
    _fmpz_vec_content(g, A->zpoly->coeffs, A->zpoly->length);

    if (fmpz_sgn(A->zpoly->coeffs + 0) < 0)
        fmpz_neg(g, g);

    if (fmpz_is_zero(g))
    {
        /* bail if A->zpoly has only zeros stored (A->zpoly not canonical) */
        fmpq_one(A->content);
    }
    else if (fmpz_is_pm1(g))
    {
        if (!fmpz_is_one(g))
        {
            fmpq_neg(A->content, A->content);
            _fmpz_vec_neg(A->zpoly->coeffs, A->zpoly->coeffs, A->zpoly->length);
        }
    }
    else
    {
        fmpq_mul_fmpz(A->content, A->content, g);
        _fmpz_vec_scalar_divexact_fmpz(A->zpoly->coeffs, A->zpoly->coeffs,
                                                          A->zpoly->length, g);
    }

    fmpz_clear(g);
}
