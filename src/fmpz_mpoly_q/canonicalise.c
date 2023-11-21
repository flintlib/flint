/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

void
fmpz_mpoly_q_canonicalise(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(res), ctx))
    {
        return;
    }
    else if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(res), ctx))
    {
        fmpz_mpoly_one(fmpz_mpoly_q_denref(res), ctx);
        return;
    }
    else if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(res), ctx))
    {
        fmpz_t g;
        fmpz_init(g);

        _fmpz_vec_content(g, fmpz_mpoly_q_numref(res)->coeffs, fmpz_mpoly_q_numref(res)->length);
        fmpz_gcd(g, g, fmpz_mpoly_q_denref(res)->coeffs);

        if (fmpz_sgn(fmpz_mpoly_q_denref(res)->coeffs) < 0)
            fmpz_neg(g, g);

        if (!fmpz_is_one(g))
        {
            fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), g, ctx);
            fmpz_divexact(fmpz_mpoly_q_denref(res)->coeffs, fmpz_mpoly_q_denref(res)->coeffs, g);
        }

        fmpz_clear(g);
    }
    else if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(res), ctx))
    {
        fmpz_t g;
        fmpz_init(g);

        _fmpz_vec_content(g, fmpz_mpoly_q_denref(res)->coeffs, fmpz_mpoly_q_denref(res)->length);
        fmpz_gcd(g, g, fmpz_mpoly_q_numref(res)->coeffs);

        if (fmpz_sgn(fmpz_mpoly_q_denref(res)->coeffs) < 0)
            fmpz_neg(g, g);

        if (!fmpz_is_one(g))
        {
            fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_q_denref(res), fmpz_mpoly_q_denref(res), g, ctx);
            fmpz_divexact(fmpz_mpoly_q_numref(res)->coeffs, fmpz_mpoly_q_numref(res)->coeffs, g);
        }

        fmpz_clear(g);
    }
    else
    {
        fmpz_mpoly_t g;
        fmpz_mpoly_init(g, ctx);

        fmpz_mpoly_gcd_assert_successful(g, fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res), ctx);

        if (fmpz_sgn(fmpz_mpoly_q_denref(res)->coeffs) < 0)
            fmpz_mpoly_neg(g, g, ctx);

        if (!fmpz_mpoly_is_one(g, ctx))
        {
            _fmpz_mpoly_q_mpoly_divexact(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), g, ctx);
            _fmpz_mpoly_q_mpoly_divexact(fmpz_mpoly_q_denref(res), fmpz_mpoly_q_denref(res), g, ctx);
        }

        fmpz_mpoly_clear(g, ctx);
    }
}

