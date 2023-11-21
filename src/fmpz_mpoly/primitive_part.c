/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void
fmpz_mpoly_primitive_part(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)
{
    if (res != f)
        fmpz_mpoly_set(res, f, ctx);

    if (fmpz_mpoly_is_zero(res, ctx))
        return;

    if (fmpz_sgn(res->coeffs) < 0)
        fmpz_mpoly_neg(res, res, ctx);

    if (!fmpz_is_one(res->coeffs))
    {
        fmpz_t c;
        fmpz_init(c);

        _fmpz_vec_content(c, res->coeffs, res->length);
        if (!fmpz_is_one(c))
            fmpz_mpoly_scalar_divexact_fmpz(res, res, c, ctx);

        fmpz_clear(c);
    }
}
