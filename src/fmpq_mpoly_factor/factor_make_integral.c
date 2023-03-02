/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"


int fmpq_mpoly_factor_make_integral(fmpq_mpoly_factor_t f,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i;
    fmpq_t t;

    fmpq_init(t);

    for (i = 0; i < f->num; i++)
    {
        success = !fmpq_is_zero(f->poly[i].content) &&
                  fmpq_pow_fmpz(t, f->poly[i].content, f->exp + i);
        if (!success)
            goto cleanup;
        fmpq_mul(f->constant, f->constant, t);
        fmpq_one(f->poly[i].content);
    }

cleanup:

    fmpq_clear(t);

	return success;
}

