/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"


int fmpq_mpoly_factor_make_monic(fmpq_mpoly_factor_t f,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i;
    fmpq_t t;

    fmpq_init(t);

    for (i = 0; i < f->num; i++)
    {
        if (fmpq_is_zero(f->poly[i].content) || f->poly[i].zpoly->length < 1)
        {
            success = 0;
            goto cleanup;
        }

        fmpq_div_fmpz(t, f->poly[i].content, f->poly[i].zpoly->coeffs + 0);

        if (!fmpq_pow_fmpz(t, t, f->exp + i))
        {
            success = 0;
            goto cleanup;
        }

        fmpq_div(f->constant, f->constant, t);

        fmpz_one(fmpq_numref(f->poly[i].content));
        fmpz_set(fmpq_denref(f->poly[i].content), f->poly[i].zpoly->coeffs + 0);
    }

cleanup:

    fmpq_clear(t);

	return success;
}

