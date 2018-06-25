/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"


int fmpz_mpoly_gcd_is_unit(const fmpz_mpoly_t a, const fmpz_mpoly_t b,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int ret = 1;
    slong v, nvars = ctx->minfo->nvars;
    fmpz_t d, ac, bc;
    fmpz_mpoly_t t;

    fmpz_init(d);
    fmpz_init(ac);
    fmpz_init(bc);
    fmpz_mpoly_init(t, ctx);

    if (a->length == 0)
    {
        ret =  fmpz_mpoly_equal_si(b, +WORD(1), ctx)
            || fmpz_mpoly_equal_si(b, -WORD(1), ctx);
        goto done;
    }

    if (b->length == 0)
    {
        ret =  fmpz_mpoly_equal_si(a, +WORD(1), ctx)
            || fmpz_mpoly_equal_si(a, -WORD(1), ctx);
        goto done;
    }

    _fmpz_vec_content(ac, a->coeffs, a->length);
    _fmpz_vec_content(bc, b->coeffs, b->length);
    fmpz_gcd(d, ac, bc);
    if (!fmpz_is_one(d))
    {
        ret = 0;
        goto done;
    }

    for (v = 0; v < nvars; v++)
    {
        fmpz_mpoly_resultant(t, a, b, v, ctx);
        if (fmpz_mpoly_is_zero(t, ctx))
        {
            ret = 0;
            goto done;
        }
    }

done:
    fmpz_clear(d);
    fmpz_clear(ac);
    fmpz_clear(bc);
    fmpz_mpoly_clear(t, ctx);

    return ret;
}
