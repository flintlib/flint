/*
    Copyright (C) 2018 Daniel Schultz

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
#include "assert.h"


int fmpz_mpoly_discriminant(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    fmpz_mpoly_t lcfx;
    fmpz_mpoly_univar_t rx, fx, fxp;
    fmpz_mpoly_init(lcfx, ctx);
    fmpz_mpoly_univar_init(rx, ctx);
    fmpz_mpoly_univar_init(fx, ctx);
    fmpz_mpoly_univar_init(fxp, ctx);

    success = success && fmpz_mpoly_to_univar(fx, poly2, var, ctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_univar_derivative(fxp, fx, ctx);

    /* the discriminant of a constant polynomial "a" should be "1/a^2" */
    if (fxp->length == 0)
    {
        if (   fmpz_mpoly_equal_si(poly2, WORD(1), ctx)
            || fmpz_mpoly_equal_si(poly2, -WORD(1), ctx)
           )
        {
            fmpz_mpoly_set_si(poly1, WORD(1), ctx);
        }
        else
        {
            flint_throw(FLINT_IMPINV, "Non-unit constant polynomial in fmpz_mpoly_discriminant");
        }

    /* the discriminant of a linear polynomial "a*x+b" should be "1" */
    } else if (fxp->exps[0] == 0)
    {
        fmpz_mpoly_set_ui(poly1, 1, ctx);

    /* the discriminant is (-1)^(n*(n-1)/2) res(f, f')/a_n */
    } else
    {
        if (fx->exps[0] & 2)
            fmpz_mpoly_neg(lcfx, fx->coeffs + 0, ctx);
        else
            fmpz_mpoly_set(lcfx, fx->coeffs + 0, ctx);

        _fmpz_mpoly_univar_pgcd_ducos(rx, fx, fxp, ctx);
        assert (rx->length != 0);
        if (rx->length == 1 && rx->exps[0] == 0)
            fmpz_mpoly_divides_monagan_pearce(poly1, rx->coeffs + 0, lcfx, ctx);
        else
            fmpz_mpoly_zero(poly1, ctx);
    }

cleanup:

    fmpz_mpoly_clear(lcfx, ctx);
    fmpz_mpoly_univar_clear(rx, ctx);
    fmpz_mpoly_univar_clear(fx, ctx);
    fmpz_mpoly_univar_clear(fxp, ctx);
    return success;
}

