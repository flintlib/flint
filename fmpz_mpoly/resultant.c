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


int fmpz_mpoly_resultant(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
               const fmpz_mpoly_t poly3, slong var, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    int change_sign = 0;
    fmpz_mpoly_univar_t rx, fx, gx;
    fmpz_mpoly_univar_init(rx, ctx);
    fmpz_mpoly_univar_init(fx, ctx);
    fmpz_mpoly_univar_init(gx, ctx);

    success = success && fmpz_mpoly_to_univar(fx, poly2, var, ctx);
    success = success && fmpz_mpoly_to_univar(gx, poly3, var, ctx);
    if (!success)
        goto cleanup;

    if (fx->length == 0 || gx->length == 0)
    {
        fmpz_mpoly_zero(poly1, ctx);
    } else 
    {
        if (fx->exps[0] < gx->exps[0])
        {
            fmpz_mpoly_univar_swap(fx, gx, ctx);
            change_sign = 1 & fx->exps[0] & gx->exps[0];
        }

        if (gx->exps[0] == 0)
        {
            fmpz_mpoly_pow_fps(poly1, gx->coeffs + 0, fx->exps[0], ctx);
        } else {

            _fmpz_mpoly_univar_pgcd_ducos(rx, fx, gx, ctx);
            assert (rx->length != 0);
            if (rx->length == 1 && rx->exps[0] == 0)
                fmpz_mpoly_swap(poly1, rx->coeffs + 0, ctx);
            else
                fmpz_mpoly_zero(poly1, ctx);
        }

        if (change_sign)
        {
            fmpz_mpoly_neg(poly1, poly1, ctx);
        }
    }

cleanup:

    fmpz_mpoly_univar_clear(rx, ctx);
    fmpz_mpoly_univar_clear(fx, ctx);
    fmpz_mpoly_univar_clear(gx, ctx);
    return success;
}

