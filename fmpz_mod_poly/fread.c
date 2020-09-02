/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdio.h>
#include "flint.h"
#include "fmpz_mod_poly.h"
#include "fmpz.h"

int fmpz_mod_poly_fread(FILE * f, fmpz_mod_poly_t poly, fmpz_mod_ctx_t ctx)
{
    int success = 0;
    slong i, length;
    fmpz_t coeff;

    fmpz_init(coeff);

    poly->length = 0;

    if (flint_fscanf(f, "%wd", &length) != 1)
        goto cleanup;

    if (!fmpz_fread(f, coeff))
        goto cleanup;

    if (fmpz_cmp_ui(coeff, 2) < 0)
        goto cleanup;

    fmpz_mod_ctx_set_modulus(ctx, coeff);

    fmpz_mod_poly_fit_length(poly, length, ctx);

    for (i = 0; i < length; i++)
    {
        if (!fmpz_fread(f, coeff))
            goto cleanup;

        fmpz_mod_poly_set_coeff_fmpz(poly, i, coeff, ctx);
    }

    poly->length = length;
    _fmpz_mod_poly_normalise(poly);

    success = 1;

cleanup:

    fmpz_clear(coeff);

    return success;
}

