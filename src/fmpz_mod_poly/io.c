/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

/* printing *******************************************************************/

int _fmpz_mod_poly_fprint(FILE * file, const fmpz *poly, slong len,
                          const fmpz_t p)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%wd ", len);
    if (r <= 0)
        return r;

    r = fmpz_fprint(file, p);
    if (r <= 0)
        return r;

    if (len == 0)
        return r;

    r = flint_fprintf(file, " ");
    if (r <= 0)
        return r;

    for (i = 0; (r > 0) && (i < len); i++)
    {
        r = flint_fprintf(file, " ");
        if (r <= 0)
            return r;
        r = fmpz_fprint(file, poly + i);
        if (r <= 0)
            return r;
    }

    return r;
}

int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx) { return _fmpz_mod_poly_fprint(file, poly->coeffs, poly->length, fmpz_mod_ctx_modulus(ctx)); }
int _fmpz_mod_poly_print(const fmpz *poly, slong len, const fmpz_t p) { return _fmpz_mod_poly_fprint(stdout, poly, len, p); }
int fmpz_mod_poly_print(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx) { return fmpz_mod_poly_fprint(stdout, poly, ctx); }
int fmpz_mod_poly_fprint_pretty(FILE * file, const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx) { return _fmpz_poly_fprint_pretty(file, poly->coeffs, poly->length, x); }
int fmpz_mod_poly_print_pretty(const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx) { return fmpz_mod_poly_fprint_pretty(stdout, poly, x, ctx); }

/* reading ********************************************************************/

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
