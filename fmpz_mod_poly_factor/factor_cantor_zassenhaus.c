/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res,
                             const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_t h, v, g, x;
    slong i, j, num;

    fmpz_mod_poly_init(h, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(v, ctx);
    fmpz_mod_poly_init(x, ctx);

    fmpz_mod_poly_set_coeff_ui(h, 1, 1, ctx);
    fmpz_mod_poly_set_coeff_ui(x, 1, 1, ctx);

    fmpz_mod_poly_make_monic(v, f, ctx);

    i = 0;
    do
    {
        i++;
        fmpz_mod_poly_powmod_fmpz_binexp(h, h, fmpz_mod_ctx_modulus(ctx), v, ctx);

        fmpz_mod_poly_sub(h, h, x, ctx);
        fmpz_mod_poly_gcd(g, h, v, ctx);
        fmpz_mod_poly_add(h, h, x, ctx);

        if (g->length != 1)
        {
            fmpz_mod_poly_make_monic(g, g, ctx);
            num = res->num;

            fmpz_mod_poly_factor_equal_deg(res, g, i, ctx);
            for (j = num; j < res->num; j++)
                res->exp[j] = fmpz_mod_poly_remove(v, res->poly + j, ctx);
        }
    }
    while (v->length >= 2 * i + 3);

    if (v->length > 1)
        fmpz_mod_poly_factor_insert(res, v, 1, ctx);

    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(h, ctx);
    fmpz_mod_poly_clear(v, ctx);
    fmpz_mod_poly_clear(x, ctx);
}
