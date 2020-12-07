/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

void fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res,
                             const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_poly_t t, h, v, g, x;
    fmpz_mod_poly_factor_t tfac;

    res->num = 0;

    fmpz_mod_poly_init(t, ctx);
    fmpz_mod_poly_init(h, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(v, ctx);
    fmpz_mod_poly_init(x, ctx);
    fmpz_mod_poly_factor_init(tfac, ctx);

    fmpz_mod_poly_gen(h, ctx);
    fmpz_mod_poly_gen(x, ctx);

    fmpz_mod_poly_make_monic(v, f, ctx);

    i = 0;
    do {
        i++;
        fmpz_mod_poly_powmod_fmpz_binexp(t, h, fmpz_mod_ctx_modulus(ctx), v, ctx);
        fmpz_mod_poly_swap(h, t, ctx);

        fmpz_mod_poly_sub(t, h, x, ctx);
        fmpz_mod_poly_gcd(g, t, v, ctx);

        if (g->length != 1)
        {
            FLINT_ASSERT(fmpz_mod_poly_is_monic(g, ctx));

            fmpz_mod_poly_factor_equal_deg(tfac, g, i, ctx);
            fmpz_mod_poly_factor_fit_length(res, res->num + tfac->num, ctx);
            for (j = 0; j < tfac->num; j++)
            {
                res->exp[res->num] = fmpz_mod_poly_remove(v, tfac->poly + j, ctx);
                fmpz_mod_poly_swap(res->poly + res->num, tfac->poly + j, ctx);
                res->num++;
            }
        }
    } while (v->length >= 2 * i + 3);

    if (v->length > 1)
        fmpz_mod_poly_factor_insert(res, v, 1, ctx);

    fmpz_mod_poly_clear(t, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(h, ctx);
    fmpz_mod_poly_clear(v, ctx);
    fmpz_mod_poly_clear(x, ctx);
    fmpz_mod_poly_factor_clear(tfac, ctx);
}

