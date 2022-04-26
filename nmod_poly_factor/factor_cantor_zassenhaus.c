/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010, 2022 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"

void
nmod_poly_factor_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t f)
{
    nmod_poly_t h, v, g, x;
    slong i, j, num;

    nmod_poly_init_mod(h, f->mod);
    nmod_poly_init_mod(g, f->mod);
    nmod_poly_init_mod(v, f->mod);
    nmod_poly_init_mod(x, f->mod);

    nmod_poly_set_coeff_ui(h, 1, 1);
    nmod_poly_set_coeff_ui(x, 1, 1);

    nmod_poly_make_monic(v, f);

    i = 0;
    do
    {
        i++;
        nmod_poly_powmod_ui_binexp(h, h, f->mod.n, v);

        nmod_poly_sub(h, h, x);
        nmod_poly_gcd(g, h, v);
        nmod_poly_add(h, h, x);

        if (g->length != 1)
        {
            nmod_poly_make_monic(g, g);
            num = res->num;
            nmod_poly_factor_equal_deg(res, g, i);

            for (j = num; j < res->num; j++)
                res->exp[j] = nmod_poly_remove(v, res->p + j);
        }
    } while (v->length >= 2*i + 3);

    if (v->length > 1)
        nmod_poly_factor_insert(res, v, 1);

    nmod_poly_clear(g);
    nmod_poly_clear(h);
    nmod_poly_clear(v);
    nmod_poly_clear(x);
}
