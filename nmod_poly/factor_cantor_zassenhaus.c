/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include "nmod_poly.h"

void
nmod_poly_factor_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t f)
{
    nmod_poly_t h, v, g, x;
    len_t i, j, num;

    nmod_poly_init_preinv(h, f->mod.n, f->mod.ninv);
    nmod_poly_init_preinv(g, f->mod.n, f->mod.ninv);
    nmod_poly_init_preinv(v, f->mod.n, f->mod.ninv);
    nmod_poly_init_preinv(x, f->mod.n, f->mod.ninv);

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
    }
    while (v->length >= 2*i + 3);

    if (v->length > 1)
        nmod_poly_factor_insert(res, v, 1);

    nmod_poly_clear(g);
    nmod_poly_clear(h);
    nmod_poly_clear(v);
    nmod_poly_clear(x);
}
