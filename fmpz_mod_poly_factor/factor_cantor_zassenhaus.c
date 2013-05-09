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
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include "fmpz_mod_poly_factor.h"

void
fmpz_mod_poly_factor_cantor_zassenhaus(fmpz_mod_poly_factor_t res,
                                       const fmpz_mod_poly_t f)
{
    fmpz_mod_poly_t h, v, g, x;
    len_t i, j, num;

    fmpz_mod_poly_init(h, &f->p);
    fmpz_mod_poly_init(g, &f->p);
    fmpz_mod_poly_init(v, &f->p);
    fmpz_mod_poly_init(x, &f->p);

    fmpz_mod_poly_set_coeff_ui(h, 1, 1);
    fmpz_mod_poly_set_coeff_ui(x, 1, 1);

    fmpz_mod_poly_make_monic(v, f);

    i = 0;
    do
    {
        i++;
        fmpz_mod_poly_powmod_fmpz_binexp(h, h, &f->p, v);

        fmpz_mod_poly_sub(h, h, x);
        fmpz_mod_poly_gcd(g, h, v);
        fmpz_mod_poly_add(h, h, x);

        if (g->length != 1)
        {
            fmpz_mod_poly_make_monic(g, g);
            num = res->num;

            fmpz_mod_poly_factor_equal_deg(res, g, i);
            for (j = num; j < res->num; j++)
                res->exp[j] = fmpz_mod_poly_remove(v, res->poly + j);
        }
    }
    while (v->length >= 2 * i + 3);

    if (v->length > 1)
        fmpz_mod_poly_factor_insert(res, v, 1);

    fmpz_mod_poly_clear(g);
    fmpz_mod_poly_clear(h);
    fmpz_mod_poly_clear(v);
    fmpz_mod_poly_clear(x);
}
