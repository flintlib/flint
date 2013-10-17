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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

void
fq_poly_factor_cantor_zassenhaus(fq_poly_factor_t res, const fq_poly_t f,
                                 const fq_ctx_t ctx)
{
    fq_poly_t h, v, g, x;
    fmpz_t q;
    slong i, j, num;

    fmpz_init(q);
    fq_ctx_order(q, ctx);

    fq_poly_init(h, ctx);
    fq_poly_init(g, ctx);
    fq_poly_init(v, ctx);
    fq_poly_init(x, ctx);

    fq_poly_gen(h, ctx);
    fq_poly_gen(x, ctx);

    fq_poly_make_monic(v, f, ctx);

    i = 0;
    do
    {
        i++;
        fq_poly_powmod_fmpz_binexp(h, h, q, v, ctx);

        fq_poly_sub(h, h, x, ctx);
        fq_poly_gcd(g, h, v, ctx);
        fq_poly_add(h, h, x, ctx);

        if (g->length != 1)
        {
            fq_poly_make_monic(g, g, ctx);
            num = res->num;

            fq_poly_factor_equal_deg(res, g, i, ctx);
            for (j = num; j < res->num; j++)
                res->exp[j] = fq_poly_remove(v, res->poly + j, ctx);
        }
    }
    while (v->length >= 2 * i + 3);

    if (v->length > 1)
        fq_poly_factor_insert(res, v, 1, ctx);

    fq_poly_clear(g, ctx);
    fq_poly_clear(h, ctx);
    fq_poly_clear(v, ctx);
    fq_poly_clear(x, ctx);
    fmpz_clear(q);
}
