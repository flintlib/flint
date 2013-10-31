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
fq_poly_randtest_irreducible(fq_poly_t f, flint_rand_t state,
                             long len, const fq_ctx_t ctx)
{
    fq_poly_t xq, xqi, x, g_i, finv;
    fmpz_t q;
    slong i, restart;

    /* Compute q */
    fmpz_init_set(q, fq_ctx_prime(ctx));
    fmpz_pow_ui(q, q, fq_ctx_degree(ctx));

    fq_poly_init(x, ctx);
    fq_poly_gen(x, ctx);
    fq_poly_init(xq, ctx);
    fq_poly_init(xqi, ctx);
    fq_poly_init(g_i, ctx);
    fq_poly_init(finv, ctx);

    while (1)
    {
        restart = 0;

        /* Generate random monic polynomial of length len */
        fq_poly_randtest_monic(f, state, len, ctx);

        fq_poly_reverse(finv, f, f->length, ctx);
        fq_poly_inv_series_newton(finv, finv, f->length, ctx);

        /* Compute xq = x^q mod f */
        fq_poly_powmod_fmpz_binexp_preinv(xq, x, q, f, finv, ctx);
        fq_poly_set(xqi, xq, ctx);

        for (i = 1; i <= (len - 1) / 2; i++)
        {
            fq_poly_sub(xqi, xqi, x, ctx);
            fq_poly_gcd(g_i, xqi, f, ctx);
            fq_poly_add(xqi, xqi, x, ctx);
            if (!fq_poly_is_one(g_i, ctx))
            {
                restart = 1;
                break;
            }
            fq_poly_compose_mod_brent_kung_preinv(xqi, xqi, xq, f, finv, ctx);

        }
        if (!restart)
        {
            break;
        }
    }

    fq_poly_clear(x, ctx);
    fq_poly_clear(xq, ctx);
    fq_poly_clear(xqi, ctx);
    fq_poly_clear(g_i, ctx);
    fq_poly_clear(finv, ctx);
    fmpz_clear(q);
}
