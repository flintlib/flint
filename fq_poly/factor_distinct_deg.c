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
#include <math.h>

void
fq_poly_factor_distinct_deg(fq_poly_factor_t res, const fq_poly_t poly,
                            slong * const *degs, const fq_ctx_t ctx)
{
    fq_poly_t f, g, s, reducedH0, v, vinv, tmp;
    fq_poly_t *h, *H, *I;
    fmpz_t q;
    slong i, j, l, m, n, index, d;
    double beta;

    n = fq_poly_degree(poly, ctx);
    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    fmpz_init(q);
    fq_ctx_order(q, ctx);

    fq_poly_init(f, ctx);
    fq_poly_init(g, ctx);
    fq_poly_init(s, ctx);
    fq_poly_init(reducedH0, ctx);
    fq_poly_init(v, ctx);
    fq_poly_init(vinv, ctx);
    fq_poly_init(tmp, ctx);

    if (!(h = flint_malloc((2 * m + l + 1) * sizeof(fq_poly_struct))))
    {
        printf("Exception (fq_poly_factor_distinct_deg):\n");
        printf("Not enough memory.\n");
        abort();
    }
    H = h + (l + 1);
    I = H + m;
    for (i = 0; i < l + 1; i++)
        fq_poly_init(h[i], ctx);
    for (i = 0; i < m; i++)
    {
        fq_poly_init(H[i], ctx);
        fq_poly_init(I[i], ctx);
    }

    fq_poly_make_monic(v, poly, ctx);

    fq_poly_reverse(vinv, v, v->length, ctx);
    fq_poly_inv_series_newton(vinv, vinv, v->length, ctx);

    /* compute baby steps: h[i]=x^{q^i}mod v */
    /*     h[0] = x */
    /* TODO: Better cutoff */
    if (FQ_POLY_ITERATED_FROBENIUS_CUTOFF(ctx, v->length))
    {
        fq_poly_gen(h[0], ctx);
        for (i = 1; i < l + 1; i++)
            fq_poly_powmod_fmpz_binexp_preinv(h[i], h[i-1], q, v, vinv, ctx);
    }
    else
    {
        fq_poly_gen(h[0], ctx);
        fq_poly_powmod_fmpz_binexp_preinv(h[1], h[0], q, v, vinv, ctx);
        for (i = 2; i < l + 1; i++)
            fq_poly_compose_mod_preinv(h[i], h[i-1], h[1], v, vinv, ctx);
    }

    /* compute coarse distinct-degree factorisation */
    index = 0;
    fq_poly_set(s, v, ctx);
    fq_poly_set(H[0], h[l], ctx);
    fq_poly_set(reducedH0, H[0], ctx);
    d = 1;
    for (j = 0; j < m; j++)
    {
        /* compute giant steps: H[j]=x^{q^(lj)}mod s */
        if (j > 0)
        {
            fq_poly_rem (reducedH0, reducedH0, s, ctx);
            fq_poly_rem(tmp, H[j-1], s, ctx);
            fq_poly_compose_mod_preinv(H[j], tmp, reducedH0, s, vinv, ctx);
        }
        
        /* compute interval polynomials */
        fq_poly_one(I[j], ctx);
        for (i = l-1; (i >= 0) && (2*d <= s->length-1); i--, d++)
        {
            fq_poly_rem(tmp, h[i], s, ctx);
            fq_poly_sub(tmp, H[j], tmp, ctx);
            fq_poly_mulmod_preinv (I[j], tmp, I[j], s, vinv, ctx);
        }

        /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
        /* F_j is stored on the place of I_j */
        fq_poly_gcd(I[j], s, I[j], ctx);
        if (I[j]->length > 1)
        {
            fq_poly_remove(s, I[j], ctx);
            fq_poly_reverse (vinv, s, s->length, ctx);
            fq_poly_inv_series_newton (vinv, vinv, s->length, ctx);
        }
        if (s->length-1 < 2*d)
        {
            break;
        }
    }
    if (s->length > 1)
    {
        fq_poly_factor_insert(res, s, 1, ctx);
        (*degs)[index++] = s->length - 1;
    }


    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j]->length-1 > (j+1)*l || j == 0)
        {
            fq_poly_set(g, I[j], ctx);
            for (i = l - 1; i >= 0 && (g->length > 1); i-- )
            {
                /* compute f^{[l*(j+1)-i]} */
                fq_poly_sub(tmp, H[j], h[i], ctx);
                fq_poly_gcd(f, g, tmp, ctx);
                if (f->length > 1)
                {
                    /* insert f^{[l*(j+1)-i]} into res */
                    fq_poly_make_monic(f, f, ctx);
                    fq_poly_factor_insert(res, f, 1, ctx);
                    (*degs)[index++] = l * (j + 1) - i;

                    fq_poly_remove(g, f, ctx);
                }
            }
        }
        else if (I[j]->length > 1)
        {
            fq_poly_make_monic(I[j], I[j], ctx);
            fq_poly_factor_insert(res, I[j], 1, ctx);
            (*degs)[index++] = I[j]->length-1;
        }
    }

    /* cleanup */
    fmpz_clear(q);
    fq_poly_clear(f, ctx);
    fq_poly_clear(g, ctx);
    fq_poly_clear(s, ctx);
    fq_poly_clear(reducedH0, ctx);
    fq_poly_clear(v, ctx);
    fq_poly_clear(vinv, ctx);
    fq_poly_clear(tmp, ctx);

    for (i = 0; i < l + 1; i++)
        fq_poly_clear(h[i], ctx);
    for (i = 0; i < m; i++)
    {
        fq_poly_clear(H[i], ctx);
        fq_poly_clear(I[i], ctx);
    }
    flint_free(h);
}
