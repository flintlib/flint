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

    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee

******************************************************************************/

#include <math.h>
#include "nmod_poly.h"

void nmod_poly_factor_distinct_deg(nmod_poly_factor_t res,
                                   const nmod_poly_t poly, slong * const *degs)
{
    nmod_poly_t f, g, s, v, vinv, reducedH0, tmp;
    nmod_poly_t *h, *H, *I;
    slong i, j, l, m, n, index, d;
    double beta;

    n = nmod_poly_degree(poly);
    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow (n, beta));
    m = ceil(0.5 * n / l);

    /* initialization */
    nmod_poly_init_preinv(f, poly->mod.n, poly->mod.ninv);
    nmod_poly_init_preinv(g, poly->mod.n, poly->mod.ninv);
    nmod_poly_init_preinv(s, poly->mod.n, poly->mod.ninv);
    nmod_poly_init_preinv(v, poly->mod.n, poly->mod.ninv);
    nmod_poly_init_preinv(vinv, poly->mod.n, poly->mod.ninv);
    nmod_poly_init_preinv(reducedH0, poly->mod.n, poly->mod.ninv);
    nmod_poly_init_preinv(tmp, poly->mod.n, poly->mod.ninv);

    if (!(h = flint_malloc((2 * m + l + 1) * sizeof(nmod_poly_struct))))
    {
        flint_printf("Exception (nmod_poly_factor_distinct_deg):\n");
        flint_printf("Not enough memory.\n");
        abort();
    }
    H = h + (l + 1);
    I = H + m;
    for (i = 0; i < l + 1; i++)
        nmod_poly_init_preinv(h[i], poly->mod.n, poly->mod.ninv);
    for (i = 0; i < m; i++)
    {
        nmod_poly_init_preinv(H[i], poly->mod.n, poly->mod.ninv);
        nmod_poly_init_preinv(I[i], poly->mod.n, poly->mod.ninv);
    }

    nmod_poly_make_monic(v, poly);

    nmod_poly_reverse(vinv, v, v->length);
    nmod_poly_inv_series(vinv, vinv, v->length);
    /* compute baby steps: h[i]=x^{p^i}mod v */
    nmod_poly_set_coeff_ui(h[0], 1, 1);
    for (i = 1; i < l + 1; i++)
        nmod_poly_powmod_ui_binexp_preinv(h[i], h[i - 1], poly->mod.n, v, vinv); /* may be for large l use compose_mod instead */

    /* compute coarse distinct-degree factorisation */
    index = 0;
    nmod_poly_set(s, v);
    nmod_poly_set(H[0], h[l]);
    nmod_poly_set(reducedH0, H[0]);
    d= 1;
    for (j = 0; j < m; j++)
    {
        /* compute giant steps: H[j]=x^{p^(lj)}mod s */
        if (j > 0)
        {
            nmod_poly_rem (reducedH0, reducedH0, s);
            nmod_poly_rem (tmp, H[j-1], s);
            nmod_poly_compose_mod_brent_kung_preinv(H[j], tmp, reducedH0, s, vinv);
        }
        /* compute interval polynomials */
        nmod_poly_set_coeff_ui(I[j], 0, 1);
        for (i = l-1; (i >= 0) && (2*d <= s->length-1); i--, d++)
        {
            nmod_poly_rem(tmp, h[i], s);
            nmod_poly_sub(tmp, H[j], tmp);
            nmod_poly_mulmod_preinv (I[j], tmp, I[j], s, vinv);
        }

        /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
        /* F_j is stored on the place of I_j */
        nmod_poly_gcd(I[j], s, I[j]);
        if (I[j]->length > 1)
        {
            nmod_poly_remove(s, I[j]);
            nmod_poly_reverse (vinv, s, s->length);
            nmod_poly_inv_series (vinv, vinv, s->length);
        }
        if (s->length-1 < 2*d)
        {
            break;
        }
    }
    if (s->length > 1)
    {
        nmod_poly_factor_insert(res, s, 1);
        (*degs)[index++] = s->length - 1;
    }

    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j]->length-1 > (j+1)*l || j == 0)
        {
            nmod_poly_set(g, I[j]);
            for (i = l - 1; i >= 0 && (g->length > 1); i-- )
            {
                /* compute f^{[l*(j+1)-i]} */
                nmod_poly_sub(tmp, H[j], h[i]);
                nmod_poly_gcd(f, g, tmp);
                if (f->length > 1)
                {
                    /* insert f^{[l*(j+1)-i]} into res */
                    nmod_poly_make_monic(f, f);
                    nmod_poly_factor_insert(res, f, 1);
                    (*degs)[index++] = l * (j + 1) - i;

                    nmod_poly_remove(g, f);
                }
            }
        }
        else if (I[j]->length > 1)
        {
            nmod_poly_make_monic(I[j], I[j]);
            nmod_poly_factor_insert(res, I[j], 1);
            (*degs)[index++] = I[j]->length-1;
        }
    }

    /* cleanup */
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(s);
    nmod_poly_clear(reducedH0);
    nmod_poly_clear(v);
    nmod_poly_clear(vinv);
    nmod_poly_clear(tmp);

    for (i = 0; i < l + 1; i++)
        nmod_poly_clear(h[i]);
    for (i = 0; i < m; i++)
    {
        nmod_poly_clear(H[i]);
        nmod_poly_clear(I[i]);
    }
    flint_free(h);
}
