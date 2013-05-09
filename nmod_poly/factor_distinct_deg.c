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

******************************************************************************/

#include <math.h>
#include "nmod_poly.h"

void nmod_poly_factor_distinct_deg(nmod_poly_factor_t res,
                                   const nmod_poly_t poly, len_t **degs)
{
    nmod_poly_t f, g, s, v, tmp;
    nmod_poly_t *h, *H, *I;
    len_t i, j, l, m, n, index;
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
    nmod_poly_init_preinv(tmp, poly->mod.n, poly->mod.ninv);

    if (!(h = flint_malloc((2 * m + l + 1) * sizeof(nmod_poly_struct))))
    {
        printf("Exception (nmod_poly_factor_distinct_deg):\n");
        printf("Not enough memory.\n");
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

    /* compute baby steps: h[i]=x^{p^i}mod v */
    nmod_poly_set_coeff_ui(h[0], 1, 1);
    for (i = 1; i < l + 1; i++)
        nmod_poly_powmod_ui_binexp(h[i], h[i - 1], poly->mod.n, v);

    /* compute giant steps: H[i]=x^{p^(li)}mod v */
    nmod_poly_set(H[0], h[l]);
    for (j = 1; j < m; j++)
        nmod_poly_compose_mod(H[j], H[j - 1], H[0], v);

    /* compute interval polynomials I[j] = (H_j-h_0)*...*(H_j-h_{l-1})*/
    for (j = 0; j < m; j++)
    {
        nmod_poly_set_coeff_ui(I[j], 0, 1);
        for (i = 0; i < l; i++)
        {
            nmod_poly_sub(tmp, H[j], h[i]);
            nmod_poly_mulmod(I[j], tmp, I[j], v);
        }
    }

    /* compute coarse distinct-degree factorisation */
    index = 0;
    nmod_poly_set(s, v);
    for (j = 0; j < m; j++)
    {
        /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
        /* F_j is stored on the place of I_j */
        nmod_poly_gcd(I[j], s, I[j]);
        if (I[j]->length > 1)
            nmod_poly_remove(s, I[j]);
    }
    if (s->length > 1)
    {
        nmod_poly_factor_insert(res, s, 1);
        (*degs)[index++] = s->length - 1;
    }

    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j]->length > 1)
        {
            nmod_poly_set(g, I[j]);
            for (i = l - 1; i >= 0; i--)
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
    }

    /* cleanup */
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(s);
    nmod_poly_clear(v);
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
