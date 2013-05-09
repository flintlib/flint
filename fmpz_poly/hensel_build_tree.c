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

    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_hensel_build_tree(len_t *link, fmpz_poly_t *v, fmpz_poly_t *w, 
                                 const nmod_poly_factor_t fac)
{
    const len_t r = fac->num;
    const nmod_t mod = (fac->p + 0)->mod;

    len_t i, j;

    nmod_poly_t d;
    nmod_poly_t *V = flint_malloc((2*r - 2)*sizeof(nmod_poly_t));
    nmod_poly_t *W = flint_malloc((2*r - 2)*sizeof(nmod_poly_t));

    nmod_poly_init_preinv(d, mod.n, mod.ninv);

    for (i = 0; i < 2*r - 2; i++)
    {
        nmod_poly_init_preinv(V[i], mod.n, mod.ninv);
        nmod_poly_init_preinv(W[i], mod.n, mod.ninv);
    }

    for (i = 0; i < r; i++)
    {
        nmod_poly_set(V[i], fac->p + i);
        link[i] = - i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        len_t s;
        len_t minp, mind;
        len_t tmp;

        minp = j;
        mind = nmod_poly_degree(V[j]);
        for (s = j+1; s < i; s++)
        {
            if (nmod_poly_degree(V[s]) < mind)
            {
                minp = s;
                mind = nmod_poly_degree(V[s]);
            }
        }

        nmod_poly_swap(V[j], V[minp]);

        /* Swap link[j] and V[minp] */
        tmp = link[j];
        link[j] = link[minp];
        link[minp] = tmp; 

        minp = j+1;
        mind = nmod_poly_degree(V[j+1]);

        for (s = j + 2; s < i; s++)
        {
            if (nmod_poly_degree(V[s]) < mind)
            {
                minp = s;
                mind = nmod_poly_degree(V[s]);
            }
        }

        nmod_poly_swap(V[j + 1], V[minp]);

        /* Swap link[j+1] and V[minp] */
        tmp = link[j + 1];
        link[j+1] = link[minp];
        link[minp] = tmp;

        nmod_poly_mul(V[i], V[j], V[j+1]);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        /* N.B.  d == 1 */
        nmod_poly_xgcd(d, W[j], W[j+1], V[j], V[j+1]);
    }

    for (j = 0; j < 2*r - 2; j++)
    {
        fmpz_poly_set_nmod_poly(v[j], V[j]);
        fmpz_poly_set_nmod_poly(w[j], W[j]);
    }

    for (i = 0; i < 2*r - 2; i++)
    {
        nmod_poly_clear(V[i]);
        nmod_poly_clear(W[i]);
    }

    nmod_poly_clear(d);
    flint_free(V);
    flint_free(W);
}

