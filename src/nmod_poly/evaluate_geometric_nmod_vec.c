/*
    Copyright (C) 2025, Vincent Neiger, Éric Schost
    Copyright (C) 2025, Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(nn_ptr vs, nn_srcptr poly, 
    slong plen, const nmod_geometric_progression_t G, slong len)
{
    nmod_poly_t a, b;
    slong i, d;

    d = G->d;
    FLINT_ASSERT(len <= d);

    if (plen == 0)
    {
        for (i = 0; i < d; i++)
        {
            vs[i] = 0;
        }
        return;
    }

    nmod_poly_init2(a, G->mod.n, d);
    nmod_poly_init(b, G->mod.n);

    for (i = 0; i < plen; i++)
    {
        nmod_poly_set_coeff_ui(a, d - 1 - i, nmod_mul(G->x[i], poly[i], G->mod));
    }

    for (i = plen + 1; i < d; i++)
    {
        nmod_poly_set_coeff_ui(a, d - 1 - i, 0);
    }

    nmod_poly_mul(b, a, G->f);

    for (i = 0; i < len; i++)
    {
        vs[i] = nmod_mul(G->x[i], nmod_poly_get_coeff_ui(b, i+d-1), G->mod);
    }

    nmod_poly_clear(b);
    nmod_poly_clear(a);
}

void 
_nmod_poly_evaluate_geometric_nmod_vec_fast(nn_ptr ys, nn_srcptr poly, 
    slong plen, ulong r, slong n, nmod_t mod)
{
    nmod_geometric_progression_t G;
    nmod_geometric_progression_init(G, r, n, mod);

    _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(ys, poly, plen, G, n);
    nmod_geometric_progression_clear(G);
}

void
nmod_poly_evaluate_geometric_nmod_vec_fast(nn_ptr ys,
        const nmod_poly_t poly, ulong r, slong n)
{
    _nmod_poly_evaluate_geometric_nmod_vec_fast(ys, poly->coeffs,
                                                poly->length, r, n, poly->mod);
}
