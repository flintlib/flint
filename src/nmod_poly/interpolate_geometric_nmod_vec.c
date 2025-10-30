/*
    Copyright (C) 2025, Vincent Neiger
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
nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nmod_poly_t poly, nn_srcptr v,
    const nmod_geometric_progression_t G, slong len)
{
    slong i, N;
    nmod_poly_t f, h;
    nmod_t mod;

    N = G->d;
    FLINT_ASSERT(len == N);
    nmod_poly_fit_length(poly, N);
    nmod_poly_zero(poly);
    
    if (N == 1)
    {
        nmod_poly_set_coeff_ui(poly, 0, v[0]);
        return;
    }

    mod = G->mod;
    nmod_poly_init2(f, mod.n, N);
    nmod_poly_init2(h, mod.n, N);
    
    for(i = 0; i < N; i++)
    {
        nmod_poly_set_coeff_ui(f, i, nmod_mul(v[i], G->w[i], mod));
    }

    nmod_poly_mullow(h, f, G->g1, N);

    for (i = 0; i < N; i++)
    {
        nmod_poly_set_coeff_ui(f, N - 1 - i, nmod_mul(nmod_poly_get_coeff_ui(h, i), G->y[i], mod));
    }

    nmod_poly_mullow(h, f, G->g2, N);

    for (i = 0; i < N; i++)
    {
        nmod_poly_set_coeff_ui(poly, i, nmod_mul(nmod_poly_get_coeff_ui(h, N - 1 - i), G->z[i], mod));
    }

    nmod_poly_clear(f);
    nmod_poly_clear(h);
}

void
nmod_poly_interpolate_geometric_nmod_vec_fast(nmod_poly_t poly, ulong r, 
    nn_srcptr ys, slong n)
{
    nmod_geometric_progression_t G;
    nmod_geometric_progression_init(G, r, n, poly->mod);

    nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly, ys, G, n);
    nmod_geometric_progression_clear(G);
}
