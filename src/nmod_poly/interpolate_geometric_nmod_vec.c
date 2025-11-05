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
#include "nmod_poly.h"
#include "nmod_vec.h"

void
nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nmod_poly_t poly, nn_srcptr v,
    const nmod_geometric_progression_t G, slong len)
{
    slong i, N;
    nmod_poly_t f, h;
    nmod_t mod;

    N = G->len;

    if (len > N) 
    {
        flint_abort();
    }

    nmod_poly_fit_length(poly, len);
    poly->length = len;
    
    if (len == 1)
    {
        poly->coeffs[0] = v[0];
        _nmod_poly_normalise(poly);
        return;
    }

    mod = G->mod;
    nmod_poly_init2(f, mod.n, N); f->length = N;
    nmod_poly_init2(h, mod.n, N); //h->length = N;
    
    for(i = 0; i < N; i++)
    {
        f->coeffs[i] = nmod_mul(v[i], G->w[i], mod);
    }
    _nmod_poly_normalise(f);
    nmod_poly_mullow(h, f, G->g1, N);

    f->length = N;
    for (i = 0; i < h->length; i++)
    {
       f->coeffs[N - 1 - i] = nmod_mul(h->coeffs[i], G->y[i], mod);
    }
    _nmod_vec_zero(f->coeffs, N - h->length);
    _nmod_poly_normalise(f);
    nmod_poly_mullow(h, f, G->g2, N);

    for (i = 0; i < len; i++)
    {
        poly->coeffs[i] = nmod_mul(nmod_poly_get_coeff_ui(h, N - 1 - i), G->z[i], mod);
    }
    _nmod_poly_normalise(poly);

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
