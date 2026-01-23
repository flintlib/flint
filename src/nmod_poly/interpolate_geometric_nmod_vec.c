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
_nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nn_ptr poly, nn_srcptr v,
    const nmod_geometric_progression_t G, slong len, nmod_t mod)
{
    slong i, N, f1_len, f2_len, h1_len, h2_min;
    nn_ptr f, h;

    N = G->len;

    FLINT_ASSERT(len == N);

    if (len == 1)
    {
        poly[0] = v[0];
        return;
    }

    f = _nmod_vec_init(N);
    h = _nmod_vec_init(N);

    for (i = 0; i < N; i++)
    {
        if (v[N - i - 1] != 0)
        {
            break;
        }
    }

    f1_len = N - i;
    h1_len = FLINT_MIN(G->g1->length + f1_len - 1, N);

    for (i = 0; i < f1_len; i++)
    {
        f[i] = nmod_mul(v[i], G->w[i], mod);
    }
    _nmod_poly_mullow(h, G->g1->coeffs, G->g1->length, f, f1_len, N, mod);

    while (h1_len > 0 && h[h1_len - 1] == 0)
    {
        h1_len--;
    }

    for (i = 0; i < h1_len; i++)
    {
        if (h[i] != 0)
        {
            break;
        }
    }
    h2_min = i;

    for (i = h2_min; i < h1_len; i++)
    {
       f[N - 1 - i] = nmod_mul(h[i], G->y[i], mod);
    }

    f2_len = N - h2_min;
    _nmod_vec_zero(f, N - h1_len);
    _nmod_poly_mullow(h, G->g2->coeffs, G->g2->length, f, f2_len, N, mod);

    for (i = 0; i < len; i++)
    {
        poly[i] = nmod_mul(h[N - 1 - i], G->z[i], mod);
    }

    _nmod_vec_clear(f);
    _nmod_vec_clear(h);
}

void
nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nmod_poly_t poly, nn_srcptr v,
    const nmod_geometric_progression_t G, slong len)
{
    nmod_poly_fit_length(poly, len);
    _nmod_poly_set_length(poly, len);
    _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly->coeffs, v, G, len, G->mod);
    _nmod_poly_normalise(poly);
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
