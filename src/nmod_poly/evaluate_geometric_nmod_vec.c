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
_nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(nn_ptr vs, nn_srcptr poly, 
    slong plen, const nmod_geometric_progression_t G, slong len, nmod_t mod)
{
    nn_ptr a, b;
    slong i, i_min, G_len, a_len, b_len;

    G_len = G->len;
    FLINT_ASSERT(len <= G_len);

    if (plen == 0)
    {
        _nmod_vec_zero(vs, G_len);
        return;
    }

    for (i_min = 0; i_min < plen; i_min++)
    {
        if (poly[i_min] != 0) 
        {
            break;
        }
    }
    
    a_len = plen - i_min;
    b_len = G->f->length + a_len - 1;
    a = _nmod_vec_init(a_len);
    b = _nmod_vec_init(b_len);
    
    for (i = i_min; i < plen; i++)
    {
        a[plen - 1 - i] = nmod_mul(G->x[i], poly[i], mod);
    }

    // this is a temporary replacement to _nmod_poly_mulhigh which is not yet optimised
    nn_ptr Gfr;
    Gfr = _nmod_vec_init(G->f->length);
    _nmod_poly_reverse(Gfr, G->f->coeffs, G->f->length, G->f->length);
    _nmod_poly_reverse(a, a, a_len, a_len);
    _nmod_poly_mullow(b, Gfr, G->f->length, a, a_len, G->f->length - i_min, mod);
    _nmod_poly_reverse(b, b, b_len, b_len);
    _nmod_vec_clear(Gfr);
    //_nmod_poly_mulhigh(b, G->f->coeffs, G->f->length, a, a_len, plen - 1, mod);
 
    for (i = 0; i < len; i++)
    {
        vs[i] = nmod_mul(G->x[i], b[plen - 1 + i], G->mod);
    }

    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
}

void 
_nmod_poly_evaluate_geometric_nmod_vec_fast(nn_ptr ys, nn_srcptr poly, 
    slong plen, ulong r, slong n, nmod_t mod)
{
    if (n == 0)
        return;

    nmod_geometric_progression_t G;
    nmod_geometric_progression_init(G, r, FLINT_MAX(n, plen), mod);
    _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(ys, poly, plen, G, n, mod);
    nmod_geometric_progression_clear(G);
}

void
nmod_poly_evaluate_geometric_nmod_vec_fast(nn_ptr ys,
        const nmod_poly_t poly, ulong r, slong n)
{
    _nmod_poly_evaluate_geometric_nmod_vec_fast(ys, poly->coeffs,
                                                poly->length, r, n, poly->mod);
}

void
_nmod_poly_evaluate_geometric_nmod_vec_iter(nn_ptr ys, nn_srcptr coeffs, slong len,
    ulong r, slong n, nmod_t mod)
{
    slong i;
    ulong rpow = 1;
    
    ulong r2 = nmod_mul(r, r, mod);

    for (i = 0; i < n; i++)
    {
        ys[i] = _nmod_poly_evaluate_nmod(coeffs, len, rpow, mod);
        rpow = nmod_mul(rpow, r2, mod);
    }
}

void
nmod_poly_evaluate_geometric_nmod_vec_iter(nn_ptr ys,
    const nmod_poly_t poly, ulong r, slong n)
{
    _nmod_poly_evaluate_geometric_nmod_vec_iter(ys, poly->coeffs,
                                        poly->length, r, n, poly->mod);
}
