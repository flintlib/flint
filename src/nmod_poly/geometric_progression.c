/*
    Copyright (C) 2025, Vincent Neiger, Éric Schost
    Copyright (C) 2025, Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"

void 
nmod_geometric_progression_init(nmod_geometric_progression_t G, ulong r, slong d, nmod_t mod)
{
    ulong q, inv_r, inv_q, tmp, qk, inv_qk, qq, s;
    nn_ptr diff, inv_diff, prod_diff;
    slong i;
    
    G->d = d;
    G->mod = mod;

    nmod_poly_init2(G->f, mod.n, 2*d-1);
    nmod_poly_init(G->g1, mod.n);
    nmod_poly_init(G->g2, mod.n);
    nmod_poly_set_coeff_ui(G->g1, 0, 1);
    nmod_poly_set_coeff_ui(G->g2, 0, 1);
    
    G->x = _nmod_vec_init(d);
    G->w = _nmod_vec_init(d);
    G->z = _nmod_vec_init(d);
    G->y = _nmod_vec_init(d);

    G->x[0] = 1;
    G->y[0] = 1;
    G->w[0] = 1;
    G->z[0] = 1;

    q = nmod_mul(r, r, mod);
    inv_r = nmod_inv(r, mod);
    inv_q = nmod_mul(inv_r, inv_r, mod);

    nmod_poly_set_coeff_ui(G->f, 0, 1);
    tmp = r;
    for (i = 1; i < 2*d-1; i++)
    {
        nmod_poly_set_coeff_ui(G->f, i, nmod_mul(nmod_poly_get_coeff_ui(G->f, i-1), tmp, mod));
        tmp = nmod_mul(tmp, q, mod);
    }

    tmp = inv_r;
    for (i = 1; i < d; i++)
    {
        G->x[i] = nmod_mul(G->x[i-1], tmp, mod);
        tmp = nmod_mul(tmp, inv_q, mod);
    }
    
    inv_diff = _nmod_vec_init(d);
    diff = _nmod_vec_init(d);
    prod_diff = _nmod_vec_init(d);
    inv_diff[0] = 1;
    diff[0] = 1;
    prod_diff[0] = 1;
 
    qk = q;  // montgomery inversion
    for (i = 1; i < d; i++)
    {
        diff[i] = qk - 1;
        inv_diff[i] = diff[i];
        qk = nmod_mul(qk, q, mod);
        prod_diff[i] = nmod_mul(diff[i], prod_diff[i-1], mod);
    }
    
    tmp = nmod_inv(prod_diff[d-1], mod);
    for (i = d-1; i > 0; i--)
    {
        inv_diff[i] = nmod_mul(prod_diff[i-1], tmp, mod);
        tmp = nmod_mul(tmp, diff[i], mod);
    }
    inv_diff[0] = tmp;
    // end montgomery inversion

    // sets sequences w, y, z and polynomials g1, g2
    qk = 1;
    inv_qk = 1;
    qq = 1;
    s = 1;

    for (i = 1; i < d; i++)
    {
        qq = nmod_mul(qq, qk, mod);   // prod q^i
        s = nmod_mul(s, inv_qk, mod); // prod 1/q^i
        G->w[i] = nmod_mul(G->w[i-1], inv_diff[i], mod); // prod 1/(q^i-1)
        tmp = nmod_mul(qq, G->w[i], mod); // prod q^i/(q^i-1)
        nmod_poly_set_coeff_ui(G->g2, i, tmp);

        if ((i & 1) == 1)
        {
            nmod_poly_set_coeff_ui(G->g1, i, nmod_neg(tmp, mod));
            G->y[i] = nmod_neg(prod_diff[i], mod);
            G->z[i] = nmod_neg(G->w[i], mod);
        }
        else
        {
            nmod_poly_set_coeff_ui(G->g1, i, tmp);
            G->y[i] = prod_diff[i];
            G->z[i] = G->w[i];
        }
        G->y[i] = nmod_mul(G->y[i], s, mod);

        qk = nmod_mul(qk, q, mod);
        inv_qk = nmod_mul(inv_qk, inv_q, mod);
    }

    _nmod_vec_clear(prod_diff);
    _nmod_vec_clear(inv_diff);
    _nmod_vec_clear(diff);
}

void 
nmod_geometric_progression_clear(nmod_geometric_progression_t G)
{
    nmod_poly_clear(G->f);
    nmod_poly_clear(G->g2);
    nmod_poly_clear(G->g1);
    _nmod_vec_clear(G->x);
    _nmod_vec_clear(G->z);
    _nmod_vec_clear(G->y);
    _nmod_vec_clear(G->w);
}