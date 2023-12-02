/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_mul_interpolate(nmod_poly_mat_t C, const nmod_poly_mat_t A,
    const nmod_poly_mat_t B)
{
    slong i, j, k;
    slong A_len, B_len, len;

    nmod_mat_t *C_mod, *A_mod, *B_mod;

    mp_ptr xs;
    mp_ptr tt, uu;
    mp_ptr * tree;
    mp_ptr weights;
    nmod_t mod;

    if (B->r == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    A_len = nmod_poly_mat_max_length(A);
    B_len = nmod_poly_mat_max_length(B);

    if (A_len == 0 || B_len == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    len = A_len + B_len - 1;
    nmod_init(&mod, nmod_poly_mat_modulus(A));

    if (mod.n < len)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_mat_mul_interpolate): Characteristic is too small.\n");
    }

    xs = _nmod_vec_init(len);
    tt = _nmod_vec_init(len);
    uu = _nmod_vec_init(len);
    weights = _nmod_vec_init(len);

    A_mod = flint_malloc(sizeof(nmod_mat_t) * len);
    B_mod = flint_malloc(sizeof(nmod_mat_t) * len);
    C_mod = flint_malloc(sizeof(nmod_mat_t) * len);

    for (i = 0; i < len; i++)
    {
        xs[i] = i;
        nmod_mat_init(A_mod[i], A->r, A->c, mod.n);
        nmod_mat_init(B_mod[i], B->r, B->c, mod.n);
        nmod_mat_init(C_mod[i], C->r, C->c, mod.n);
    }

    tree = _nmod_poly_tree_alloc(len);
    _nmod_poly_tree_build(tree, xs, len, mod);
    _nmod_poly_interpolation_weights(weights, tree, len, mod);

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            _nmod_poly_evaluate_nmod_vec_fast_precomp(tt,
                nmod_poly_mat_entry(A, i, j)->coeffs,
                nmod_poly_mat_entry(A, i, j)->length,
                tree, len, mod);

            for (k = 0; k < len; k++)
                A_mod[k]->rows[i][j] = tt[k];
        }
    }

    for (i = 0; i < B->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            _nmod_poly_evaluate_nmod_vec_fast_precomp(tt,
                nmod_poly_mat_entry(B, i, j)->coeffs,
                nmod_poly_mat_entry(B, i, j)->length,
                tree, len, mod);

            for (k = 0; k < len; k++)
                B_mod[k]->rows[i][j] = tt[k];
        }
    }

    for (i = 0; i < len; i++)
        nmod_mat_mul(C_mod[i], A_mod[i], B_mod[i]);

    for (i = 0; i < C->r; i++)
    {
        for (j = 0; j < C->c; j++)
        {
            nmod_poly_struct * poly;

            for (k = 0; k < len; k++)
                tt[k] = C_mod[k]->rows[i][j];

            poly = nmod_poly_mat_entry(C, i, j);
            nmod_poly_fit_length(poly, len);
            _nmod_poly_interpolate_nmod_vec_fast_precomp(poly->coeffs,
                tt, tree, weights, len, mod);
            poly->length = len;
            _nmod_poly_normalise(poly);
        }
    }

    _nmod_poly_tree_free(tree, len);

    for (i = 0; i < len; i++)
    {
        nmod_mat_clear(A_mod[i]);
        nmod_mat_clear(B_mod[i]);
        nmod_mat_clear(C_mod[i]);
    }

    flint_free(A_mod);
    flint_free(B_mod);
    flint_free(C_mod);

    _nmod_vec_clear(xs);
    _nmod_vec_clear(tt);
    _nmod_vec_clear(uu);
    _nmod_vec_clear(weights);
}
