/*
    Copyright (C) 2011, 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fmpz_mat.h"
#include "nmod_poly_mat.h"

#define KS_MIN_DIM 10
#define INTERPOLATE_MIN_DIM 80
#define KS_MAX_LENGTH 128

void
nmod_poly_mat_sqr(nmod_poly_mat_t C, const nmod_poly_mat_t A)
{
    slong dim = A->r;

    if (dim < KS_MIN_DIM)
    {
        nmod_poly_mat_sqr_classical(C, A);
    }
    else
    {
        slong Alen;
        mp_limb_t mod = nmod_poly_mat_modulus(A);

        Alen = nmod_poly_mat_max_length(A);

        if ((FLINT_BIT_COUNT(mod) > FLINT_BITS / 4)
            && (dim > INTERPOLATE_MIN_DIM + n_sqrt(Alen))
            && (mod >= 2 * Alen - 1) && n_is_prime(mod))
            nmod_poly_mat_sqr_interpolate(C, A);

        if (Alen > KS_MAX_LENGTH)
            nmod_poly_mat_sqr_classical(C, A);
        else
            nmod_poly_mat_sqr_KS(C, A);
    }
}

static inline void
nmod_poly_sqr(nmod_poly_t y, const nmod_poly_t x)
{
    nmod_poly_mul(y, x, x);
}

#define E nmod_poly_mat_entry

void
nmod_poly_mat_sqr_classical(nmod_poly_mat_t B, const nmod_poly_mat_t A)
{
    slong n = A->r;

    if (n == 0)
        return;

    if (n == 1)
    {
        nmod_poly_sqr(E(B, 0, 0), E(A, 0, 0));
        return;
    }

    if (n == 2)
    {
        nmod_poly_t t, u;

        nmod_poly_init(t, nmod_poly_mat_modulus(A));
        nmod_poly_init(u, nmod_poly_mat_modulus(A));

        nmod_poly_add(t, E(A, 0, 0), E(A, 1, 1));
        nmod_poly_mul(u, E(A, 0, 1), E(A, 1, 0));

        nmod_poly_sqr(E(B, 0, 0), E(A, 0, 0));
        nmod_poly_add(E(B, 0, 0), E(B, 0, 0), u);

        nmod_poly_sqr(E(B, 1, 1), E(A, 1, 1));
        nmod_poly_add(E(B, 1, 1), E(B, 1, 1), u);

        nmod_poly_mul(E(B, 0, 1), E(A, 0, 1), t);
        nmod_poly_mul(E(B, 1, 0), E(A, 1, 0), t);

        nmod_poly_clear(t);
        nmod_poly_clear(u);
        return;
    }

    nmod_poly_mat_mul_classical(B, A, A);
}

void
nmod_poly_mat_sqr_interpolate(nmod_poly_mat_t C, const nmod_poly_mat_t A)
{
    slong i, j, k;
    slong A_len, len;

    nmod_mat_t *C_mod, *A_mod;

    mp_ptr xs;
    mp_ptr tt, uu;
    mp_ptr * tree;
    mp_ptr weights;
    nmod_t mod;

    if (A->c == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    A_len = nmod_poly_mat_max_length(A);

    if (A_len == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    len = 2 * A_len - 1;
    nmod_init(&mod, nmod_poly_mat_modulus(A));

    if (mod.n < len)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_mat_sqr_interpolate): Characteristic is too small.\n");
    }

    xs = _nmod_vec_init(len);
    tt = _nmod_vec_init(len);
    uu = _nmod_vec_init(len);
    weights = _nmod_vec_init(len);

    A_mod = flint_malloc(sizeof(nmod_mat_t) * len);
    C_mod = flint_malloc(sizeof(nmod_mat_t) * len);

    for (i = 0; i < len; i++)
    {
        xs[i] = i;
        nmod_mat_init(A_mod[i], A->r, A->c, mod.n);
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

    /* should be nmod_mat_sqr */
    for (i = 0; i < len; i++)
        nmod_mat_mul(C_mod[i], A_mod[i], A_mod[i]);

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
        nmod_mat_clear(C_mod[i]);
    }

    flint_free(A_mod);
    flint_free(C_mod);

    _nmod_vec_clear(xs);
    _nmod_vec_clear(tt);
    _nmod_vec_clear(uu);
    _nmod_vec_clear(weights);
}

void
nmod_poly_mat_sqr_KS(nmod_poly_mat_t B, const nmod_poly_mat_t A)
{
    slong i, j, n;
    slong A_len;
    flint_bitcnt_t bit_size;
    fmpz_mat_t AA, BB;

    n = A->r;

    if (n == 0)
    {
        nmod_poly_mat_zero(B);
        return;
    }

    A_len = nmod_poly_mat_max_length(A);

    bit_size = 2 * FLINT_BIT_COUNT(nmod_poly_mat_modulus(A));
    bit_size += FLINT_BIT_COUNT(A_len);
    bit_size += FLINT_BIT_COUNT(n);

    fmpz_mat_init(AA, n, n);
    fmpz_mat_init(BB, n, n);

    for (i = 0; i < n; i++)
        for (j = 0; j < A->c; j++)
            nmod_poly_bit_pack(fmpz_mat_entry(AA, i, j),
                               nmod_poly_mat_entry(A, i, j), bit_size);

    /* Should use fmpz_mat_sqr */
    fmpz_mat_mul(BB, AA, AA);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
                nmod_poly_bit_unpack(nmod_poly_mat_entry(B, i, j),
                    fmpz_mat_entry(BB, i, j), bit_size);

    fmpz_mat_clear(AA);
    fmpz_mat_clear(BB);
}
