/*
    Copyright (C) 2011, 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_interpolate_nmod_vec(nn_ptr poly,
                            nn_srcptr xs, nn_srcptr ys, slong n, nmod_t mod)
{
    if (n < 6)
        _nmod_poly_interpolate_nmod_vec_newton(poly, xs, ys, n, mod);
    else if (n < 16)
        _nmod_poly_interpolate_nmod_vec_barycentric(poly, xs, ys, n, mod);
    else
        _nmod_poly_interpolate_nmod_vec_fast(poly, xs, ys, n, mod);
}

void
nmod_poly_interpolate_nmod_vec(nmod_poly_t poly,
                                    nn_srcptr xs, nn_srcptr ys, slong n)
{
    if (n == 0)
    {
        nmod_poly_zero(poly);
    }
    else
    {
        nmod_poly_fit_length(poly, n);
        poly->length = n;
        _nmod_poly_interpolate_nmod_vec(poly->coeffs,
            xs, ys, n, poly->mod);
        _nmod_poly_normalise(poly);
    }
}

void
_nmod_poly_interpolate_nmod_vec_barycentric(nn_ptr poly,
                            nn_srcptr xs, nn_srcptr ys, slong n, nmod_t mod)
{
    nn_ptr P, Q, w;
    slong i, j;

    if (n == 1)
    {
        poly[0] = ys[0];
        return;
    }

    P = _nmod_vec_init(n + 1);
    Q = _nmod_vec_init(n);
    w = _nmod_vec_init(n);

    _nmod_poly_product_roots_nmod_vec(P, xs, n, mod);

    for (i = 0; i < n; i++)
    {
        w[i] = UWORD(1);
        for (j = 0; j < n; j++)
        {
            if (i != j)
                w[i] = nmod_mul(w[i], nmod_sub(xs[i], xs[j], mod), mod);
        }
        w[i] = n_invmod(w[i], mod.n);
    }

    _nmod_vec_zero(poly, n);

    for (i = 0; i < n; i++)
    {
        _nmod_poly_div_root(Q, P, n + 1, xs[i], mod);
        _nmod_vec_scalar_addmul_nmod(poly, Q, n,
            nmod_mul(w[i], ys[i], mod), mod);
    }

    _nmod_vec_clear(P);
    _nmod_vec_clear(Q);
    _nmod_vec_clear(w);
}

void
nmod_poly_interpolate_nmod_vec_barycentric(nmod_poly_t poly,
                                    nn_srcptr xs, nn_srcptr ys, slong n)
{
    if (n == 0)
    {
        nmod_poly_zero(poly);
    }
    else
    {
        nmod_poly_fit_length(poly, n);
        poly->length = n;
        _nmod_poly_interpolate_nmod_vec_barycentric(poly->coeffs,
            xs, ys, n, poly->mod);
        _nmod_poly_normalise(poly);
    }
}

void
_nmod_poly_interpolation_weights(nn_ptr w, const nn_ptr * tree, slong len, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_interpolation_weights(w, (const gr_ptr *) tree, len, ctx));
}

void
_nmod_poly_interpolate_nmod_vec_fast_precomp(nn_ptr poly, nn_srcptr ys,
    const nn_ptr * tree, nn_srcptr weights, slong len, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_interpolate_fast_precomp(poly, ys, (const gr_ptr *) tree, weights, len, ctx));
}

void
_nmod_poly_interpolate_nmod_vec_fast(nn_ptr poly,
                            nn_srcptr xs, nn_srcptr ys, slong len, nmod_t mod)
{
    nn_ptr * tree;
    nn_ptr w;

    tree = _nmod_poly_tree_alloc(len);
    _nmod_poly_tree_build(tree, xs, len, mod);

    w = _nmod_vec_init(len);
    _nmod_poly_interpolation_weights(w, tree, len, mod);

    _nmod_poly_interpolate_nmod_vec_fast_precomp(poly, ys, tree, w, len, mod);

    _nmod_vec_clear(w);
    _nmod_poly_tree_free(tree, len);
}

void
nmod_poly_interpolate_nmod_vec_fast(nmod_poly_t poly,
                                    nn_srcptr xs, nn_srcptr ys, slong n)
{
    if (n == 0)
    {
        nmod_poly_zero(poly);
    }
    else
    {
        nmod_poly_fit_length(poly, n);
        poly->length = n;
        _nmod_poly_interpolate_nmod_vec_fast(poly->coeffs,
            xs, ys, n, poly->mod);
        _nmod_poly_normalise(poly);
    }
}

static void
_interpolate_newton(nn_ptr ys, nn_srcptr xs, slong n, nmod_t mod)
{
    ulong p, q, t;
    slong i, j;

    for (i = 1; i < n; i++)
    {
        t = ys[i - 1];

        for (j = i; j < n; j++)
        {
            p = nmod_sub(ys[j], t, mod);
            q = nmod_sub(xs[j], xs[j - i], mod);
            t = ys[j];
            q = n_invmod(q, mod.n);
            ys[j] = n_mulmod2_preinv(p, q, mod.n, mod.ninv);
        }
    }
}

static void
_newton_to_monomial(nn_ptr ys, nn_srcptr xs, slong n, nmod_t mod)
{
    ulong t;
    slong i, j;

    for (i = n - 2; i >= 0; i--)
    {
        t = ys[i];
        ys[i] = ys[i + 1];

        for (j = i + 1; j < n - 1; j++)
        {
            ys[j] = nmod_sub(ys[j + 1],
                n_mulmod2_preinv(ys[j], xs[i], mod.n, mod.ninv), mod);
        }

        ys[n - 1] = nmod_sub(t,
            n_mulmod2_preinv(ys[n - 1], xs[i], mod.n, mod.ninv), mod);
    }

    _nmod_poly_reverse(ys, ys, n, n);
}

void
_nmod_poly_interpolate_nmod_vec_newton(nn_ptr poly, nn_srcptr xs,
    nn_srcptr ys, slong n, nmod_t mod)
{
    if (n == 1)
    {
        poly[0] = ys[0];
    }
    else
    {
        _nmod_vec_set(poly, ys, n);
        _interpolate_newton(poly, xs, n, mod);
        while (n > 0 && !poly[n-1]) n--;
        _newton_to_monomial(poly, xs, n, mod);
    }
}

void
nmod_poly_interpolate_nmod_vec_newton(nmod_poly_t poly,
                                    nn_srcptr xs, nn_srcptr ys, slong n)
{
    if (n == 0)
    {
        nmod_poly_zero(poly);
    }
    else
    {
        nmod_poly_fit_length(poly, n);
        poly->length = n;
        _nmod_poly_interpolate_nmod_vec_newton(poly->coeffs,
            xs, ys, n, poly->mod);
        _nmod_poly_normalise(poly);
    }
}
