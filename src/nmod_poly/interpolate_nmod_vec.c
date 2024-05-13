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
    nn_ptr tmp;
    slong i, n, height;

    if (len == 0)
        return;

    if (len == 1)
    {
        w[0] = 1;
        return;
    }

    tmp = _nmod_vec_init(len + 1);
    height = FLINT_CLOG2(len);
    n = WORD(1) << (height - 1);

    _nmod_poly_mul(tmp, tree[height-1], n + 1,
                        tree[height-1] + (n + 1), (len - n + 1), mod);

    _nmod_poly_derivative(tmp, tmp, len + 1, mod);
    _nmod_poly_evaluate_nmod_vec_fast_precomp(w, tmp, len, tree, len, mod);

    for (i = 0; i < len; i++)
        w[i] = n_invmod(w[i], mod.n);

    _nmod_vec_clear(tmp);
}

void
_nmod_poly_interpolate_nmod_vec_fast_precomp(nn_ptr poly, nn_srcptr ys,
    const nn_ptr * tree, nn_srcptr weights, slong len, nmod_t mod)
{
    nn_ptr t, u, pa, pb;
    slong i, pow, left;

    if (len == 0)
        return;

    t = _nmod_vec_init(len);
    u = _nmod_vec_init(len);

    for (i = 0; i < len; i++)
        poly[i] = nmod_mul(weights[i], ys[i], mod);

    for (i = 0; i < (slong) FLINT_CLOG2(len); i++)
    {
        pow = (WORD(1) << i);
        pa = tree[i];
        pb = poly;
        left = len;

        while (left >= 2 * pow)
        {
            _nmod_poly_mul(t, pa, pow + 1, pb + pow, pow, mod);
            _nmod_poly_mul(u, pa + pow + 1, pow + 1, pb, pow, mod);
            _nmod_vec_add(pb, t, u, 2 * pow, mod);

            left -= 2 * pow;
            pa += 2 * pow + 2;
            pb += 2 * pow;
        }

        if (left > pow)
        {
            _nmod_poly_mul(t, pa, pow + 1, pb + pow, left - pow, mod);
            _nmod_poly_mul(u, pb, pow, pa + pow + 1, left - pow + 1, mod);
            _nmod_vec_add(pb, t, u, left, mod);
        }
    }

    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
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
