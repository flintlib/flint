/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: multithread */

gr_ptr * _gr_poly_tree_alloc(slong len, gr_ctx_t ctx)
{
    gr_ptr * tree = NULL;

    if (len)
    {
        slong i, height = FLINT_CLOG2(len);

        tree = flint_malloc(sizeof(gr_ptr) * (height + 1));
        for (i = 0; i <= height; i++)
        {
            slong n = len + (len >> i) + 1;

            tree[i] = flint_malloc(n * ctx->sizeof_elem);
            _gr_vec_init(tree[i], n, ctx);
        }
    }

    return tree;
}

void _gr_poly_tree_free(gr_ptr * tree, slong len, gr_ctx_t ctx)
{
    if (len)
    {
        slong i, height = FLINT_CLOG2(len);

        for (i = 0; i <= height; i++)
        {
            slong n = len + (len >> i) + 1;

            _gr_vec_clear(tree[i], n, ctx);
            flint_free(tree[i]);
        }

        flint_free(tree);
    }
}

int
_gr_poly_mul_monic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len1 + len2 - 2 > 0)
        status |= _gr_poly_mullow(res, poly1, len1, poly2, len2, len1 + len2 - 2, ctx);

    status |= gr_one(GR_ENTRY(res, len1 + len2 - 2, ctx->sizeof_elem), ctx);

    return status;
}

int
_gr_poly_tree_build(gr_ptr * tree, gr_srcptr roots, slong len, gr_ctx_t ctx)
{
    slong height, pow, left, i;
    gr_ptr pa, pb;
    gr_srcptr a, b;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len == 0)
        return GR_SUCCESS;

    height = FLINT_CLOG2(len);

    /* zeroth level, (x-a) */
    for (i = 0; i < len; i++)
    {
        status |= gr_one(GR_ENTRY(tree[0], 2 * i + 1, sz), ctx);
        status |= gr_neg(GR_ENTRY(tree[0], 2 * i, sz), GR_ENTRY(roots, i, sz), ctx);
    }

    /* first level, (x-a)(x-b) = x^2 + (-a-b)*x + a*b */
    if (height > 1)
    {
        pa = tree[1];

        for (i = 0; i < len / 2; i++)
        {
            a = GR_ENTRY(roots, 2 * i, sz);
            b = GR_ENTRY(roots, 2 * i + 1, sz);

            status |= gr_mul(GR_ENTRY(pa, (3 * i), sz), a, b, ctx);
            status |= gr_add(GR_ENTRY(pa, (3 * i + 1), sz), a, b, ctx);
            status |= gr_neg(GR_ENTRY(pa, (3 * i + 1), sz), GR_ENTRY(pa, (3 * i + 1), sz), ctx);
            status |= gr_one(GR_ENTRY(pa, (3 * i + 2), sz), ctx);
        }

        if (len & 1)
        {
            status |= gr_neg(GR_ENTRY(pa, 3 * (len / 2), sz), GR_ENTRY(roots, len - 1, sz), ctx);
            status |= gr_one(GR_ENTRY(pa, 3 * (len / 2) + 1, sz), ctx);
        }
    }

    for (i = 1; i < height - 1; i++)
    {
        left = len;
        pow = WORD(1) << i;
        pa = tree[i];
        pb = tree[i + 1];

        while (left >= 2 * pow)
        {
            status |= _gr_poly_mul_monic(pb, pa, pow + 1, GR_ENTRY(pa, pow + 1, sz), pow + 1, ctx);
            left -= 2 * pow;
            pa = GR_ENTRY(pa, 2 * pow + 2, sz);
            pb = GR_ENTRY(pb, 2 * pow + 1, sz);
        }

        if (left > pow)
        {
            status |= _gr_poly_mul_monic(pb, pa, pow + 1, GR_ENTRY(pa, pow + 1, sz), left - pow + 1, ctx);
        }
        else if (left > 0)
        {
            status |= _gr_vec_set(pb, pa, left + 1, ctx);
        }
    }

    return status;
}

/* This gives some speedup for small lengths. */
static inline int
_gr_poly_rem_2(gr_ptr q, gr_ptr r, gr_srcptr a, slong al,
    gr_srcptr b, slong bl, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (al == 2)
    {
        status |= gr_mul(r, GR_ENTRY(a, 1, sz), b, ctx);
        status |= gr_sub(r, a, r, ctx);
    }
    else
    {
        status |= _gr_poly_divrem(q, r, a, al, b, bl, ctx);
    }

    return status;
}

int
_gr_poly_evaluate_vec_fast_precomp(gr_ptr vs, gr_srcptr poly,
    slong plen, gr_ptr * tree, slong len, gr_ctx_t ctx)
{
    slong height, i, j, pow, left;
    slong tree_height;
    slong tlen, alloc, qalloc;
    gr_ptr tmp, q, t, u, swap, pa, pb, pc;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    /* avoid worrying about some degenerate cases */
    if (len < 2 || plen < 2)
    {
        if (len == 1)
        {
            gr_ptr tmp;
            GR_TMP_INIT(tmp, ctx);
            status |= gr_neg(tmp, tree[0], ctx);
            status |= _gr_poly_evaluate(vs, poly, plen, tmp, ctx);
            GR_TMP_CLEAR(tmp, ctx);
        }
        else if (len != 0 && plen == 0)
        {
            status |= _gr_vec_zero(vs, len, ctx);
        }
        else if (len != 0 && plen == 1)
        {
            for (i = 0; i < len; i++)
                status |= gr_set(GR_ENTRY(vs, i, sz), poly, ctx);
        }
        return status;
    }

    left = len;

    /* Initial reduction. We allow the polynomial to be larger
        or smaller than the number of points. */
    height = FLINT_BIT_COUNT(plen - 1) - 1;
    tree_height = FLINT_CLOG2(len);
    while (height >= tree_height)
        height--;
    pow = WORD(1) << height;

    /* Maximum temp space for q (todo: could use dedicated
       rem instead of divrem methods where more efficient). */
    qalloc = pow;

    for (i = j = 0; i < len; i += pow, j += (pow + 1))
    {
        tlen = ((i + pow) <= len) ? pow : len % pow;
        qalloc = FLINT_MAX(qalloc, plen - tlen);
    }

    alloc = 2 * len + qalloc;

    GR_TMP_INIT_VEC(tmp, alloc, ctx);

    t = tmp;
    u = GR_ENTRY(t, len, sz);
    q = GR_ENTRY(u, len, sz);

    for (i = j = 0; i < len; i += pow, j += (pow + 1))
    {
        tlen = ((i + pow) <= len) ? pow : len % pow;
        status |= _gr_poly_divrem(q, GR_ENTRY(t, i, sz), poly, plen, GR_ENTRY(tree[height], j, sz), tlen + 1, ctx);
    }

    for (i = height - 1; i >= 0; i--)
    {
        pow = WORD(1) << i;
        left = len;
        pa = tree[i];
        pb = t;
        pc = u;

        while (left >= 2 * pow)
        {
            status |= _gr_poly_rem_2(q, pc, pb, 2 * pow, pa, pow + 1, ctx);
            status |= _gr_poly_rem_2(q, GR_ENTRY(pc, pow, sz), pb, 2 * pow, GR_ENTRY(pa, pow + 1, sz), pow + 1, ctx);

            pa = GR_ENTRY(pa, 2 * pow + 2, sz);
            pb = GR_ENTRY(pb, 2 * pow, sz);
            pc = GR_ENTRY(pc, 2 * pow, sz);
            left -= 2 * pow;
        }

        if (left > pow)
        {
            status |= _gr_poly_divrem(q, pc, pb, left, pa, pow + 1, ctx);
            status |= _gr_poly_divrem(q, GR_ENTRY(pc, pow, sz), pb, left, GR_ENTRY(pa, pow + 1, sz), left - pow + 1, ctx);
        }
        else if (left > 0)
            status |= _gr_vec_set(pc, pb, left, ctx);

        swap = t;
        t = u;
        u = swap;
    }

    _gr_vec_swap(vs, t, len, ctx);

    GR_TMP_CLEAR_VEC(tmp, alloc, ctx);

    return status;
}

int _gr_poly_evaluate_vec_fast(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)
{
    gr_ptr * tree;
    int status = GR_SUCCESS;

    tree = _gr_poly_tree_alloc(n, ctx);
    status |= _gr_poly_tree_build(tree, xs, n, ctx);
    status |= _gr_poly_evaluate_vec_fast_precomp(ys, poly, plen, tree, n, ctx);
    _gr_poly_tree_free(tree, n, ctx);

    return status;
}

int
gr_poly_evaluate_vec_fast(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)
{
    gr_vec_set_length(ys, xs->length, ctx);
    return _gr_poly_evaluate_vec_fast(ys->entries, poly->coeffs, poly->length, xs->entries, xs->length, ctx);
}
