/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

fmpz ** _fmpz_poly_tree_alloc(slong len)
{
    fmpz ** tree = NULL;

    if (len)
    {
        slong i, height = FLINT_CLOG2(len);

        tree = flint_malloc(sizeof(fmpz_t) * (height + 1));
        for (i = 0; i <= height; i++)
            tree[i] = _fmpz_vec_init(len + (len >> i) + 1);
    }

    return tree;
}

void _fmpz_poly_tree_free(fmpz_t * tree, slong len)
{
    if (len)
    {
        slong i, height = FLINT_CLOG2(len);

        for (i = 0; i <= height; i++)
            flint_free(tree[i]);

        flint_free(tree);
    }
}

/* Given (a1/b1,...,an/bn), construct the subproduct tree of (b1*x-a1)...(bn*x-an) */
void
_fmpz_poly_tree_build_fmpq_vec(fmpz ** tree, const fmpq * roots, slong len)
{
    slong height, pow, left, i;
    fmpz * pa, * pb;
    fmpq_t ab, cd;

    if (len == 0)
        return;

    height = FLINT_CLOG2(len);

    /* zeroth level, (b*x-a) */
    for (i = 0; i < len; i++)
    {
        fmpz_neg(tree[0] + 2 * i,     fmpq_numref(roots + i));
        fmpz_set(tree[0] + 2 * i + 1, fmpq_denref(roots + i));
    }

    /* first level, (b*x-a)(d*x-c) = bd*x^2 + (-ad-bc)*x + ac */
    if (height > 1)
    {
        fmpq_init(ab); fmpq_init(cd);
        pa = tree[1];

        for (i = 0; i < len / 2; i++)
        {
            fmpq_set(ab, roots + 2 * i);
            fmpq_set(cd, roots + 2 * i + 1);

            // deg 0
            fmpz_mul(pa + 3 * i, fmpq_numref(ab), fmpq_numref(cd));
            // deg 1
            fmpz_mul(   pa + 3 * i + 1, fmpq_numref(ab), fmpq_denref(cd));
            fmpz_addmul(pa + 3 * i + 1, fmpq_denref(ab), fmpq_numref(cd));
            fmpz_neg(pa + 3 * i + 1, pa + 3 * i + 1);
            // deg 2
            fmpz_mul(pa + 3 * i + 2, fmpq_denref(ab), fmpq_denref(cd));
        }

        if (len & 1)
        {
            fmpz_neg(pa + 3 * (len / 2),     fmpq_numref(roots + len - 1));
            fmpz_set(pa + 3 * (len / 2) + 1, fmpq_denref(roots + len - 1));
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
            _fmpz_poly_mul(pb, pa, pow + 1, pa + pow + 1, pow + 1);
            left -= 2 * pow;
            pa = pa + 2 * pow + 2;
            pb = pb + 2 * pow + 1;
        }

        if (left > pow)
        {
            _fmpz_poly_mul(pb, pa, pow + 1, pa + pow + 1, left - pow + 1);
        }
        else if (left > 0)
        {
            _fmpz_vec_set(pb, pa, left + 1);
        }
    }
}

/* This gives some speedup for small lengths. */
void
_fmpz_poly_rem_2(gr_ptr q, gr_ptr r, gr_srcptr a, slong al,
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

void
_fmpq_poly_evaluate_vec_fast_precomp(fmpq * vs, const fmpz * poly, fmpz_t den,
    slong plen, const fmpz ** tree, slong len)
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
