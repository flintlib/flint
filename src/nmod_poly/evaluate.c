/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011, 2012 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fmpz.h"

/* TODO: Add a divconquer method */

void
_nmod_poly_evaluate_fmpz(fmpz_t rop, mp_srcptr poly, const slong len, const fmpz_t c)
{
    fmpz_t t;
    slong m;

    if (len == 0)
    {
        fmpz_zero(rop);
        return;
    }

    if (len == 1 || c == 0)
    {
        fmpz_set_ui(rop, poly[0]);
        return;
    }

    m = len - 1;


    fmpz_init(t);
    fmpz_set_ui(rop, poly[m]);
    m--;

    for ( ; m >= 0; m--)
    {
        fmpz_mul(t, rop, c);
        fmpz_add_ui(rop, t, poly[m]);
    }
    fmpz_clear(t);
}

void
nmod_poly_evaluate_fmpz(fmpz_t rop, const nmod_poly_t poly, const fmpz_t c)
{
    _nmod_poly_evaluate_fmpz(rop, poly->coeffs, poly->length, c);
}

void
nmod_mat_one_addmul(nmod_mat_t dest, const nmod_mat_t mat, mp_limb_t c)
{
    slong i, j;

    if (dest == mat)
    {
        for (i = 0; i < mat->r; i++)
        {
            nmod_mat_entry(dest, i, i) =
                n_addmod(nmod_mat_entry(mat, i, i), c, mat->mod.n);
        }
        return;
    }
    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
        {
            nmod_mat_entry(dest, i, j) = nmod_mat_entry(mat, i, j);
            if (i == j)
            {
                nmod_mat_entry(dest, i, i) =
                    n_addmod(nmod_mat_entry(dest, i, i), c, mat->mod.n);
            }
        }
}

void
_nmod_poly_evaluate_mat_horner(nmod_mat_t dest, mp_srcptr poly, slong len, const nmod_mat_t c)
{
    slong m = len-1;
    nmod_mat_t temp;

    nmod_mat_zero(dest);

    if (len == 0)
    {
        return;
    }

    if (len == 1 || nmod_mat_is_zero(c))
    {
        nmod_mat_one_addmul(dest, dest, poly[0]);
        return;
    }

    nmod_mat_init_set(temp, c);
    nmod_mat_one_addmul(dest, dest, poly[m]);

    for( m-- ; m >= 0 ; m--)
    {
        nmod_mat_mul(temp, dest, c);
        nmod_mat_one_addmul(dest, temp, poly[m]);
    }
    nmod_mat_clear(temp);
}

void
nmod_poly_evaluate_mat_horner(nmod_mat_t dest, const nmod_poly_t poly, const nmod_mat_t c)
{
    nmod_mat_t temp;
    if (dest == c)
    {
        nmod_mat_init_set(temp, c);
        _nmod_poly_evaluate_mat_horner(dest, poly->coeffs, poly->length, temp);
        nmod_mat_clear(temp);
    }
    else
    {
        _nmod_poly_evaluate_mat_horner(dest, poly->coeffs, poly->length, c);
    }
}

void
nmod_poly_evaluate_mat_paterson_stockmeyer(nmod_mat_t dest, const nmod_poly_t poly,
                       const nmod_mat_t c)
{
    slong lim = n_sqrt(poly->length), i, j, rem, quo, curr;
    nmod_mat_t tmat, *temp;

    nmod_mat_zero(dest);

    if (poly->length == 0)
    {
        return;
    }

    if (poly->length == 1 || nmod_mat_is_zero(c))
    {
        nmod_mat_one_addmul(dest, dest, poly->coeffs[0]);
        return;
    }

    temp = flint_malloc((lim + 1) * sizeof(nmod_mat_t));
    nmod_mat_init(temp[0], c->r, c->c, c->mod.n);
    nmod_mat_one(temp[0]);
    nmod_mat_init(temp[1], c->r, c->c, c->mod.n);
    nmod_mat_set(temp[1], c);
    nmod_mat_init(tmat, c->r, c->c, c->mod.n);

    for (i = 2; i <= lim; i++)
    {
        nmod_mat_init(temp[i], c->r, c->c, c->mod.n);
        nmod_mat_mul(temp[i], temp[i - 1], c);
    }

    rem = (poly->length % lim);
    quo = poly->length / lim;
    curr = poly->length - rem - 1;

    for (i = 0; i < rem; i++)
    {
        nmod_mat_scalar_addmul_ui(dest, dest,
                                temp[i], poly->coeffs[poly->length - rem + i]);
    }

    for (i = 0; i < quo; i++)
    {
        nmod_mat_mul(tmat, dest, temp[lim]);
        nmod_mat_scalar_addmul_ui(dest, tmat, temp[lim - 1], poly->coeffs[curr--]);

        for (j = 1; j < lim; j++)
        {
            nmod_mat_scalar_addmul_ui(dest, dest, temp[lim - 1 - j],
                                                          poly->coeffs[curr--]);
        }
    }

    for (i = 0; i <= lim; i++)
    {
        nmod_mat_clear(temp[i]);
    }
    nmod_mat_clear(tmat);
    flint_free(temp);
}

mp_limb_t
_nmod_poly_evaluate_nmod(mp_srcptr poly, slong len, mp_limb_t c, nmod_t mod)
{
    slong m;
    mp_limb_t val;

    if (len == 0)
        return 0;

    if (len == 1 || c == 0)
        return poly[0];

    m = len - 1;

    val = poly[m];
    m--;

    for ( ; m >= 0; m--)
    {
        val = n_mulmod2_preinv(val, c, mod.n, mod.ninv);
        val = n_addmod(val, poly[m], mod.n);
    }

    return val;
}

mp_limb_t
nmod_poly_evaluate_nmod(const nmod_poly_t poly, mp_limb_t c)
{
    return _nmod_poly_evaluate_nmod(poly->coeffs, poly->length, c, poly->mod);
}

void
_nmod_poly_evaluate_nmod_vec(mp_ptr ys, mp_srcptr coeffs, slong len,
    mp_srcptr xs, slong n, nmod_t mod)
{
    if (len < 32)
        _nmod_poly_evaluate_nmod_vec_iter(ys, coeffs, len, xs, n, mod);
    else
        _nmod_poly_evaluate_nmod_vec_fast(ys, coeffs, len, xs, n, mod);
}

void
nmod_poly_evaluate_nmod_vec(mp_ptr ys,
        const nmod_poly_t poly, mp_srcptr xs, slong n)
{
    _nmod_poly_evaluate_nmod_vec(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}

/* This gives some speedup for small lengths. */
static __inline__ void _nmod_poly_rem_2(mp_ptr r, mp_srcptr a, slong al,
    mp_srcptr b, slong bl, nmod_t mod)
{
    if (al == 2)
        r[0] = nmod_sub(a[0], nmod_mul(a[1], b[0], mod), mod);
    else
        _nmod_poly_rem(r, a, al, b, bl, mod);
}

void
_nmod_poly_evaluate_nmod_vec_fast_precomp(mp_ptr vs, mp_srcptr poly,
    slong plen, const mp_ptr * tree, slong len, nmod_t mod)
{
    slong height, i, j, pow, left;
    slong tree_height;
    slong tlen;
    mp_ptr t, u, swap, pa, pb, pc;

    /* avoid worrying about some degenerate cases */
    if (len < 2 || plen < 2)
    {
        if (len == 1)
            vs[0] = _nmod_poly_evaluate_nmod(poly, plen,
                nmod_neg(tree[0][0], mod), mod);
        else if (len != 0 && plen == 0)
            _nmod_vec_zero(vs, len);
        else if (len != 0 && plen == 1)
            for (i = 0; i < len; i++)
                vs[i] = poly[0];
        return;
    }

    t = _nmod_vec_init(len);
    u = _nmod_vec_init(len);

    left = len;

    /* Initial reduction. We allow the polynomial to be larger
       or smaller than the number of points. */
    height = FLINT_BIT_COUNT(plen - 1) - 1;
    tree_height = FLINT_CLOG2(len);
    while (height >= tree_height)
        height--;
    pow = WORD(1) << height;

    for (i = j = 0; i < len; i += pow, j += (pow + 1))
    {
        tlen = ((i + pow) <= len) ? pow : len % pow;
        _nmod_poly_rem(t + i, poly, plen, tree[height] + j, tlen + 1, mod);
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
            _nmod_poly_rem_2(pc, pb, 2 * pow, pa, pow + 1, mod);
            _nmod_poly_rem_2(pc + pow, pb, 2 * pow, pa + pow + 1, pow + 1, mod);

            pa += 2 * pow + 2;
            pb += 2 * pow;
            pc += 2 * pow;
            left -= 2 * pow;
        }

        if (left > pow)
        {
            _nmod_poly_rem(pc, pb, left, pa, pow + 1, mod);
            _nmod_poly_rem(pc + pow, pb, left, pa + pow + 1, left - pow + 1, mod);
        }
        else if (left > 0)
            _nmod_vec_set(pc, pb, left);

        swap = t;
        t = u;
        u = swap;
    }

    _nmod_vec_set(vs, t, len);
    _nmod_vec_clear(t);
    _nmod_vec_clear(u);
}

void _nmod_poly_evaluate_nmod_vec_fast(mp_ptr ys, mp_srcptr poly, slong plen,
    mp_srcptr xs, slong n, nmod_t mod)
{
    mp_ptr * tree;

    tree = _nmod_poly_tree_alloc(n);
    _nmod_poly_tree_build(tree, xs, n, mod);
    _nmod_poly_evaluate_nmod_vec_fast_precomp(ys, poly, plen, tree, n, mod);
    _nmod_poly_tree_free(tree, n);
}

void
nmod_poly_evaluate_nmod_vec_fast(mp_ptr ys,
        const nmod_poly_t poly, mp_srcptr xs, slong n)
{
    _nmod_poly_evaluate_nmod_vec_fast(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}

void
_nmod_poly_evaluate_nmod_vec_iter(mp_ptr ys, mp_srcptr coeffs, slong len,
    mp_srcptr xs, slong n, nmod_t mod)
{
    slong i;
    for (i = 0; i < n; i++)
        ys[i] = _nmod_poly_evaluate_nmod(coeffs, len, xs[i], mod);
}

void
nmod_poly_evaluate_nmod_vec_iter(mp_ptr ys,
        const nmod_poly_t poly, mp_srcptr xs, slong n)
{
    _nmod_poly_evaluate_nmod_vec_iter(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}
