/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"

#define MAC(h, m, l, a, b)                          \
{                                                   \
    mp_limb_t p1, p0;                               \
    umul_ppmm(p1, p0, a, b);                        \
    add_sssaaaaaa(h, m, l, h, m, l, 0, p1, p0);     \
}

/********************** p = 1 mod 4 ****************************************/

static slong _find_eval_points4(
    mp_limb_t * list,
    slong d,
    nmod_t ctx)
{
    slong i, len;
    mp_limb_t p = ctx.n;
    mp_limb_t n;

    FLINT_ASSERT(d > 0);
    FLINT_ASSERT((p & UWORD(3)) == 1);

    list[0] = 1;
    len = 1;

    for (n = 2; len < d && n <= (p - 1)/2; n++)
    {
        int ok = 1;
        mp_limb_t mn2 = p - nmod_mul(n, n, ctx);
        for (i = 0; ok && i < len; i++)
            ok = (nmod_mul(list[i], list[i], ctx) != mn2);
        if (ok)
            list[len++] = n;
    }
    return len;
}

static int _fill_matrices4(
    mp_limb_t * M,          /* length d by 4d */
    mp_limb_t * Q,          /* length d by 4d+1 */
    slong d,
    nmod_t ctx)
{
    slong i, j;
    n_poly_t g, h;
    mp_limb_t * list;
    mp_limb_t g0i, c;

    list = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    if (d != _find_eval_points4(list, d, ctx))
    {
        flint_free(list);
        return 0;
    }

    n_poly_init2(g, 4*d + 4);
    n_poly_init2(h, 4*d + 4);

    n_poly_one(g);
    for (i = 0; i < d; i++)
    {
        n_poly_mod_shift_left_scalar_addmul(g, 4,
                             nmod_neg(nmod_pow_ui(list[i], 4, ctx), ctx), ctx);
    }

    g0i = nmod_inv(g->coeffs[0], ctx);
    for (i = 0; i < d; i++)
    {
        FLINT_ASSERT(4*i+4 < g->length);
        Q[i*(4*d+1) + 0] = nmod_mul(g0i, g->coeffs[4*i+4], ctx);
        n_poly_mod_div_root(h, g, list[i], ctx);
        c = n_poly_mod_evaluate_nmod(h, list[i], ctx);
        c = nmod_mul(list[i], c, ctx);
        c = nmod_inv(c, ctx);
        for (j = 0; j < 4*d; j++)
        {
            M[i*(4*d) + j] = nmod_pow_ui(list[i], 1+j, ctx);
            Q[(j/4)*(4*d+1) + 4*i + (j%4) + 1] = nmod_mul(h->coeffs[j], c, ctx);
        }
    }

    n_poly_clear(g);
    n_poly_clear(h);
    flint_free(list);
    return 1;
}


static void _from_coeffs4(
    mp_limb_t * v,          /* length 4d+1 */
    const mp_limb_t * a,
    slong alen,
    const mp_limb_t * M,    /* length d by 4d */
    slong d,
    mp_limb_t w,
    nmod_t ctx)
{
    slong i, j;

    FLINT_ASSERT(0 <= alen);
    FLINT_ASSERT(alen <= 1 + 4*d);

    if (alen <= 1)
    {
        mp_limb_t t = (alen == 1) ? a[0] : 0;
        for (i = 0; i < 4*d+1; i++)
            v[i] = t;
        return;
    }

    v[0] = a[0];
    for (i = 0; i < d; i++)
    {
        mp_limb_t t1, t2, t3, t4;
        mp_limb_t c1h, c1m, c1;
        mp_limb_t c2h, c2m, c2;
        mp_limb_t c3h, c3m, c3;
        mp_limb_t c4h, c4m, c4;
        c1h = c1m = c1 = 0;
        c2h = c2m = c2 = 0;
        c3h = c3m = c3 = 0;
        c4h = c4m = c4 = 0;

        for (j = 0; j + 4 < alen; j += 4)
        {
            FLINT_ASSERT(j + 3 < 4*d);
            MAC(c1h, c1m, c1, a[j + 1], M[j + 0]);
            MAC(c2h, c2m, c2, a[j + 2], M[j + 1]);
            MAC(c3h, c3m, c3, a[j + 3], M[j + 2]);
            MAC(c4h, c4m, c4, a[j + 4], M[j + 3]);
        }

        MAC(c1h, c1m, c1, j + 1 < alen ? a[j + 1] : 0, M[j + 0]);
        MAC(c2h, c2m, c2, j + 2 < alen ? a[j + 2] : 0, M[j + 1]);
        MAC(c3h, c3m, c3, j + 3 < alen ? a[j + 3] : 0, M[j + 2]);
        FLINT_ASSERT(j + 4 >= alen);

        NMOD_RED3(c4, c4h, c4m, c4, ctx);
        NMOD_RED3(c1, c1h, c1m, c1, ctx);
        NMOD_RED3(c2, c2h, c2m, c2, ctx);
        NMOD_RED3(c3, c3h, c3m, c3, ctx);

        M += 4*d;

        c4 = nmod_add(c4, a[0], ctx);
        t1 = nmod_add(c4, c2, ctx);
        t2 = nmod_sub(c4, c2, ctx);
        t3 = nmod_add(c1, c3, ctx);
        t4 = nmod_mul(nmod_sub(c1, c3, ctx), w, ctx);

        v[4*i + 1] = nmod_add(t1, t3, ctx);
        v[4*i + 2] = nmod_add(t2, t4, ctx);
        v[4*i + 3] = nmod_sub(t1, t3, ctx);
        v[4*i + 4] = nmod_sub(t2, t4, ctx);
    }
}

static void _to_coeffs4(
    mp_limb_t * a,          /* length 4d+1 */
    const mp_limb_t * v,    /* length 4d+1 */
    mp_limb_t * t,          /* length 4d   */
    const mp_limb_t * Q,    /* length d by 4d+1 */
    slong d,
    mp_limb_t w,
    nmod_t ctx)
{
    slong i, j;

    a[0] = v[0];

    for (i = 0; i < d; i++)
    {
        mp_limb_t t2 = nmod_add(v[1+4*i+0], v[1+4*i+2], ctx);
        mp_limb_t t1 = nmod_sub(v[1+4*i+0], v[1+4*i+2], ctx);
        mp_limb_t t3 = nmod_add(v[1+4*i+1], v[1+4*i+3], ctx);
        mp_limb_t t4 = nmod_mul(nmod_sub(v[1+4*i+1], v[1+4*i+3], ctx), w, ctx);
        t[4*i+0] = nmod_sub(t1, t4, ctx);
        t[4*i+1] = nmod_sub(t2, t3, ctx);
        t[4*i+2] = nmod_add(t1, t4, ctx);
        t[4*i+3] = nmod_add(t2, t3, ctx);
    }

    for (i = 0; i < d; i++)
    {
        mp_limb_t c1h, c1m, c1;
        mp_limb_t c2h, c2m, c2;
        mp_limb_t c3h, c3m, c3;
        mp_limb_t c4h, c4m, c4;
        c1h = c1m = c1 = 0;
        c2h = c2m = c2 = 0;
        c3h = c3m = c3 = 0;
        c4h = c4m = c4 = 0;
        umul_ppmm(c4m, c4, Q[0], v[0]);
        for (j = 0; j < d; j++)
        {
            MAC(c1h, c1m, c1, t[4*j + 0], Q[4*j + 1]);
            MAC(c2h, c2m, c2, t[4*j + 1], Q[4*j + 2]);
            MAC(c3h, c3m, c3, t[4*j + 2], Q[4*j + 3]);
            MAC(c4h, c4m, c4, t[4*j + 3], Q[4*j + 4]);
        }

        Q += 4*d + 1;

        NMOD_RED3(a[4*i + 1], c1h, c1m, c1, ctx);
        NMOD_RED3(a[4*i + 2], c2h, c2m, c2, ctx);
        NMOD_RED3(a[4*i + 3], c3h, c3m, c3, ctx);
        NMOD_RED3(a[4*i + 4], c4h, c4m, c4, ctx);
    }
}



static int _fill_matrices2(
    mp_limb_t * M,          /* length d by 2d */
    mp_limb_t * Q,          /* length d by 2d+1 */
    slong d,
    nmod_t ctx)
{
    slong i, j;
    n_poly_t g, h;
    mp_limb_t g0i, c;

    if (2*d >= ctx.n)
        return 0;

    n_poly_init2(g, 2*d + 2);
    n_poly_init2(h, 2*d + 2);

    n_poly_one(g);
    for (i = 0; i < d; i++)
    {
        n_poly_mod_shift_left_scalar_addmul(g, 2,
                             nmod_neg(nmod_pow_ui(i + 1, 2, ctx), ctx), ctx);
    }

    g0i = nmod_inv(g->coeffs[0], ctx);
    for (i = 0; i < d; i++)
    {
        FLINT_ASSERT(2*(i+1) < g->length);
        Q[i*(2*d+1) + 0] = nmod_mul(g0i, g->coeffs[2*(i+1)], ctx);
        n_poly_mod_div_root(h, g, i + 1, ctx);
        c = n_poly_mod_evaluate_nmod(h, i + 1, ctx);
        c = nmod_mul(i + 1, c, ctx);
        c = nmod_inv(c, ctx);
        for (j = 0; j < 2*d; j++)
        {
            M[i*(2*d) + j] = nmod_pow_ui(i + 1, 1+j, ctx);
            Q[(j/2)*(2*d+1) + 2*i + (j%2) + 1] = nmod_mul(h->coeffs[j], c, ctx);
        }
    }

    n_poly_clear(g);
    n_poly_clear(h);
    return 1;
}


static void _from_coeffs2(
    mp_limb_t * v,          /* length 2d+1 */
    const mp_limb_t * a,    /* length alen <= 2d+1 */
    slong alen,
    const mp_limb_t * M,    /* length d by 2d */
    slong d,
    nmod_t ctx)
{
    slong i, j;

    FLINT_ASSERT(0 <= alen);
    FLINT_ASSERT(alen <= 1 + 2*d);

    if (alen <= 1)
    {
        mp_limb_t t = (alen == 1) ? a[0] : 0;
        for (i = 0; i < 2*d+1; i++)
            v[i] = t;
        return;
    }

    v[0] = a[0];
    for (i = 0; i < d; i++)
    {
        mp_limb_t c1h, c1m, c1;
        mp_limb_t c2h, c2m, c2;
        c1h = c1m = c1 = 0;
        c2h = c2m = c2 = 0;

        for (j = 0; j + 2 < alen; j += 2)
        {
            FLINT_ASSERT(j + 1 < 4*d);
            MAC(c1h, c1m, c1, a[j + 1], M[j + 0]);
            MAC(c2h, c2m, c2, a[j + 2], M[j + 1]);
        }

        MAC(c1h, c1m, c1, j + 1 < alen ? a[j + 1] : 0, M[j + 0]);
        FLINT_ASSERT(j + 2 >= alen);

        NMOD_RED3(c2, c2h, c2m, c2, ctx);
        NMOD_RED3(c1, c1h, c1m, c1, ctx);

        M += 2*d;

        c2 = nmod_add(c2, a[0], ctx);

        v[2*i + 1] = nmod_add(c2, c1, ctx);
        v[2*i + 2] = nmod_sub(c2, c1, ctx);
    }
}

static void _to_coeffs2(
    mp_limb_t * a,          /* length 2d+1 */
    const mp_limb_t * v,    /* length 2d+1 */
    mp_limb_t * t,          /* length 2d   */
    const mp_limb_t * Q,    /* length d by 2d+1 */
    slong d,
    nmod_t ctx)
{
    slong i, j;

    a[0] = v[0];

    for (i = 0; i < d; i++)
    {
        t[2*i+0] = nmod_sub(v[1+2*i+0], v[1+2*i+1], ctx);
        t[2*i+1] = nmod_add(v[1+2*i+0], v[1+2*i+1], ctx);
    }

    for (i = 0; i < d; i++)
    {
        mp_limb_t c1h, c1m, c1;
        mp_limb_t c2h, c2m, c2;
        c1h = c1m = c1 = 0;
        c2h = c2m = c2 = 0;
        umul_ppmm(c2m, c2, Q[0], v[0]);
        for (j = 0; j < d; j++)
        {
            MAC(c1h, c1m, c1, t[2*j + 0], Q[2*j + 1]);
            MAC(c2h, c2m, c2, t[2*j + 1], Q[2*j + 2]);
        }

        Q += 2*d + 1;

        NMOD_RED3(a[2*i + 1], c1h, c1m, c1, ctx);
        NMOD_RED3(a[2*i + 2], c2h, c2m, c2, ctx);
    }
}


void nmod_eval_interp_init(nmod_eval_interp_t E)
{
    E->M = NULL;
    E->T = NULL;
    E->Q = NULL;
    E->array = NULL;
    E->alloc = 0;
    E->d = 0;
    E->radix = 0;
    E->w = 0;
}

void nmod_eval_interp_clear(nmod_eval_interp_t E)
{
    if (E->alloc > 0)
        flint_free(E->array);
}

int nmod_eval_interp_set_degree_modulus(
    nmod_eval_interp_t E,
    slong deg,
    nmod_t ctx)
{
    slong d, new_alloc;
    mp_limb_t p = ctx.n;

    FLINT_ASSERT(deg >= 0);

    if (p < 3 || (p % 2) == 0 || deg >= p)
        return 0;

    if ((p % 4) == 1)
    {
        d = (deg + 3)/4;
        d = FLINT_MAX(d, 1);

        new_alloc = d*(4*d) + 4*d + d*(4*d + 1);

        if (E->alloc > 0)
            E->array = flint_realloc(E->array, new_alloc*sizeof(mp_limb_t));
        else
            E->array = flint_malloc(new_alloc*sizeof(mp_limb_t));

        E->radix = 4;
        E->alloc = new_alloc;
        E->d = d;
        E->M = E->array;
        E->T = E->M + d*(4*d);
        E->Q = E->T + 4*d;
        E->w = n_sqrtmod(p - 1, p);

        return _fill_matrices4(E->M, E->Q, d, ctx);
    }
    else
    {
        d = (deg + 1)/2;
        d = FLINT_MAX(d, 1);

        new_alloc = d*(2*d) + 2*d + d*(2*d + 1);

        if (E->alloc > 0)
            E->array = flint_realloc(E->array, new_alloc*sizeof(mp_limb_t));
        else
            E->array = flint_malloc(new_alloc*sizeof(mp_limb_t));

        E->radix = 2;
        E->alloc = new_alloc;
        E->d = d;
        E->M = E->array;
        E->T = E->M + d*(2*d);
        E->Q = E->T + 2*d;
        E->w = -UWORD(1);

        return _fill_matrices2(E->M, E->Q, d, ctx);
    }
}

void nmod_eval_interp_to_coeffs(
    mp_limb_t * a,
    const mp_limb_t * v,
    nmod_eval_interp_t E,
    nmod_t ctx)
{
    if (E->radix == 4)
        _to_coeffs4(a, v, E->T, E->Q, E->d, E->w, ctx);
    else
        _to_coeffs2(a, v, E->T, E->Q, E->d, ctx);
}

void nmod_eval_interp_from_coeffs(
    mp_limb_t * v,
    const mp_limb_t * a,
    slong alen,
    nmod_eval_interp_t E,
    nmod_t ctx)
{
    if (E->radix == 4)
        _from_coeffs4(v, a, alen, E->M, E->d, E->w, ctx);
    else
        _from_coeffs2(v, a, alen, E->M, E->d, ctx);
}
