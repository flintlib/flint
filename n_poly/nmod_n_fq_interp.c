/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"

/*
    conversion between polynomials in coefficient form and point-value form
    and artihmetic in point-value form
*/

#define MAC(h, m, l, a, b)                          \
{                                                   \
    mp_limb_t p1, p0;                               \
    umul_ppmm(p1, p0, a, b);                        \
    add_sssaaaaaa(h, m, l, h, m, l, 0, p1, p0);     \
}


/**************** conversion ************************************************/

/* p = 1 mod 4 */

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


static void _from_coeffs4_n_fq(
    mp_limb_t * v,          /* length 4d+1 */
    const mp_limb_t * a,
    slong alen,
    const mp_limb_t * M_,    /* length d by 4d */
    slong D,
    mp_limb_t w,
    slong d,
    nmod_t ctx)
{
    slong i, j, k;
    const mp_limb_t * Mrow;

    FLINT_ASSERT(0 <= alen);
    FLINT_ASSERT(alen <= 1 + 4*D);

    if (alen <= 1)
    {
        if (alen == 1)
        {
            for (i = 0; i < 4*D + 1; i++)
                _n_fq_set(v + d*i, a + d*0, d);
        }
        else
        {
            _nmod_vec_zero(v, d*(4*D+1));
        }

        return;
    }

    _n_fq_set(v + d*0, a + d*0, d);

    for (k = 0; k < d; k++)
    {
        Mrow = M_;
        for (i = 0; i < D; i++)
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
                FLINT_ASSERT(j + 3 < 4*D);
                MAC(c1h, c1m, c1, a[d*(j+1)+k], Mrow[j+0]);
                MAC(c2h, c2m, c2, a[d*(j+2)+k], Mrow[j+1]);
                MAC(c3h, c3m, c3, a[d*(j+3)+k], Mrow[j+2]);
                MAC(c4h, c4m, c4, a[d*(j+4)+k], Mrow[j+3]);
            }

            MAC(c1h, c1m, c1, j + 1 < alen ? a[d*(j+1)+k] : 0, Mrow[j+0]);
            MAC(c2h, c2m, c2, j + 2 < alen ? a[d*(j+2)+k] : 0, Mrow[j+1]);
            MAC(c3h, c3m, c3, j + 3 < alen ? a[d*(j+3)+k] : 0, Mrow[j+2]);
            FLINT_ASSERT(j + 4 >= alen);

            NMOD_RED3(c4, c4h, c4m, c4, ctx);
            NMOD_RED3(c1, c1h, c1m, c1, ctx);
            NMOD_RED3(c2, c2h, c2m, c2, ctx);
            NMOD_RED3(c3, c3h, c3m, c3, ctx);

            Mrow += 4*D;

            c4 = nmod_add(c4, a[d*0+k], ctx);
            t1 = nmod_add(c4, c2, ctx);
            t2 = nmod_sub(c4, c2, ctx);
            t3 = nmod_add(c1, c3, ctx);
            t4 = nmod_mul(nmod_sub(c1, c3, ctx), w, ctx);

            v[d*(4*i+1)+k] = nmod_add(t1, t3, ctx);
            v[d*(4*i+2)+k] = nmod_add(t2, t4, ctx);
            v[d*(4*i+3)+k] = nmod_sub(t1, t3, ctx);
            v[d*(4*i+4)+k] = nmod_sub(t2, t4, ctx);
        }
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

static void _to_coeffs4_n_fq(
    mp_limb_t * a,          /* length 4D+1 */
    const mp_limb_t * v,    /* length 4D+1 */
    mp_limb_t * t,          /* length 4D   */
    const mp_limb_t * Q_,   /* length D by 4D+1 */
    slong D,
    mp_limb_t w,
    slong d,
    nmod_t ctx)
{
    slong i, j, k;
    const mp_limb_t * Qrow;

    _n_fq_set(a + d*0, v + d*0, d);

    for (k = 0; k < d; k++)
    {

        for (i = 0; i < D; i++)
        {
            mp_limb_t t2 = nmod_add(v[d*(1+4*i+0)+k], v[d*(1+4*i+2)+k], ctx);
            mp_limb_t t1 = nmod_sub(v[d*(1+4*i+0)+k], v[d*(1+4*i+2)+k], ctx);
            mp_limb_t t3 = nmod_add(v[d*(1+4*i+1)+k], v[d*(1+4*i+3)+k], ctx);
            mp_limb_t t4 = nmod_mul(nmod_sub(v[d*(1+4*i+1)+k], v[d*(1+4*i+3)+k], ctx), w, ctx);
            t[4*i+0] = nmod_sub(t1, t4, ctx);
            t[4*i+1] = nmod_sub(t2, t3, ctx);
            t[4*i+2] = nmod_add(t1, t4, ctx);
            t[4*i+3] = nmod_add(t2, t3, ctx);
        }

        Qrow = Q_;
        for (i = 0; i < D; i++)
        {
            mp_limb_t c1h, c1m, c1;
            mp_limb_t c2h, c2m, c2;
            mp_limb_t c3h, c3m, c3;
            mp_limb_t c4h, c4m, c4;
            c1h = c1m = c1 = 0;
            c2h = c2m = c2 = 0;
            c3h = c3m = c3 = 0;
            c4h = c4m = c4 = 0;
            umul_ppmm(c4m, c4, Qrow[0], v[d*0+k]);
            for (j = 0; j < D; j++)
            {
                MAC(c1h, c1m, c1, t[4*j+0], Qrow[4*j+1]);
                MAC(c2h, c2m, c2, t[4*j+1], Qrow[4*j+2]);
                MAC(c3h, c3m, c3, t[4*j+2], Qrow[4*j+3]);
                MAC(c4h, c4m, c4, t[4*j+3], Qrow[4*j+4]);
            }

            Qrow += 4*D + 1;

            NMOD_RED3(a[d*(4*i+1)+k], c1h, c1m, c1, ctx);
            NMOD_RED3(a[d*(4*i+2)+k], c2h, c2m, c2, ctx);
            NMOD_RED3(a[d*(4*i+3)+k], c3h, c3m, c3, ctx);
            NMOD_RED3(a[d*(4*i+4)+k], c4h, c4m, c4, ctx);
        }
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

static void _from_coeffs2_n_fq(
    mp_limb_t * v,          /* length 4D+1 */
    const mp_limb_t * a,    /* length alen <= 2D+1 */
    slong alen,
    const mp_limb_t * M_,    /* length D by 4D */
    slong D,
    slong d,
    nmod_t ctx)
{
    slong i, j, k;
    const mp_limb_t * Mrow;

    FLINT_ASSERT(0 <= alen);
    FLINT_ASSERT(alen <= 1 + 2*D);

    if (alen <= 1)
    {
        if (alen == 1)
        {
            for (i = 0; i < 2*D+1; i++)
                _n_fq_set(v + d*i, a + d*0, d);
        }
        else
        {
            _nmod_vec_zero(v, d*(2*D+1));
        }

        return;
    }

    _n_fq_set(v + d*0, a + d*0, d);

    for (k = 0; k < d; k++)
    {
        Mrow = M_;
        for (i = 0; i < D; i++)
        {
            mp_limb_t c1h, c1m, c1;
            mp_limb_t c2h, c2m, c2;
            c1h = c1m = c1 = 0;
            c2h = c2m = c2 = 0;

            for (j = 0; j + 2 < alen; j += 2)
            {
                FLINT_ASSERT(j + 1 < 4*D);
                MAC(c1h, c1m, c1, a[d*(j+1)+k], Mrow[j+0]);
                MAC(c2h, c2m, c2, a[d*(j+2)+k], Mrow[j+1]);
            }

            MAC(c1h, c1m, c1, j + 1 < alen ? a[d*(j+1)+k] : 0, Mrow[j+0]);
            FLINT_ASSERT(j + 2 >= alen);

            NMOD_RED3(c2, c2h, c2m, c2, ctx);
            NMOD_RED3(c1, c1h, c1m, c1, ctx);

            Mrow += 2*D;

            c2 = nmod_add(c2, a[d*0+k], ctx);

            v[d*(2*i+1)+k] = nmod_add(c2, c1, ctx);
            v[d*(2*i+2)+k] = nmod_sub(c2, c1, ctx);
        }
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

static void _to_coeffs2_n_fq(
    mp_limb_t * a,          /* length 2d+1 */
    const mp_limb_t * v,    /* length 2d+1 */
    mp_limb_t * t,          /* length 2d   */
    const mp_limb_t * Q_,    /* length d by 2d+1 */
    slong D,
    slong d,
    nmod_t ctx)
{
    slong i, j, k;
    const mp_limb_t * Qrow;

    _n_fq_set(a + d*0, v + d*0, d);

    for (k = 0; k < d; k++)
    {
        for (i = 0; i < D; i++)
        {
            t[2*i+0] = nmod_sub(v[d*(1+2*i+0)+k], v[d*(1+2*i+1)+k], ctx);
            t[2*i+1] = nmod_add(v[d*(1+2*i+0)+k], v[d*(1+2*i+1)+k], ctx);
        }

        Qrow = Q_;
        for (i = 0; i < D; i++)
        {
            mp_limb_t c1h, c1m, c1;
            mp_limb_t c2h, c2m, c2;
            c1h = c1m = c1 = 0;
            c2h = c2m = c2 = 0;
            umul_ppmm(c2m, c2, Qrow[0], v[d*0+k]);
            for (j = 0; j < D; j++)
            {
                MAC(c1h, c1m, c1, t[2*j+0], Qrow[2*j+1]);
                MAC(c2h, c2m, c2, t[2*j+1], Qrow[2*j+2]);
            }

            Qrow += 2*D + 1;

            NMOD_RED3(a[d*(2*i+1)+k], c1h, c1m, c1, ctx);
            NMOD_RED3(a[d*(2*i+2)+k], c2h, c2m, c2, ctx);
        }
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

static void nmod_eval_interp_to_coeffs(
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

static void nmod_eval_interp_from_coeffs(
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

static void nmod_eval_interp_to_coeffs_n_fq(
    mp_limb_t * a,
    const mp_limb_t * v,
    nmod_eval_interp_t E,
    slong d,
    nmod_t ctx)
{
    if (E->radix == 4)
        _to_coeffs4_n_fq(a, v, E->T, E->Q, E->d, E->w, d, ctx);
    else
        _to_coeffs2_n_fq(a, v, E->T, E->Q, E->d, d, ctx);
}

static void nmod_eval_interp_from_coeffs_n_fq(
    mp_limb_t * v,
    const mp_limb_t * a,
    slong alen,
    nmod_eval_interp_t E,
    slong d,
    nmod_t ctx)
{
    if (E->radix == 4)
        _from_coeffs4_n_fq(v, a, alen, E->M, E->d, E->w, d, ctx);
    else
        _from_coeffs2_n_fq(v, a, alen, E->M, E->d, d, ctx);
}


/********** conversion over Fp **********/

void nmod_eval_interp_to_coeffs_poly(
    n_poly_t a,
    const n_poly_t v,
    nmod_eval_interp_t E,
    nmod_t ctx)
{
    slong l = nmod_eval_interp_eval_length(E);
    if (v->length == 0)
    {
        a->length = 0;
        return;
    }
    FLINT_ASSERT(v->length == l);
    n_poly_fit_length(a, l);
    nmod_eval_interp_to_coeffs(a->coeffs, v->coeffs, E, ctx);
    a->length = l;
    _n_poly_normalise(a);
}

void nmod_eval_interp_from_coeffs_poly(
    n_poly_t v,
    const n_poly_t a,
    nmod_eval_interp_t E,
    nmod_t ctx)
{
    slong l = nmod_eval_interp_eval_length(E);
    if (a->length == 0)
    {
        v->length = 0;
        return;
    }
    n_poly_fit_length(v, l);
    v->length = l;
    nmod_eval_interp_from_coeffs(v->coeffs, a->coeffs, a->length, E, ctx);
}

/********** conversion over Fq **********/

void nmod_eval_interp_to_coeffs_n_fq_poly(
    n_fq_poly_t a,
    const n_fq_poly_t v,
    nmod_eval_interp_t E,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong l = nmod_eval_interp_eval_length(E);

    if (v->length == 0)
    {
        a->length = 0;
        return;
    }

    FLINT_ASSERT(v->length == l);
    n_poly_fit_length(a, d*l);
    nmod_eval_interp_to_coeffs_n_fq(a->coeffs, v->coeffs, E, d, ctx->mod);
    a->length = l;
    _n_fq_poly_normalise(a, d);
}

void nmod_eval_interp_from_coeffs_n_fq_poly(
    n_fq_poly_t v,
    const n_fq_poly_t a,
    nmod_eval_interp_t E,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong l = nmod_eval_interp_eval_length(E);

    if (a->length == 0)
    {
        v->length = 0;
        return;
    }

    n_poly_fit_length(v, d*l);
    v->length = l;
    nmod_eval_interp_from_coeffs_n_fq(v->coeffs, a->coeffs, a->length, E, d, ctx->mod);
}


/*********** arithmetic *****************************************************/

/* a += b */
void nmod_evals_add_inplace(
    n_poly_t a,
    n_poly_t b,
    slong len,
    nmod_t ctx)
{
    slong i;

    if (b->length == 0)
        return;

    n_poly_fit_length(a, len);

    if (a->length == 0)
    {
        _nmod_vec_set(a->coeffs, b->coeffs, len);
        a->length = len;
        return;
    }

    for (i = 0; i < len; i++)
        a->coeffs[i] = nmod_add(a->coeffs[i], b->coeffs[i], ctx);

    a->length = _nmod_vec_is_zero(a->coeffs, len) ? 0 : len;
}


/* a = b*c */
void nmod_evals_mul(
    n_poly_t a,
    n_poly_t b,
    n_poly_t c,
    slong len,
    nmod_t ctx)
{
    slong i;

    if (b->length == 0 || c->length == 0)
    {
        a->length = 0;
        return;
    }

    n_poly_fit_length(a, len);

    for (i = 0; i < len; i++)
        a->coeffs[i] = nmod_mul(b->coeffs[i], c->coeffs[i], ctx);

    a->length = _nmod_vec_is_zero(a->coeffs, len) ? 0 : len;
}

/* a += b*c */
void nmod_evals_addmul(
    n_poly_t a,
    n_poly_t b,
    n_poly_t c,
    slong len,
    nmod_t ctx)
{
    slong i;

    if (b->length == 0 || c->length == 0)
        return;

    if (a->length == 0)
    {
        nmod_evals_mul(a, b, c, len, ctx);
        return;
    }

    for (i = 0; i < len; i++)
        NMOD_ADDMUL(a->coeffs[i], b->coeffs[i], c->coeffs[i], ctx);

    a->length = _nmod_vec_is_zero(a->coeffs, len) ? 0 : len;
}

void nmod_evals_fmma(
    n_poly_t a,
    n_poly_t b,
    n_poly_t c,
    n_poly_t d,
    n_poly_t e,
    slong len,
    nmod_t ctx)
{
    slong i;

    if (b->length == 0 || c->length == 0)
    {
        nmod_evals_mul(a, d, e, len, ctx);
        return;
    }

    if (d->length == 0 || e->length == 0)
    {
        nmod_evals_mul(a, b, c, len, ctx);
        return;
    }

    n_poly_fit_length(a, len);

    for (i = 0; i < len; i++)
    {
        mp_limb_t t = nmod_mul(b->coeffs[i], c->coeffs[i], ctx);
        NMOD_ADDMUL(t, d->coeffs[i], e->coeffs[i], ctx);
        a->coeffs[i] = t;
    }

    a->length = _nmod_vec_is_zero(a->coeffs, len) ? 0 : len;
}



/* a += b */
void n_fq_evals_add_inplace(
    n_fq_poly_t a,
    n_fq_poly_t b,
    slong len,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (b->length == 0)
        return;

    n_poly_fit_length(a, d*len);

    if (a->length == 0)
    {
        _nmod_vec_set(a->coeffs, b->coeffs, d*len);
        a->length = len;
        return;
    }

    _nmod_vec_add(a->coeffs, a->coeffs, b->coeffs, d*len, ctx->mod);

    a->length = _nmod_vec_is_zero(a->coeffs, d*len) ? 0 : len;
}

/* a = b*c */
void n_fq_evals_mul(
    n_fq_poly_t a,
    n_fq_poly_t b,
    n_fq_poly_t c,
    slong len,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * tmp;
    TMP_INIT;

    if (b->length == 0 || c->length == 0)
    {
        a->length = 0;
        return;
    }

    n_poly_fit_length(a, d*len);

    TMP_START;

    tmp = (mp_limb_t *) TMP_ALLOC(d*N_FQ_MUL_ITCH*sizeof(mp_limb_t));

    for (i = 0; i < len; i++)
        _n_fq_mul(a->coeffs + d*i, b->coeffs + d*i, c->coeffs + d*i, ctx, tmp);

    a->length = _nmod_vec_is_zero(a->coeffs, d*len) ? 0 : len;

    TMP_END;
}

/* a += b*c */
void n_fq_evals_addmul(
    n_fq_poly_t a,
    n_fq_poly_t b,
    n_fq_poly_t c,
    slong len,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * tmp;
    TMP_INIT;

    if (b->length == 0 || c->length == 0)
        return;

    if (a->length == 0)
    {
        n_fq_evals_mul(a, b, c, len, ctx);
        return;
    }

    TMP_START;

    tmp = (mp_limb_t *) TMP_ALLOC(d*N_FQ_MUL_ITCH*sizeof(mp_limb_t));

    for (i = 0; i < len; i++)
        _n_fq_addmul(a->coeffs + d*i, a->coeffs + d*i,
                            b->coeffs + d*i, c->coeffs + d*i, ctx, tmp);

    a->length = _nmod_vec_is_zero(a->coeffs, d*len) ? 0 : len;

    TMP_END;
}

void n_fq_evals_fmma(
    n_fq_poly_t a,
    n_fq_poly_t b,
    n_fq_poly_t c,
    n_fq_poly_t f,
    n_fq_poly_t e,
    slong len,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * tmp, * t;
    TMP_INIT;

    if (b->length == 0 || c->length == 0)
    {
        n_fq_evals_mul(a, f, e, len, ctx);
        return;
    }

    if (f->length == 0 || e->length == 0)
    {
        n_fq_evals_mul(a, b, c, len, ctx);
        return;
    }

    n_poly_fit_length(a, d*len);

    TMP_START;

    tmp = (mp_limb_t *) TMP_ALLOC(d*(1 + N_FQ_MUL_ITCH)*sizeof(mp_limb_t));
    t = tmp + d*N_FQ_MUL_ITCH;

    for (i = 0; i < len; i++)
    {
        _n_fq_mul(t, b->coeffs + d*i, c->coeffs + d*i, ctx, tmp);
        _n_fq_addmul(a->coeffs + d*i, t, f->coeffs + d*i,
                                          e->coeffs + d*i, ctx, tmp);
    }

    a->length = _nmod_vec_is_zero(a->coeffs, d*len) ? 0 : len;

    TMP_END;
}

