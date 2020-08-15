/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

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

    list = (mp_limb_t *) flint_malloc(d*sizeof(mp_limb_t));
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

    if (alen < 1)
    {
        flint_mpn_zero(v, 4*d+1);
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
    mp_limb_t * v,          /* length 4d+1 */
    const mp_limb_t * a,    /* length alen <= 2d+1 */
    slong alen,
    const mp_limb_t * M,    /* length d by 4d */
    slong d,
    nmod_t ctx)
{
    slong i, j;

    FLINT_ASSERT(0 <= alen);
    FLINT_ASSERT(alen <= 1 + 2*d);

    if (alen < 1)
    {
        flint_mpn_zero(v, 2*d+1);
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

typedef struct {
    mp_limb_t * M;
    mp_limb_t * T;
    mp_limb_t * Q;
    mp_limb_t * array;
    slong alloc;
    slong d;
    slong radix;
    mp_limb_t w;
} nmod_eval_interp_struct;

typedef nmod_eval_interp_struct nmod_eval_interp_t[1];


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

int nmod_eval_interp_eval_length(nmod_eval_interp_t E)
{
    return 1 + E->radix*E->d;
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


/*****************************************************************************/

/* M is a d x d array */
void computeM(
    mp_limb_t * M,
    slong d,
    const nmod_t ctx)
{
    slong i, j;
    mp_limb_t i2, i2j;

    for (i = 1; i <= d; i++)
    {
        i2 = nmod_mul(i, i, ctx);
        i2j = 1;
        for (j = 1; j <= d; j++)
        {
            i2j = nmod_mul(i2j, i2, ctx);
            M[j - 1 + d*(i - 1)] = i2j;
        }
    }
}

void coeffs2ptvals(
    mp_limb_t * Aptvals,
    slong d,
    const mp_limb_t * Acoeffs,
    slong Alength,
    const mp_limb_t * M,
    const nmod_t ctx)
{
    slong i, j;
    mp_limb_t vp0, vp1, vp2, vm0, vm1, vm2, pp1, pp0;

    if (Alength <= 1)
    {
        vp0 = Alength == 1 ? Acoeffs[0] : 0;
        for (i = 0; i <= 2*d; i++)
            Aptvals[i] = vp0;
        return;
    }

    Aptvals[0] = Acoeffs[0];

    for (i = 1; i <= d; i++)
    {
        vp2 = 0; vp1 = 0; vp0 = Acoeffs[0];
        vm2 = 0; vm1 = 0; vm0 = Acoeffs[1];
        for (j = 1; j < Alength/2; j++)
        {
            umul_ppmm(pp1, pp0, Acoeffs[2*j + 0], M[j - 1]);
            add_sssaaaaaa(vp2, vp1, vp0, vp2, vp1, vp0, 0, pp1, pp0);
            umul_ppmm(pp1, pp0, Acoeffs[2*j + 1], M[j - 1]);
            add_sssaaaaaa(vm2, vm1, vm0, vm2, vm1, vm0, 0, pp1, pp0);
        }
        NMOD_RED3(vm0, vm2, vm1, vm0, ctx);
        if ((Alength % 2) != 0)
        {
            umul_ppmm(pp1, pp0, Acoeffs[2*j + 0], M[j - 1]);
            add_sssaaaaaa(vp2, vp1, vp0, vp2, vp1, vp0, 0, pp1, pp0);
        }
        NMOD_RED3(vp0, vp2, vp1, vp0, ctx);
        vm0 = nmod_mul(i, vm0, ctx);
        Aptvals[2*i - 1] = nmod_add(vp0, vm0, ctx);
        Aptvals[2*i - 0] = nmod_sub(vp0, vm0, ctx);
        M += d;
    }
}

/* L is a (2*d + 1) x (d + 1) array */
void computeL(
    mp_limb_t * L,
    slong d,
    const nmod_t ctx)
{
    slong i, j;
    nmod_poly_t f, g;
    mp_limb_t c;

    nmod_poly_init_mod(f, ctx);
    nmod_poly_init_mod(g, ctx);

    nmod_poly_set_coeff_ui(f, 1, 1);
    nmod_poly_set_coeff_ui(g, 2, 1);
    for (i = 1; i <= d; i++)
    {
        nmod_poly_set_coeff_ui(g, 0, nmod_neg(nmod_mul(i, i, ctx), ctx));
        nmod_poly_mul(f, f, g);
    }
    for (i = 0; i <= d; i++)
    {
        nmod_poly_zero(g);
        nmod_poly_set_coeff_ui(g, 1, 1);
        nmod_poly_set_coeff_ui(g, 0, nmod_neg(i, ctx));
        nmod_poly_div(g, f, g);
        c = nmod_poly_evaluate_nmod(g, i);
        c = nmod_inv(c, ctx);
        nmod_poly_scalar_mul_nmod(g, g, c);
        for (j = 0; j <= 2*d; j++)
        {
            L[j*(d + 1) + i] = nmod_poly_get_coeff_ui(g, j);
        }
    }

    nmod_poly_clear(f);
    nmod_poly_clear(g);
}

void ptvals2coeffs(
    mp_limb_t * Acoeffs,
    slong Alength,
    const mp_limb_t * Aptvals,
    slong d,
    const mp_limb_t * L,
    mp_limb_t * t,
    const nmod_t ctx)
{
    slong i, j;
    mp_limb_t vp2, vp1, vp0, vm2, vm1, vm0, pp1, pp0;

    t[0] = Aptvals[0];
    for (i = 1; i <= d; i++)
    {
        vp0 = Aptvals[2*i - 1];
        vm0 = Aptvals[2*i - 0];
        t[2*i - 1] = nmod_add(vp0, vm0, ctx);
        t[2*i - 0] = nmod_sub(vp0, vm0, ctx);
    }

    for (i = 0; i + 2 <= Alength; i += 2)
    {
        vp2 = 0; umul_ppmm(vp1, vp0, L[0], t[0]);
        vm2 = 0; umul_ppmm(vm1, vm0, L[d + 1], t[0]);
        for (j = 1; j <= d; j++)
        {
            umul_ppmm(pp1, pp0, t[2*j - 1], L[j]);
            add_sssaaaaaa(vp2, vp1, vp0, vp2, vp1, vp0, 0, pp1, pp0);
            umul_ppmm(pp1, pp0, t[2*j - 0], L[j + d + 1]);
            add_sssaaaaaa(vm2, vm1, vm0, vm2, vm1, vm0, 0, pp1, pp0);
        }
        NMOD_RED3(Acoeffs[i + 0], vp2, vp1, vp0, ctx);
        NMOD_RED3(Acoeffs[i + 1], vm2, vm1, vm0, ctx);
        L += 2*(d + 1);
    }

    if (i < Alength)
    {
        vp2 = 0; umul_ppmm(vp1, vp0, L[0], t[0]);
        for (j = 1; j <= d; j++)
        {
            umul_ppmm(pp1, pp0, t[2*j - 1], L[j]);
            add_sssaaaaaa(vp2, vp1, vp0, vp2, vp1, vp0, 0, pp1, pp0);
        }
        NMOD_RED3(Acoeffs[i], vp2, vp1, vp0, ctx);
    }
}

void nmod_poly_to_ptvals(
    nmod_poly_t v,
    slong d,
    const nmod_poly_t a,
    mp_limb_t * M,
    nmod_t ctx)
{
    nmod_poly_fit_length(v, 2*d + 1);
    v->length = 2*d + 1;
    coeffs2ptvals(v->coeffs, d, a->coeffs, a->length, M, ctx);
}

void nmod_ptvals_to_poly(
    nmod_poly_t a,
    const mp_limb_t * v,
    slong d,
    mp_limb_t * L,
    nmod_t ctx)
{
    mp_limb_t * t = flint_malloc((2*d + 1)*sizeof(mp_limb_t));
    nmod_poly_fit_length(a, 2*d + 1);
    ptvals2coeffs(a->coeffs, 2*d + 1, v, d, L, t, ctx);
    a->length = 2*d + 1;
    _nmod_poly_normalise(a);
    flint_free(t);
}


void _nmod_vec_addmul(
    mp_limb_t * a,
    const mp_limb_t * b,
    const mp_limb_t * c,
    slong length,
    nmod_t ctx)
{
    slong i;
    for (i = 0; i < length; i++)
        a[i] = nmod_add(a[i], nmod_mul(b[i], c[i], ctx), ctx);
}
