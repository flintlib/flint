/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"

#define FLINT_TMP_ARRAY_ALLOC(n, T) (T *) TMP_ALLOC(n*sizeof(T))


static void n_bpoly_reverse_gens(n_bpoly_t a, const n_bpoly_t b)
{
    slong i, j;
    n_bpoly_zero(a);
    for (i = 0; i < b->length; i++)
    {
        const n_fq_poly_struct * bi = b->coeffs + i;
        for (j = 0; j < bi->length; j++)
        {
            n_bpoly_set_coeff(a, j, i, bi->coeffs[j]);
        }
    }
}

static void n_bpoly_mod_make_monic_series(
    n_bpoly_t A,
    const n_bpoly_t B,
    slong order,
    nmod_t ctx)
{
    slong i;
    n_poly_t lcinv;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B, ctx));

    n_poly_init(lcinv);
    n_poly_mod_inv_series(lcinv, B->coeffs + B->length - 1, order, ctx);

    n_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
        n_poly_mod_mullow(A->coeffs + i, B->coeffs + i, lcinv, order, ctx);

    A->length = B->length;
    n_bpoly_normalise(A);

    n_poly_clear(lcinv);
}


static void n_bpoly_mod_eval(
    nmod_poly_t E,
    const n_bpoly_t A,
    mp_limb_t alpha,
    nmod_t ctx)
{
    slong i;
    n_poly_t alphapow;

    nmod_poly_zero(E);

    if (alpha == 0)
    {
        for (i = A->length - 1; i >= 0; i--)
            nmod_poly_set_coeff_ui(E, i, n_poly_get_coeff(A->coeffs + i, 0));
        return;
    }

    n_poly_init2(alphapow, 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    for (i = A->length - 1; i >= 0; i--)
        nmod_poly_set_coeff_ui(E, i, n_poly_mod_eval_pow(A->coeffs + i, alphapow, ctx));

    n_poly_clear(alphapow);
}


/****************** lifting **************************************************/


static void _hensel_build_tree(
    slong * link,
    n_bpoly_struct * v,
    n_bpoly_struct * w,
    const nmod_poly_struct * local_facs,
    slong r,
    nmod_t ctx)
{
    slong i, j;
    n_poly_t d;
    n_poly_struct * V;
    n_poly_struct * W;

    V = FLINT_ARRAY_ALLOC((2*r - 2), n_poly_struct);
    W = FLINT_ARRAY_ALLOC((2*r - 2), n_poly_struct);

    n_poly_init(d);
    for (i = 0; i < 2*r - 2; i++)
    {
        n_poly_init(V + i);
        n_poly_init(W + i);
    }

    for (i = 0; i < r; i++)
    {
        n_poly_set_nmod_poly(V + i, local_facs + i);
        link[i] = -i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        slong s, minp, mind;

        minp = j;
        mind = n_poly_degree(V + j);
        for (s = j + 1; s < i; s++)
        {
            if (n_poly_degree(V + s) < mind)
            {
                minp = s;
                mind = n_poly_degree(V + s);
            }
        }
        n_poly_swap(V + j, V + minp);
        SLONG_SWAP(link[j], link[minp]);

        minp = j + 1;
        mind = n_poly_degree(V + j + 1);
        for (s = j + 2; s < i; s++)
        {
            if (n_poly_degree(V + s) < mind)
            {
                minp = s;
                mind = n_poly_degree(V + s);
            }
        }
        n_poly_swap(V + j + 1, V + minp);
        SLONG_SWAP(link[j + 1], link[minp]);

        n_poly_mod_mul(V + i, V + j, V + j + 1, ctx);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        n_poly_mod_xgcd(d, W + j, W + j + 1, V + j, V + j + 1, ctx);
        FLINT_ASSERT(n_poly_is_one(d));
    }

    for (j = 0; j < 2*r - 2; j++)
    {
        n_bpoly_set_poly_gen0(v + j, V + j);
        n_bpoly_set_poly_gen0(w + j, W + j);
    }

    n_poly_clear(d);
    for (i = 0; i < 2*r - 2; i++)
    {
        n_poly_clear(V + i);
        n_poly_clear(W + i);
    }
    flint_free(V);
    flint_free(W);
}

static void _hensel_lift_fac(
    n_bpoly_t G,
    n_bpoly_t H,
    const n_bpoly_t f,
    n_bpoly_t g,
    n_bpoly_t h,
    const n_bpoly_t a,
    const n_bpoly_t b,
    slong p0,
    slong p1,
    nmod_t ctx)
{
    slong i, j;
    n_bpoly_t c, t1, t2, q, r;

    n_bpoly_init(c);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(q);
    n_bpoly_init(r);

    n_bpoly_mod_mul(t1, g, h, ctx);
    n_bpoly_mod_sub(c, f, t1, ctx);
    for (i = 0; i < c->length; i++)
    {
        for (j = 0; j < p0; j++)
            FLINT_ASSERT(n_poly_get_coeff(c->coeffs + i, j) == 0);
        n_poly_shift_right(c->coeffs + i, c->coeffs + i, p0);
        n_poly_truncate(c->coeffs + i, p1);
    }

    n_bpoly_mod_mul_series(t1, c, b, p1, ctx);
    n_bpoly_mod_divrem_series(q, r, t1, g, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_shift_left(r->coeffs + i, r->coeffs + i, p0);
    for (i = 0; i < g->length; i++)
        n_poly_truncate(g->coeffs + i, p0);
    n_bpoly_mod_add(t1, r, g, ctx);

    n_bpoly_mod_mul_series(t2, c, a, p1, ctx);
    n_bpoly_mod_divrem_series(q, r, t2, h, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_shift_left(r->coeffs + i, r->coeffs + i, p0);
    for (i = 0; i < h->length; i++)
        n_poly_truncate(h->coeffs + i, p0);
    n_bpoly_mod_add(t2, r, h, ctx);

    n_bpoly_swap(G, t1);
    n_bpoly_swap(H, t2);

#if FLINT_WANT_ASSERT
    n_bpoly_mod_mul(t1, G, H, ctx);
    n_bpoly_mod_sub(c, f, t1, ctx);
    for (i = 0; i < c->length; i++)
        for (j = 0; j < p0 + p1; j++)
            FLINT_ASSERT(n_poly_get_coeff(c->coeffs + i, j) == 0);
#endif

    n_bpoly_clear(c);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_bpoly_clear(q);
    n_bpoly_clear(r);
}

static void _hensel_lift_inv(
    n_bpoly_t A,
    n_bpoly_t B,
    const n_bpoly_t G,
    const n_bpoly_t H,
    n_bpoly_t a,
    n_bpoly_t b,
    slong p0,
    slong p1,
    nmod_t ctx)
{
    slong i, j;
    n_bpoly_t c, t1, t2, q, r;

    n_bpoly_init(c);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(q);
    n_bpoly_init(r);

    for (i = 0; i < a->length; i++)
        n_poly_truncate(a->coeffs + i, p0);
    for (i = 0; i < b->length; i++)
        n_poly_truncate(b->coeffs + i, p0);

    n_bpoly_mod_mul(t1, G, a, ctx);
    n_bpoly_mod_mul(t2, H, b, ctx);
    n_bpoly_mod_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        n_poly_mod_neg(c->coeffs + i, c->coeffs + i, ctx);
    n_poly_mod_add_ui(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    n_bpoly_normalise(c);

    for (i = 0; i < c->length; i++)
    {
        for (j = 0; j < p0; j++)
            FLINT_ASSERT(n_poly_get_coeff(c->coeffs + i, j) == 0);
        n_poly_shift_right(c->coeffs + i, c->coeffs + i, p0);
        n_poly_truncate(c->coeffs + i, p1);
    }

    n_bpoly_mod_mul_series(t1, c, b, p1, ctx);
    n_bpoly_mod_divrem_series(q, r, t1, G, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_shift_left(r->coeffs + i, r->coeffs + i, p0);
    n_bpoly_mod_add(t1, r, b, ctx);

    n_bpoly_mod_mul_series(t2, c, a, p1, ctx);
    n_bpoly_mod_divrem_series(q, r, t2, H, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_shift_left(r->coeffs + i, r->coeffs + i, p0);
    n_bpoly_mod_add(t2, r, a, ctx);

    n_bpoly_swap(t1, B);
    n_bpoly_swap(t2, A);

#if FLINT_WANT_ASSERT
    n_bpoly_mod_mul(t1, G, A, ctx);
    n_bpoly_mod_mul(t2, H, B, ctx);
    n_bpoly_mod_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        n_poly_mod_neg(c->coeffs + i, c->coeffs + i, ctx);
    n_poly_mod_add_ui(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    n_bpoly_normalise(c);

    for (i = 0; i < c->length; i++)
        for (j = 0; j < p0 + p1; j++)
            FLINT_ASSERT(n_poly_get_coeff(c->coeffs + i, j) == 0);
#endif

    n_bpoly_clear(c);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_bpoly_clear(q);
    n_bpoly_clear(r); 
}

static void _hensel_lift_tree(
    int opt,
    slong * link,
    n_bpoly_struct * v,
    n_bpoly_struct * w,
    const n_bpoly_t f,
    slong j,
    slong p0,
    slong p1,
    nmod_t ctx)
{
    FLINT_ASSERT(p1 <= p0);

    if (j < 0)
        return;

    if (opt >= 0)
        _hensel_lift_fac(v + j, v + j + 1,
                           f, v + j, v + j + 1, w + j, w + j + 1, p0, p1, ctx);

    if (opt <= 0)
        _hensel_lift_inv(w + j, w + j + 1,
                              v + j, v + j + 1, w + j, w + j + 1, p0, p1, ctx);

    _hensel_lift_tree(opt, link, v, w, v + j, link[j], p0, p1, ctx);
    _hensel_lift_tree(opt, link, v, w, v + j + 1, link[j + 1], p0, p1, ctx);
}

typedef struct {
    slong * link;
    n_bpoly_struct ** lifted_fac;
    n_tpoly_t tmp;
    n_bpoly_t bmp;
    slong r;
    slong fac_lift_order;
    slong inv_lift_order;
    nmod_eval_interp_t E;
    int Eok;
    int use_linear;
} n_bpoly_mod_lift_struct;

typedef n_bpoly_mod_lift_struct n_bpoly_mod_lift_t[1];

static void n_bpoly_mod_lift_init(n_bpoly_mod_lift_t L)
{
    L->link = NULL;
    L->lifted_fac = NULL;
    n_tpoly_init(L->tmp);
    n_bpoly_init(L->bmp);
    L->r = 0;
    L->fac_lift_order = 0;
    L->inv_lift_order = 0;
    nmod_eval_interp_init(L->E);
}

static void n_bpoly_mod_lift_clear(n_bpoly_mod_lift_t L)
{
    flint_free(L->link);
    flint_free(L->lifted_fac);
    n_tpoly_clear(L->tmp);
    n_bpoly_clear(L->bmp);
    nmod_eval_interp_clear(L->E);
}

static void n_bpoly_mod_lift_start(
    n_bpoly_mod_lift_t L,
    nmod_poly_struct * local_facs,
    slong r,    
    nmod_t ctx)
{
    slong i, k, degx;
    n_bpoly_struct * A, * Bfinal, * U, * Ue, * Be, * B;
    n_poly_struct * s, * Binv;
    n_poly_struct * c, * t;

    degx = 0;
    for (i = 0; i < r; i++)
        degx += local_facs[i].length - 1;

    L->r = r;
    L->lifted_fac = (n_bpoly_struct **) flint_realloc(L->lifted_fac,
                                                   r*sizeof(n_bpoly_struct *));

    /* linear lifting has large memory requirements wrt r */
    if (r < 30 + 5*FLINT_BIT_COUNT(degx))
        L->use_linear = 1;
    else
        L->use_linear = 0;

    L->fac_lift_order = 1;
    L->inv_lift_order = 1;

    if (!L->use_linear)
    {
        n_bpoly_struct * v, * w;

        L->link = (slong *) flint_realloc(L->link, (2*r - 2)*sizeof(slong));

        n_tpoly_fit_length(L->tmp, 2*(2*r - 2));
        v = L->tmp->coeffs + 0;
        w = L->tmp->coeffs + (2*r - 2);

        _hensel_build_tree(L->link, v, w, local_facs, r, ctx);
        for (i = 0; i < 2*r - 2; i++)
            if (-L->link[i] - 1 >= 0)
                L->lifted_fac[-L->link[i] - 1] = v + i;

        return;
    }

    n_tpoly_fit_length(L->tmp, 4*r + 1);
    A = L->tmp->coeffs;
    Bfinal = A + 1;
    U = Ue = Bfinal + r;
    B = U + r;
    Be = B + r;

    n_bpoly_fit_length(A, 1);
    A->length = 1;
    n_poly_one(A->coeffs + 0);
    for (k = 0; k < r; k++)
    {
        n_bpoly_fit_length(B + k, 1);
        B[k].length = 1;
        n_poly_set_nmod_poly(B[k].coeffs + 0, local_facs + k);
        n_poly_mod_mul(A->coeffs + 0, A->coeffs + 0, B[k].coeffs + 0, ctx);

        L->lifted_fac[k] = Bfinal + k;
        n_bpoly_reverse_gens(L->lifted_fac[k], B + k);
    }

    FLINT_ASSERT(degx == n_poly_degree(A->coeffs + 0));

    /* try evaluation when not too many local factors */
    if (r < 20 + 2*FLINT_BIT_COUNT(degx))
        L->Eok = nmod_eval_interp_set_degree_modulus(L->E, degx, ctx);
    else
        L->Eok = 0;

    n_bpoly_fit_length(L->bmp, 2*r + 5);
    s = L->bmp->coeffs;
    Binv = s + r;
    c = Binv + r;
    t = c + 1;

    for (k = 0; k < r; k++)
    {
        /* s[k] = (prod_{i!=k} B[i].coeffs[0])^-1 (mod B[k].coeffs[0]) */
        n_poly_mod_div(t, A->coeffs + 0, B[k].coeffs + 0, ctx);
        if (!n_poly_mod_invmod(s + k, t, B[k].coeffs + 0, ctx))
            flint_throw(FLINT_IMPINV, "n_bpoly_mod_lift_start: bad inverse");

        /* set up mul (mod B[k].coeffs[0]) */
        n_poly_reverse(t, B[k].coeffs + 0, B[k].coeffs[0].length);
        n_poly_mod_inv_series(Binv + k, t, B[k].coeffs[0].length, ctx);

        if (L->Eok)
        {
            n_bpoly_fit_length(Be + k, 1);
            nmod_eval_interp_from_coeffs_poly(Be[k].coeffs + 0,
                                              B[k].coeffs + 0, L->E, ctx);
            /* Ue[0] is not used */
            if (k > 0)
            {
                n_bpoly_fit_length(Ue + k, 1);
                Ue[k].length = 1;
                nmod_evals_zero(Ue[k].coeffs + 0);
            }
        }
        else
        {
            /* U[0] is not used */
            if (k > 0)
            {
                n_bpoly_fit_length(U + k, 1);
                U[k].length = 1;
                n_poly_zero(U[k].coeffs + 0);
            }
        }
    }

    if (r > 2)
    {
        if (L->Eok)
        {
            slong len = nmod_eval_interp_eval_length(L->E);
            k = r - 2;
            nmod_evals_mul(Ue[k].coeffs + 0, Be[k].coeffs + 0,
                                             Be[k + 1].coeffs + 0, len, ctx);
            for (k--; k > 0; k--)
                nmod_evals_mul(Ue[k].coeffs + 0, Be[k].coeffs + 0,
                                               Ue[k + 1].coeffs + 0, len, ctx);
        }
        else
        {
            k = r - 2;
            n_poly_mod_mul(U[k].coeffs + 0, B[k].coeffs + 0,
                                            B[k + 1].coeffs + 0, ctx);
            for (k--; k > 0; k--)
                n_poly_mod_mul(U[k].coeffs + 0, B[k].coeffs + 0,
                                                U[k + 1].coeffs + 0, ctx);
        }
    }
}

static void n_bpoly_mod_lift_continue(
    n_bpoly_mod_lift_t L,
    const n_bpoly_t monicA,
    slong order,
    nmod_t ctx)
{
    slong i, j, k;
    slong r = L->r;
    n_bpoly_struct * A, * Bfinal, * U, * Ue, * Be, * B;
    n_poly_struct * s, * Binv;
    n_poly_struct * c, * t, * ce, * vk, * vek;

    if (order <= L->fac_lift_order)
        return;

    if (!L->use_linear)
    {
        n_bpoly_struct * v = L->tmp->coeffs + 0;
        n_bpoly_struct * w = L->tmp->coeffs + (2*r - 2);
        slong e[FLINT_BITS+1];

        FLINT_ASSERT(1 <= L->inv_lift_order);
        FLINT_ASSERT(L->fac_lift_order <= 2*L->inv_lift_order);

        for (i = 0, e[i] = order; e[i] > L->fac_lift_order; i++)
            e[i+1] = (e[i] + 1)/2;
        e[i]   = L->fac_lift_order;
        e[i+1] = L->inv_lift_order;

        if (e[i+1] < e[i])
            _hensel_lift_tree(-1, L->link, v, w, monicA, 2*r-4, e[i+1], e[i]-e[i+1], ctx);

        for (i--; i > 0; i--)
            _hensel_lift_tree(0, L->link, v, w, monicA, 2*r-4, e[i+1], e[i]-e[i+1], ctx);

        _hensel_lift_tree(1, L->link, v, w, monicA, 2*r-4, e[1], e[0]-e[1], ctx);

        L->fac_lift_order = e[0];
        L->inv_lift_order = e[1];
        return;
    }

    A = L->tmp->coeffs;
    Bfinal = A + 1;
    U = Ue = Bfinal + r;
    B = U + r;
    Be = B + r;

    s = L->bmp->coeffs;
    Binv = s + r;
    c = Binv + r;
    t = c + 1;
    ce = t + 1;
    vk = ce + 1;
    vek = vk + 1;

    /* tack on reversal of monicA */
    for (i = 0; i < monicA->length; i++)
    {
        n_poly_struct * Bi = monicA->coeffs + i;
        j = FLINT_MIN(Bi->length, order);
        for (j--; j >= L->fac_lift_order; j--)
            n_bpoly_set_coeff(A, j, i, Bi->coeffs[j]);
    }

    /* tack on zeros to the B[k] */
    for (k = 0; k < r; k++)
    {
        n_bpoly_fit_length(B + k, order);
        if (L->Eok)
            n_bpoly_fit_length(Be + k, order);

        for (i = B[k].length; i < order; i++)
        {
            n_poly_zero(B[k].coeffs + i);
            if (L->Eok)
                nmod_evals_zero(Be[k].coeffs + i);
        }

        /* U[0] is not used, zero out both U and Ue */
        if (k > 0)
        {
            n_bpoly_fit_length(U + k, order);
            for (i = U[k].length; i < order; i++)
                U[k].coeffs[i].length = 0;
            U[k].length = order;
        }
    }

    if (L->Eok && r > 2)
    {
        slong len = nmod_eval_interp_eval_length(L->E);

        for (j = L->fac_lift_order; j < order; j++)
        {
            k = r - 2;
            nmod_evals_zero(Ue[k].coeffs + j);
            for (i = 0; i <= j; i++)
                nmod_evals_addmul(Ue[k].coeffs + j, Be[k].coeffs + i,
                                              Be[k + 1].coeffs + j - i, len, ctx);
            for (k--; k > 0; k--)
            {
                nmod_evals_zero(Ue[k].coeffs + j);
                for (i = 0; i <= j; i++)
                    nmod_evals_addmul(Ue[k].coeffs + j, Be[k].coeffs + i,
                                              Ue[k + 1].coeffs + j - i, len, ctx);
            }

            nmod_evals_zero(ce);
            for (i = 0; i <= j; i++)
                nmod_evals_addmul(ce, Be[0].coeffs + i,
                                      Ue[1].coeffs + j - i, len, ctx);

            nmod_eval_interp_to_coeffs_poly(c, ce, L->E, ctx);

            if (j < A->length)
                n_poly_mod_sub(c, A->coeffs + j, c, ctx);
            else
                n_poly_mod_neg(c, c, ctx);

            if (n_poly_is_zero(c))
                continue;

            for (k = r - 1; k >= 0; k--)
            {
                n_poly_mod_rem(t, c, B[k].coeffs + 0, ctx);
                n_poly_mod_mulmod_preinv(vk, s + k, t, B[k].coeffs + 0, Binv + k, ctx);
                nmod_eval_interp_from_coeffs_poly(vek, vk, L->E, ctx);
                if (!n_poly_is_zero(vk))
                {
                    nmod_evals_add_inplace(Be[k].coeffs + j, vek, len, ctx);
                    n_poly_mod_add(B[k].coeffs + j, B[k].coeffs + j, vk, ctx);
                    if (!n_poly_is_zero(B[k].coeffs + j))
                        B[k].length = FLINT_MAX(B[k].length, j + 1);
                }

                /* correct the U's */
                if (k > r - 2)
                {
                    n_poly_swap(ce, vek);
                }
                else if (k > 0)
                {
                    n_poly_struct * p;
                    p = (k == r - 2) ? Be[k + 1].coeffs : Ue[k + 1].coeffs;
                    nmod_evals_fmma(ce, Be[k].coeffs + 0, ce, p, vek, len, ctx);
                    nmod_evals_add_inplace(Ue[k].coeffs + j, ce, len, ctx);
                }
            }
        }
    }
    else if (L->Eok)
    {
        slong len = nmod_eval_interp_eval_length(L->E);

        FLINT_ASSERT(r == 2);

        for (j = L->fac_lift_order; j < order; j++)
        {
            nmod_evals_zero(ce);
            for (i = 0; i <= j; i++)
                nmod_evals_addmul(ce, Be[0].coeffs + i,
                                      Be[1].coeffs + j - i, len, ctx);

            nmod_eval_interp_to_coeffs_poly(c, ce, L->E, ctx);

            if (j < A->length)
                n_poly_mod_sub(c, A->coeffs + j, c, ctx);
            else
                n_poly_mod_neg(c, c, ctx);

            if (n_poly_is_zero(c))
                continue;

            for (k = 0; k < r; k++)
            {
                n_poly_mod_rem(t, c, B[k].coeffs + 0, ctx);
                n_poly_mod_mulmod_preinv(vk, s + k, t, B[k].coeffs + 0, Binv + k, ctx);
                nmod_eval_interp_from_coeffs_poly(vek, vk, L->E, ctx);
                if (!n_poly_is_zero(vk))
                {
                    nmod_evals_add_inplace(Be[k].coeffs + j, vek, len, ctx);
                    n_poly_mod_add(B[k].coeffs + j, B[k].coeffs + j, vk, ctx);
                    if (!n_poly_is_zero(B[k].coeffs + j))
                        B[k].length = FLINT_MAX(B[k].length, j + 1);
                }
            }
        }
    }
    else if (r > 2)
    {
        for (j = L->fac_lift_order; j < order; j++)
        {
            k = r - 2;
            n_poly_zero(U[k].coeffs + j);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length && j - i < B[k + 1].length)
                {
                    n_poly_mod_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
                    n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
                }
            }
            for (k--; k > 0; k--)
            {
                n_poly_zero(U[k].coeffs + j);
                for (i = 0; i <= j; i++)
                {
                    if (i < B[k].length)
                    {
                        n_poly_mod_mul(t, B[k].coeffs + i, U[k + 1].coeffs + j - i, ctx);
                        n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
                    }
                }
            }

            if (j < A->length)
                n_poly_set(c, A->coeffs + j);
            else
                n_poly_zero(c);

            for (i = 0; i <= j; i++)
            {
                if (i < B[0].length)
                {
                    n_poly_mod_mul(t, B[0].coeffs + i, U[1].coeffs + j - i, ctx);
                    n_poly_mod_sub(c, c, t, ctx);
                }
            }

            if (n_poly_is_zero(c))
                continue;

            for (k = r - 1; k >= 0; k--)
            {
                n_poly_mod_rem(t, c, B[k].coeffs + 0, ctx);
                n_poly_mod_mulmod_preinv(vk, s + k, t, B[k].coeffs + 0, Binv + k, ctx);
                if (!n_poly_is_zero(vk))
                {
                    n_poly_mod_add(B[k].coeffs + j, B[k].coeffs + j, vk, ctx);
                    if (!n_poly_is_zero(B[k].coeffs + j))
                        B[k].length = FLINT_MAX(B[k].length, j + 1);
                }

                /* correct the U's */
                if (k > r - 2)
                {
                    n_poly_swap(ce, vk);
                }
                else if (k > 0)
                {
                    n_poly_struct * p;
                    n_poly_mod_mul(t, B[k].coeffs + 0, ce, ctx);
                    p = (k == r - 2) ? B[k + 1].coeffs : U[k + 1].coeffs;
                    n_poly_mod_mul(ce, p, vk, ctx);
                    n_poly_mod_add(ce, ce, t, ctx);
                    n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, ce, ctx);
                }
            }
        }
    }
    else
    {
        FLINT_ASSERT(r == 2);

        for (j = L->fac_lift_order; j < order; j++)
        {
            if (j < A->length)
                n_poly_set(c, A->coeffs + j);
            else
                n_poly_zero(c);

            for (i = FLINT_MIN(j, B[0].length - 1); i >= 0; i--)
            {
                n_poly_mod_mul(t, B[0].coeffs + i, B[1].coeffs + j - i, ctx);
                n_poly_mod_sub(c, c, t, ctx);
            }

            if (n_poly_is_zero(c))
                continue;

            for (k = 0; k < r; k++)
            {
                n_poly_mod_rem(t, c, B[k].coeffs + 0, ctx);
                n_poly_mod_mulmod_preinv(vk, s + k, t, B[k].coeffs + 0, Binv + k, ctx);
                if (!n_poly_is_zero(vk))
                {
                    n_poly_mod_add(B[k].coeffs + j, B[k].coeffs + j, vk, ctx);
                    if (!n_poly_is_zero(B[k].coeffs + j))
                        B[k].length = FLINT_MAX(B[k].length, j + 1);
                }
            }
        }
    }

    L->fac_lift_order = order;

    for (k = 0; k < r; k++)
        n_bpoly_reverse_gens(Bfinal + k, B + k);

#if FLINT_WANT_ASSERT
    {
        n_bpoly_t t1, t2;
        n_bpoly_init(t1);
        n_bpoly_init(t2);
        n_bpoly_set(t1, Bfinal + 0);
        for (k = 1; k < r; k++)
        {
            n_bpoly_mod_mul_series(t2, t1, Bfinal + k, order, ctx);
            n_bpoly_swap(t1, t2);
        }
        n_bpoly_mod_sub(t2, monicA, t1, ctx);

        for (i = 0; i < t2->length; i++)
        {
            for (j = 0; j < FLINT_MIN(order, t2->coeffs[i].length); j++)
            {
                FLINT_ASSERT(0 == t2->coeffs[i].coeffs[j]);
            }
        }
        n_bpoly_clear(t1);
        n_bpoly_clear(t2);
    }
#endif
}


/************* lattice reduction ********************************************/

/*
    The rows of N, if they are 0-1, correspond to combinations of the g that
    give factors of A. Compute {A/g[i]*g[i]'}_i and use the CLD bounds to
    try to make N smaller.
*/
static void _lattice(
    nmod_mat_t N,
    n_bpoly_struct * const * g,
    slong r,
    slong lift_order,
    slong * CLD,
    const n_bpoly_t A,
    nmod_t ctx)
{
    slong i, j, k;
    n_bpoly_t Q, R, dg;
    n_bpoly_struct * ld;
    nmod_mat_t M, T1, T2;
    int nlimbs;
    mp_limb_t * trow;

    nlimbs = _nmod_vec_dot_bound_limbs(r, ctx);
    trow = FLINT_ARRAY_ALLOC(r, mp_limb_t);
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(dg);
    ld = FLINT_ARRAY_ALLOC(r, n_bpoly_struct);
    for (i = 0; i < r; i++)
    {
        n_bpoly_init(ld + i);
        n_bpoly_mod_divrem_series(Q, R, A, g[i], lift_order, ctx);
        FLINT_ASSERT(R->length == 0);
        n_bpoly_mod_derivative_gen0(R, g[i], ctx);
        n_bpoly_mod_mul_series(ld + i, Q, R, lift_order, ctx);
    }

    for (k = 0; k + 1 < A->length; k++)
    {
        slong nrows = nmod_mat_nrows(N);

        FLINT_ASSERT(nrows > 0);

        if (lift_order <= CLD[k])
            continue;

        nmod_mat_init(M, lift_order - CLD[k], nrows, ctx.n);

        for (j = CLD[k]; j < lift_order; j++)
        {
            for (i = 0; i < r; i++)
                trow[i] = n_bpoly_get_coeff(ld + i, k, j);

            for (i = 0; i < nrows; i++)
                nmod_mat_entry(M, j - CLD[k], i) =
                             _nmod_vec_dot(trow, N->rows[i], r, ctx, nlimbs);
        }

        nmod_mat_init_nullspace_tr(T1, M);

        nmod_mat_init(T2, nmod_mat_nrows(T1), nmod_mat_ncols(N), ctx.n);
        nmod_mat_mul(T2, T1, N);
        nmod_mat_swap(T2, N);
        nmod_mat_rref(N);

        nmod_mat_clear(M);
        nmod_mat_clear(T1);
        nmod_mat_clear(T2);
    }

    flint_free(trow);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(dg);
    for (i = 0; i < r; i++)
        n_bpoly_clear(ld + i);
    flint_free(ld);
}


/**************** recombination **********************************************/

/*
    First multiply the g[i] into the gprod[i] according the the rows of N.
    Then, run zassenhaus on the gprod[i] as factors of A.
*/
static int _zassenhaus(
    const zassenhaus_prune_t zas,
    slong limit,
    n_tpoly_t F,
    mp_limb_t malpha,
    const nmod_mat_t N,
    n_bpoly_struct * const * g,
    slong r,
    slong order,
    const n_bpoly_t A,
    nmod_t ctx)
{
    slong total_deg;
    int success;
    n_bpoly_t Q, R, t1, t2;
    n_poly_t cont;
    slong i, j, k, len, nrows = nmod_mat_nrows(N);
    slong * subset;
    n_bpoly_struct * gprod;
    n_bpoly_struct * f;
    n_bpoly_t A_copy;

    FLINT_ASSERT(nmod_mat_ncols(N) == r);

    n_poly_init(cont);
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(A_copy);
    gprod = FLINT_ARRAY_ALLOC(nrows, n_bpoly_struct);
    subset = FLINT_ARRAY_ALLOC(nrows, slong);
    for (i = 0; i < nrows; i++)
    {
        subset[i] = i;
        n_bpoly_init(gprod + i);
        n_bpoly_one(gprod + i);
        for (j = 0; j < r; j++)
        {
            if (nmod_mat_entry(N, i, j) == 0)
                continue;
            FLINT_ASSERT(nmod_mat_entry(N, i, j) == 1);
            n_bpoly_mod_mul_series(t1, gprod + i, g[j], order, ctx);
            n_bpoly_swap(t1, gprod + i);
        }
    }

    f = (n_bpoly_struct *) A;

    len = nrows;
    for (k = 1; k <= len/2; k++)
    {
        if (k > limit)
        {
            success = 0;
            goto cleanup;
        }

        zassenhaus_subset_first(subset, len, k);
        while (1)
        {
            total_deg = 0;
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                    total_deg += gprod[subset[i]].length - 1;
            }

            if (!zassenhaus_prune_degree_is_possible(zas, total_deg))
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
                continue;
            }

            FLINT_ASSERT(f->length > 0);
            n_bpoly_set_poly_gen1(t1, f->coeffs + f->length - 1);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    n_bpoly_mod_mul_series(t2, t1, gprod + subset[i], order, ctx);
                    n_bpoly_swap(t1, t2);
                }
            }

            n_bpoly_mod_make_primitive(cont, t1, ctx);
            if (n_bpoly_mod_divides(Q, f, t1, ctx))
            {
                n_bpoly_mod_taylor_shift_gen1(t1, t1, malpha, ctx);
                n_tpoly_fit_length(F, F->length + 1);
                n_bpoly_swap(F->coeffs + F->length, t1);
                F->length++;
                f = A_copy;
                n_bpoly_swap(f, Q);

                len -= k;
                if (!zassenhaus_subset_next_disjoint(subset, len + k))
                    break;
            }
            else
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
            }
        }
    }

    if (f->length > 1)
    {
        n_tpoly_fit_length(F, F->length + 1);
        n_bpoly_mod_taylor_shift_gen1(F->coeffs + F->length, f, malpha, ctx);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(f->length == 1);
        FLINT_ASSERT(n_poly_is_one(f->coeffs + 0));
    }

    success = 1;

cleanup:

    for (i = 0; i < nrows; i++)
        n_bpoly_clear(gprod + i);
    flint_free(gprod);

    flint_free(subset);

    n_poly_clear(cont);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_bpoly_clear(A_copy);

    return success;
}


/*****************************************************************************/

/*
    x = gen(0), y = gen(1).
    A is supposed to be separable wrt x.
    Put the content of A wrt x in c, and the factors in F.
    Return 1 for success, i.e. a good evaluation point y = alpha was found.
    If allow_shift is false, only y = 0 is tried.

    TODO:
        (1) compute better CLD bounds
        (2) start with a lower initial precision
        (3) restart the lifting when nrows(N) decreases
*/
int n_bpoly_mod_factor_smprime(
    n_poly_t c,     /* poly in gen(1) */
    n_tpoly_t F,
    n_bpoly_t A,    /* clobbered */
    int allow_shift,
    nmod_t ctx)
{
    int success;
    slong i, r;
    slong Alenx, Aleny;
    slong final_order, lift_order;
    slong * CLD;
    nmod_poly_t Aeval;
    mp_limb_t alpha_best, alpha_tmp;
    nmod_poly_factor_t local_fac_best, local_fac_tmp;
    int local_fac_tries = 0;
    n_bpoly_t monicA;
    nmod_mat_t N;
    slong old_nrows;
    slong zas_limit;    
    zassenhaus_prune_t zas;
    n_bpoly_mod_lift_t L;

    n_bpoly_mod_make_primitive(c, A, ctx);

    Alenx = 1 + n_bpoly_degree0(A);
    Aleny = 1 + n_bpoly_degree1(A);

    FLINT_ASSERT(Alenx > 1);

    nmod_poly_init_mod(Aeval, ctx);
    nmod_poly_factor_init(local_fac_best);
    nmod_poly_factor_init(local_fac_tmp);
    n_bpoly_init(monicA);
    nmod_mat_init(N, 0, 0, ctx.n);
    CLD = FLINT_ARRAY_ALLOC(Alenx, slong);
    zassenhaus_prune_init(zas);
    n_bpoly_mod_lift_init(L);

    /* TODO: CLD bounds */
    for (i = 0; i < Alenx; i++)
        CLD[i] = Aleny;

    zassenhaus_prune_set_degree(zas, Alenx - 1);

    alpha_tmp = 0;
    alpha_best = 0;
    goto got_alpha;

next_alpha:

    if (!allow_shift || alpha_tmp + 1 >= ctx.n)
    {
        if (local_fac_best->num > 0)
            goto doit;
        success = 0;
        goto cleanup;
    }

    alpha_tmp++;

got_alpha:

    n_bpoly_mod_eval(Aeval, A, alpha_tmp, ctx);

    /* if killed leading coeff, get new alpha */
    if (Aeval->length != Alenx)
        goto next_alpha;

    /* note the constant term of Aeval can be zero */

    local_fac_tmp->num = 0;
    nmod_poly_factor(local_fac_tmp, Aeval);
    r = local_fac_tmp->num;

    zassenhaus_prune_start_add_factors(zas);
    for (i = 0; i < r; i++)
        zassenhaus_prune_add_factor(zas, local_fac_tmp->p[i].length - 1,
                                         local_fac_tmp->exp[i]);
    zassenhaus_prune_end_add_factors(zas);

    if ((r < 2 && local_fac_tmp->exp[0] == 1) ||
         zassenhaus_prune_must_be_irreducible(zas))
    {
        n_tpoly_fit_length(F, 1);
        F->length = 1;
        n_bpoly_swap(F->coeffs + 0, A);
        success = 1;
        goto cleanup;
    }

    /* if multiple factors, get new alpha */
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(local_fac_tmp->p[i].length > 1);
        FLINT_ASSERT(local_fac_tmp->p[i].coeffs[local_fac_tmp->p[i].length - 1] == 1);
        if (local_fac_tmp->exp[i] != 1)
            goto next_alpha;
    }

    /* done if A is constant in y */
    if (Aleny < 2)
    {
        n_poly_t mock;
        n_tpoly_fit_length(F, r);
        F->length = r;
        for (i = 0; i < r; i++)
        {
            n_poly_mock(mock, local_fac_tmp->p + i);
            n_bpoly_set_poly_gen0(F->coeffs + i, mock);
        }
        success = 1;
        goto cleanup;
    }

    /* alpha_tmp & local_fac_tmp are good; update best */
    if (local_fac_best->num < 1 || local_fac_best->num > local_fac_tmp->num)
    {
        alpha_best = alpha_tmp;
        nmod_poly_factor_swap(local_fac_best, local_fac_tmp);
    }

    if (++local_fac_tries < 2)
        goto next_alpha;

doit:

    /* local_fac_best is a factorization mod (y - alpha_best) */
    r = local_fac_best->num;

    /* precision for constructing true factors */
    final_order = Aleny;

    n_bpoly_mod_taylor_shift_gen1(A, A, alpha_best, ctx);

    n_bpoly_mod_lift_start(L, local_fac_best->p, r, ctx);

    /* precision for lifted local factors */
    lift_order = final_order + r;
    n_bpoly_mod_make_monic_series(monicA, A, lift_order, ctx);
    n_bpoly_mod_lift_continue(L, monicA, lift_order, ctx);

    /* the rows of N give the combinations of local factors */    
    nmod_mat_clear(N);
    nmod_mat_init(N, r, r, ctx.n);
    for (i = 0; i < r; i++)
        nmod_mat_entry(N, i, i) = 1;

    /* size limit on subsets in zassenhaus combination */
    zas_limit = 1;

    _lattice(N, L->lifted_fac, r, lift_order, CLD, A, ctx);
    if (!nmod_mat_is_reduced(N))
        goto more;

try_zas:

    /* zassenhaus only make sense if N is a nice 0-1 mat */
    FLINT_ASSERT(nmod_mat_is_reduced(N));

    /* combine local factors according the rows of N, then by subsets */
    F->length = 0;
    success = _zassenhaus(zas, zas_limit, F, nmod_neg(alpha_best, ctx), N,
                                        L->lifted_fac, r, final_order, A, ctx);
    if (success)
        goto cleanup;

    /* first attemp failed, try subsets of size 1 or 2 from now on */
    zas_limit = 2;

more:

    /* increase precision until N has fewer rows and is a nice 0-1 mat */ 
    old_nrows = nmod_mat_nrows(N);
    _lattice(N, L->lifted_fac, r, lift_order, CLD, A, ctx);
    if (nmod_mat_nrows(N) < old_nrows && nmod_mat_is_reduced(N))
        goto try_zas;

    lift_order += r;
    n_bpoly_mod_make_monic_series(monicA, A, lift_order, ctx);
    n_bpoly_mod_lift_continue(L, monicA, lift_order, ctx);
    goto more;

cleanup:

    n_bpoly_mod_lift_clear(L);

    flint_free(CLD);

    nmod_mat_clear(N);
    nmod_poly_clear(Aeval);
    nmod_poly_factor_clear(local_fac_best);
    nmod_poly_factor_clear(local_fac_tmp);
    n_bpoly_clear(monicA);

    zassenhaus_prune_clear(zas);

    return success;
}
