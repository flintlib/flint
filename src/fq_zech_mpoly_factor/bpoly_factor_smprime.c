/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "fq_zech.h"
#include "fq_zech_poly.h"
#include "fmpz_poly_factor.h"
#include "fq_zech_poly_factor.h"
#include "nmod_mpoly_factor.h"
#include "fq_zech_mpoly_factor.h"

int fq_zech_next(fq_zech_t x, const fq_zech_ctx_t FLINT_UNUSED(ctx))
{
    if (x->value < 1)
        return 0;
    x->value--;
    return 1;
}

void fq_zech_bpoly_eval_var1(
    fq_zech_poly_t E,
    const fq_zech_bpoly_t A,
    const fq_zech_t alpha,
    const fq_zech_ctx_t ctx)
{
    slong i;
    fq_zech_poly_fit_length(E, A->length, ctx);
    for (i = 0; i < A->length; i++)
        fq_zech_poly_evaluate_fq_zech(E->coeffs + i, A->coeffs + i, alpha, ctx);
    E->length = A->length;
    _fq_zech_poly_normalise(E, ctx);
}

void fq_zech_bpoly_make_monic_series(
    fq_zech_bpoly_t A,
    const fq_zech_bpoly_t B,
    slong order,
    const fq_zech_ctx_t ctx)
{
    slong i;
    fq_zech_poly_t lcinv;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(fq_zech_bpoly_is_canonical(B, ctx));

    fq_zech_poly_init(lcinv, ctx);
    fq_zech_poly_inv_series(lcinv, B->coeffs + B->length - 1, order, ctx);

    fq_zech_bpoly_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
        fq_zech_poly_mullow(A->coeffs + i, B->coeffs + i, lcinv, order, ctx);

    A->length = B->length;
    fq_zech_bpoly_normalise(A, ctx);

    fq_zech_poly_clear(lcinv, ctx);
}


/****************** lifting **************************************************/

static void _hensel_build_tree(
    slong * link,
    fq_zech_bpoly_struct * v,
    fq_zech_bpoly_struct * w,
    const fq_zech_poly_struct * local_facs,
    slong r,
    const fq_zech_ctx_t ctx)
{
    slong i, j;
    fq_zech_poly_t d;
    fq_zech_poly_struct * V;
    fq_zech_poly_struct * W;

    V = FLINT_ARRAY_ALLOC(2*r - 2, fq_zech_poly_struct);
    W = FLINT_ARRAY_ALLOC(2*r - 2, fq_zech_poly_struct);

    fq_zech_poly_init(d, ctx);
    for (i = 0; i < 2*r - 2; i++)
    {
        fq_zech_poly_init(V + i, ctx);
        fq_zech_poly_init(W + i, ctx);
    }

    for (i = 0; i < r; i++)
    {
        fq_zech_poly_set(V + i, local_facs + i, ctx);
        link[i] = -i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        slong s, minp, mind;

        minp = j;
        mind = fq_zech_poly_degree(V + j, ctx);
        for (s = j + 1; s < i; s++)
        {
            if (fq_zech_poly_degree(V + s, ctx) < mind)
            {
                minp = s;
                mind = fq_zech_poly_degree(V + s, ctx);
            }
        }
        fq_zech_poly_swap(V + j, V + minp, ctx);
        FLINT_SWAP(slong, link[j], link[minp]);

        minp = j + 1;
        mind = fq_zech_poly_degree(V + j + 1, ctx);
        for (s = j + 2; s < i; s++)
        {
            if (fq_zech_poly_degree(V + s, ctx) < mind)
            {
                minp = s;
                mind = fq_zech_poly_degree(V + s, ctx);
            }
        }
        fq_zech_poly_swap(V + j + 1, V + minp, ctx);
        FLINT_SWAP(slong, link[j + 1], link[minp]);

        fq_zech_poly_mul(V + i, V + j, V + j + 1, ctx);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        fq_zech_poly_xgcd(d, W + j, W + j + 1, V + j, V + j + 1, ctx);
        FLINT_ASSERT(fq_zech_poly_is_one(d, ctx));
    }

    for (j = 0; j < 2*r - 2; j++)
    {
        fq_zech_bpoly_set_fq_zech_poly_var0(v + j, V + j, ctx);
        fq_zech_bpoly_set_fq_zech_poly_var0(w + j, W + j, ctx);
    }

    fq_zech_poly_clear(d, ctx);
    for (i = 0; i < 2*r - 2; i++)
    {
        fq_zech_poly_clear(V + i, ctx);
        fq_zech_poly_clear(W + i, ctx);
    }
    flint_free(V);
    flint_free(W);
}

static void _hensel_lift_fac(
    fq_zech_bpoly_t G,
    fq_zech_bpoly_t H,
    const fq_zech_bpoly_t f,
    fq_zech_bpoly_t g,
    fq_zech_bpoly_t h,
    const fq_zech_bpoly_t a,
    const fq_zech_bpoly_t b,
    slong p0,
    slong p1,
    const fq_zech_ctx_t ctx)
{
    slong i;
    fq_zech_bpoly_t c, t1, t2, q, r;

    fq_zech_bpoly_init(c, ctx);
    fq_zech_bpoly_init(t1, ctx);
    fq_zech_bpoly_init(t2, ctx);
    fq_zech_bpoly_init(q, ctx);
    fq_zech_bpoly_init(r, ctx);

    fq_zech_bpoly_mul(t1, g, h, ctx);
    fq_zech_bpoly_sub(c, f, t1, ctx);

    for (i = 0; i < c->length; i++)
    {
    #if FLINT_WANT_ASSERT
        {
            slong j;
            for (j = 0; j < p0; j++)
            {
                FLINT_ASSERT(j >= c->coeffs[i].length ||
                                fq_zech_is_zero(c->coeffs[i].coeffs + j, ctx));
            }
        }
    #endif
        fq_zech_poly_shift_right(c->coeffs + i, c->coeffs + i, p0, ctx);
        fq_zech_poly_truncate(c->coeffs + i, p1, ctx);
    }

    fq_zech_bpoly_mul_series(t1, c, b, p1, ctx);
    fq_zech_bpoly_divrem_series(q, r, t1, g, p1, ctx);
    for (i = 0; i < r->length; i++)
        fq_zech_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < g->length; i++)
        fq_zech_poly_truncate(g->coeffs + i, p0, ctx);
    fq_zech_bpoly_add(t1, r, g, ctx);

    fq_zech_bpoly_mul_series(t2, c, a, p1, ctx);
    fq_zech_bpoly_divrem_series(q, r, t2, h, p1, ctx);
    for (i = 0; i < r->length; i++)
        fq_zech_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < h->length; i++)
        fq_zech_poly_truncate(h->coeffs + i, p0, ctx);
    fq_zech_bpoly_add(t2, r, h, ctx);

    fq_zech_bpoly_swap(G, t1, ctx);
    fq_zech_bpoly_swap(H, t2, ctx);

#if FLINT_WANT_ASSERT
    {
        slong j;
        fq_zech_bpoly_mul(t1, G, H, ctx);
        fq_zech_bpoly_sub(c, f, t1, ctx);

        for (i = 0; i < c->length; i++)
        for (j = 0; j < p0 + p1; j++)
        {
            FLINT_ASSERT(j >= c->coeffs[i].length ||
                                fq_zech_is_zero(c->coeffs[i].coeffs + j, ctx));
        }
    }
#endif

    fq_zech_bpoly_clear(c, ctx);
    fq_zech_bpoly_clear(t1, ctx);
    fq_zech_bpoly_clear(t2, ctx);
    fq_zech_bpoly_clear(q, ctx);
    fq_zech_bpoly_clear(r, ctx);
}

static void _hensel_lift_inv(
    fq_zech_bpoly_t A,
    fq_zech_bpoly_t B,
    const fq_zech_bpoly_t G,
    const fq_zech_bpoly_t H,
    fq_zech_bpoly_t a,
    fq_zech_bpoly_t b,
    slong p0,
    slong p1,
    const fq_zech_ctx_t ctx)
{
    slong i;
    fq_zech_bpoly_t c, t1, t2, q, r;

    fq_zech_bpoly_init(c, ctx);
    fq_zech_bpoly_init(t1, ctx);
    fq_zech_bpoly_init(t2, ctx);
    fq_zech_bpoly_init(q, ctx);
    fq_zech_bpoly_init(r, ctx);

    for (i = 0; i < b->length; i++)
        fq_zech_poly_truncate(b->coeffs + i, p0, ctx);
    for (i = 0; i < a->length; i++)
        fq_zech_poly_truncate(a->coeffs + i, p0, ctx);

    fq_zech_bpoly_mul(t1, G, a, ctx);
    fq_zech_bpoly_mul(t2, H, b, ctx);
    fq_zech_bpoly_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        fq_zech_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
    fq_zech_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    fq_zech_bpoly_normalise(c, ctx);

    for (i = 0; i < c->length; i++)
    {
    #if FLINT_WANT_ASSERT
        {
            slong j;
            for (j = 0; j < p0; j++)
            {
                FLINT_ASSERT(j >= c->coeffs[i].length ||
                             fq_zech_is_zero(c->coeffs[i].coeffs + j, ctx));
            }
        }
    #endif
        fq_zech_poly_shift_right(c->coeffs + i, c->coeffs + i, p0, ctx);
        fq_zech_poly_truncate(c->coeffs + i, p1, ctx);
    }

    fq_zech_bpoly_mul_series(t1, c, b, p1, ctx);
    fq_zech_bpoly_divrem_series(q, r, t1, G, p1, ctx);

    for (i = 0; i < r->length; i++)
        fq_zech_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);

    fq_zech_bpoly_add(t1, r, b, ctx);

    fq_zech_bpoly_mul_series(t2, c, a, p1, ctx);
    fq_zech_bpoly_divrem_series(q, r, t2, H, p1, ctx);

    for (i = 0; i < r->length; i++)
        fq_zech_poly_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);

    fq_zech_bpoly_add(t2, r, a, ctx);

    fq_zech_bpoly_swap(t1, B, ctx);
    fq_zech_bpoly_swap(t2, A, ctx);

#if FLINT_WANT_ASSERT
    fq_zech_bpoly_mul(t1, G, A, ctx);
    fq_zech_bpoly_mul(t2, H, B, ctx);
    fq_zech_bpoly_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        fq_zech_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
    fq_zech_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    fq_zech_bpoly_normalise(c, ctx);

    {
        slong j;
        for (i = 0; i < c->length; i++)
        for (j = 0; j < p0 + p1; j++)
        {
            FLINT_ASSERT(j >= c->coeffs[i].length ||
                         fq_zech_is_zero(c->coeffs[i].coeffs + j, ctx));
        }
    }
#endif

    fq_zech_bpoly_clear(c, ctx);
    fq_zech_bpoly_clear(t1, ctx);
    fq_zech_bpoly_clear(t2, ctx);
    fq_zech_bpoly_clear(q, ctx);
    fq_zech_bpoly_clear(r, ctx);
}


static void _hensel_lift_tree(
    int opt,
    slong * link,
    fq_zech_bpoly_struct * v,
    fq_zech_bpoly_struct * w,
    const fq_zech_bpoly_t f,
    slong j,
    slong p0,
    slong p1,
    const fq_zech_ctx_t ctx)
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
    fq_zech_bpoly_struct ** lifted_fac;
    fq_zech_tpoly_t tmp;
    fq_zech_bpoly_t bmp;
    slong r;
    slong fac_lift_order;
    slong inv_lift_order;
    nmod_eval_interp_t E;
    int Eok;
    int use_linear;
} fq_zech_bpoly_lift_struct;

typedef fq_zech_bpoly_lift_struct fq_zech_bpoly_lift_t[1];

static void fq_zech_bpoly_lift_init(
    fq_zech_bpoly_lift_t L,
    const fq_zech_ctx_t ctx)
{
    L->link = NULL;
    L->lifted_fac = NULL;
    fq_zech_tpoly_init(L->tmp, ctx);
    fq_zech_bpoly_init(L->bmp, ctx);
    L->r = 0;
    L->fac_lift_order = 0;
    L->inv_lift_order = 0;
}

static void fq_zech_bpoly_lift_clear(
    fq_zech_bpoly_lift_t L,
    const fq_zech_ctx_t ctx)
{
    flint_free(L->link);
    flint_free(L->lifted_fac);
    fq_zech_tpoly_clear(L->tmp, ctx);
    fq_zech_bpoly_clear(L->bmp, ctx);
}

static void fq_zech_bpoly_lift_start(
    fq_zech_bpoly_lift_t L,
    fq_zech_poly_struct * local_facs,
    slong r,
    const fq_zech_ctx_t ctx)
{
    slong i;

    L->r = r;
    L->lifted_fac = (fq_zech_bpoly_struct **) flint_realloc(L->lifted_fac,
                                             r*sizeof(fq_zech_bpoly_struct *));

    L->fac_lift_order = 1;
    L->inv_lift_order = 1;

    {
        fq_zech_bpoly_struct * v, * w;

        L->link = (slong *) flint_realloc(L->link, (2*r - 2)*sizeof(slong));

        fq_zech_tpoly_fit_length(L->tmp, 2*(2*r - 2), ctx);
        v = L->tmp->coeffs + 0;
        w = L->tmp->coeffs + (2*r - 2);

        _hensel_build_tree(L->link, v, w, local_facs, r, ctx);
        for (i = 0; i < 2*r - 2; i++)
            if (-L->link[i] - 1 >= 0)
                L->lifted_fac[-L->link[i] - 1] = v + i;

        return;
    }
}

static void fq_zech_bpoly_lift_continue(
    fq_zech_bpoly_lift_t L,
    const fq_zech_bpoly_t monicA,
    slong order,
    const fq_zech_ctx_t ctx)
{
    slong i;
    slong r = L->r;

    if (order <= L->fac_lift_order)
        return;

    {
        fq_zech_bpoly_struct * v = L->tmp->coeffs + 0;
        fq_zech_bpoly_struct * w = L->tmp->coeffs + (2*r - 2);
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
}


/************* lattice reduction *********************************************/

/*
    The rows of N, if they are 0-1, correspond to combinations of the g that
    give factors of A. Compute {A/g[i]*g[i]'}_i and use the CLD bounds to
    try to make N smaller.
*/
static void _lattice(
    nmod_mat_t N,
    fq_zech_bpoly_struct * const * g,
    slong r,
    slong lift_order,
    slong * starts,
    const fq_zech_bpoly_t f,
    const fq_zech_ctx_t ctx)
{
    slong deg = fq_zech_ctx_degree(ctx);
    slong i, j, k, l;
    fq_zech_bpoly_t Q, R, dg;
    fq_zech_bpoly_struct * ld;
    nmod_mat_t M, T1, T2;
    ulong * trow;

    const dot_params_t params = _nmod_vec_dot_params(r, fq_zech_ctx_mod(ctx));
    trow = (ulong *) flint_malloc(r*sizeof(ulong));
    fq_zech_bpoly_init(Q, ctx);
    fq_zech_bpoly_init(R, ctx);
    fq_zech_bpoly_init(dg, ctx);
    ld = FLINT_ARRAY_ALLOC(r, fq_zech_bpoly_struct);
    for (i = 0; i < r; i++)
        fq_zech_bpoly_init(ld + i, ctx);

    /* init done */

    for (i = 0; i < r; i++)
    {
        fq_zech_bpoly_divrem_series(Q, R, f, g[i], lift_order, ctx);
        FLINT_ASSERT(R->length == 0);
        fq_zech_bpoly_derivative(R, g[i], ctx);
        fq_zech_bpoly_mul_series(ld + i, Q, R, lift_order, ctx);
    }

    for (k = 0; k + 1 < f->length; k++)
    {
        slong d = nmod_mat_nrows(N);

        if (lift_order <= starts[k])
            continue;

        nmod_mat_init(M, deg*(lift_order - starts[k]), d, fq_zech_ctx_mod(ctx).n);

        for (j = starts[k]; j < lift_order; j++)
        for (l = 0; l < deg; l++)
        {
            for (i = 0; i < r; i++)
            {
                if (k >= ld[i].length ||
                    j >= ld[i].coeffs[k].length)
                {
                    trow[i] = 0;
                }
                else
                {
                    nmod_poly_t junk;
                    nmod_poly_init_mod(junk, fq_zech_ctx_mod(ctx));
                    fq_zech_get_nmod_poly(junk, ld[i].coeffs[k].coeffs + j, ctx);
                    trow[i] = nmod_poly_get_coeff_ui(junk, l);
                    nmod_poly_clear(junk);
                }
            }

            for (i = 0; i < d; i++)
                nmod_mat_entry(M, (j - starts[k])*deg + l, i) =
              _nmod_vec_dot(trow, N->rows[i], r, fq_zech_ctx_mod(ctx), params);
        }

        nmod_mat_init_nullspace_tr(T1, M);

        nmod_mat_init(T2, nmod_mat_nrows(T1), nmod_mat_ncols(N),
                                                       fq_zech_ctx_mod(ctx).n);
        nmod_mat_mul(T2, T1, N);
        nmod_mat_swap(T2, N);
        nmod_mat_rref(N);

        nmod_mat_clear(M);
        nmod_mat_clear(T1);
        nmod_mat_clear(T2);
    }

    flint_free(trow);
    fq_zech_bpoly_clear(Q, ctx);
    fq_zech_bpoly_clear(R, ctx);
    fq_zech_bpoly_clear(dg, ctx);
    for (i = 0; i < r; i++)
        fq_zech_bpoly_clear(ld + i, ctx);
    flint_free(ld);
}

/**************** recombination **********************************************/

/*
    First multiply the g[i] into the gprod[i] according the the rows of N.
    Then, run zassenhaus on the gprod[i] as factors of A.
*/
static int _zassenhaus(
    const zassenhaus_prune_t FLINT_UNUSED(zas),
    slong limit,
    fq_zech_tpoly_t F,
    const fq_zech_t malpha,
    const nmod_mat_t N,
    fq_zech_bpoly_struct * const * loc_fac_org,
    slong r,
    slong order,
    const fq_zech_bpoly_t B,
    const fq_zech_ctx_t ctx)
{
    int success;
    fq_zech_bpoly_t Q, R, t1, t2;
    fq_zech_poly_t g;
    slong * subset;
    slong i, j, s, len, d = nmod_mat_nrows(N);
    fq_zech_bpoly_struct * loc_fac;
    fq_zech_bpoly_struct * f;
    fq_zech_bpoly_t B_copy;

    FLINT_ASSERT(nmod_mat_ncols(N) == r);

    fq_zech_poly_init(g, ctx);
    fq_zech_bpoly_init(Q, ctx);
    fq_zech_bpoly_init(R, ctx);
    fq_zech_bpoly_init(t1, ctx);
    fq_zech_bpoly_init(t2, ctx);
    fq_zech_bpoly_init(B_copy, ctx);

    subset = FLINT_ARRAY_ALLOC(d, slong);
    loc_fac = FLINT_ARRAY_ALLOC(d, fq_zech_bpoly_struct);
    for (i = 0; i < d; i++)
    {
        subset[i] = i;
        fq_zech_bpoly_init(loc_fac + i, ctx);
    }

    for (i = 0; i < d; i++)
    {
        fq_zech_bpoly_one(loc_fac + i, ctx);
        for (j = 0; j < r; j++)
        {
            if (nmod_mat_entry(N, i, j) == 0)
                continue;
            FLINT_ASSERT(nmod_mat_entry(N, i, j) == 1);
            fq_zech_bpoly_mul_series(t1, loc_fac + i, loc_fac_org[j], order, ctx);
            fq_zech_bpoly_swap(t1, loc_fac + i, ctx);
        }
    }

    f = (fq_zech_bpoly_struct *) B;

    len = d;
    for (s = 1; s <= len/2; s++)
    {
        if (s > limit)
        {
            success = 0;
            goto cleanup;
        }

        zassenhaus_subset_first(subset, len, s);
        while (1)
        {
            FLINT_ASSERT(f->length > 0);
            fq_zech_bpoly_set_fq_zech_poly_var1(t1, f->coeffs + f->length - 1, ctx);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    fq_zech_bpoly_mul_series(t2, t1, loc_fac + subset[i], order, ctx);
                    fq_zech_bpoly_swap(t1, t2, ctx);
                }
            }

            fq_zech_bpoly_make_primitive(g, t1, ctx);
            if (fq_zech_bpoly_divides(Q, f, t1, ctx))
            {
                fq_zech_bpoly_taylor_shift_var1(t1, t1, malpha, ctx);
                fq_zech_tpoly_fit_length(F, F->length + 1, ctx);
                fq_zech_bpoly_swap(F->coeffs + F->length, t1, ctx);
                F->length++;
                f = B_copy;
                fq_zech_bpoly_swap(f, Q, ctx);
                len -= s;
                if (!zassenhaus_subset_next_disjoint(subset, len + s))
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
        fq_zech_tpoly_fit_length(F, F->length + 1, ctx);
        fq_zech_bpoly_taylor_shift_var1(F->coeffs + F->length, f, malpha, ctx);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(f->length == 1);
        FLINT_ASSERT(fq_zech_poly_is_one(f->coeffs + 0, ctx));
    }

    success = 1;

cleanup:

    fq_zech_poly_clear(g, ctx);
    fq_zech_bpoly_clear(Q, ctx);
    fq_zech_bpoly_clear(R, ctx);
    fq_zech_bpoly_clear(t1, ctx);
    fq_zech_bpoly_clear(t2, ctx);
    fq_zech_bpoly_clear(B_copy, ctx);

    for (i = 0; i < d; i++)
        fq_zech_bpoly_clear(loc_fac + i, ctx);
    flint_free(loc_fac);
    flint_free(subset);

    return success;
}


int fq_zech_bpoly_factor_smprime(
    fq_zech_poly_t c,
    fq_zech_tpoly_t F,
    fq_zech_bpoly_t A,
    int allow_shift,
    const fq_zech_ctx_t ctx)
{
    int success;
    slong i, r;
    slong Alenx, Aleny;
    slong final_order, lift_order;
    slong * CLD;
    fq_zech_t alpha, Alc;
    fq_zech_poly_t Aeval;
    fq_zech_poly_factor_t local_fac;
    fq_zech_bpoly_t monicA;
    nmod_mat_t N;
    slong old_nrows;
    slong zas_limit;
    zassenhaus_prune_t zas;
    fq_zech_bpoly_lift_t L;

    fq_zech_bpoly_make_primitive(c, A, ctx);

    Alenx = 1 + fq_zech_bpoly_degree0(A, ctx);
    Aleny = 1 + fq_zech_bpoly_degree1(A, ctx);

    FLINT_ASSERT(Alenx > 1);

    fq_zech_init(alpha, ctx);
    fq_zech_init(Alc, ctx);
    fq_zech_poly_init(Aeval, ctx);
    fq_zech_poly_factor_init(local_fac, ctx);
    fq_zech_bpoly_init(monicA, ctx);
    nmod_mat_init(N, 0, 0, fq_zech_ctx_mod(ctx).n);
    CLD = FLINT_ARRAY_ALLOC(Alenx, slong);
    zassenhaus_prune_init(zas);
    fq_zech_bpoly_lift_init(L, ctx);

    /* TODO: CLD bounds */
    for (i = 0; i < Alenx; i++)
        CLD[i] = Aleny;

    zassenhaus_prune_set_degree(zas, Alenx - 1);

    fq_zech_zero(alpha, ctx);
    goto got_alpha;

next_alpha:

    if (!allow_shift || !fq_zech_next(alpha, ctx))
    {
        success = 0;
        goto cleanup;
    }

got_alpha:

    fq_zech_bpoly_eval_var1(Aeval, A, alpha, ctx);

    /* if killed leading coeff, get new alpha */
    if (Aeval->length != Alenx)
        goto next_alpha;

    /* note the constant term of Aeval can be zero */

    fq_zech_poly_factor(local_fac, Alc, Aeval, ctx);
    r = local_fac->num;

    zassenhaus_prune_start_add_factors(zas);
    for (i = 0; i < r; i++)
        zassenhaus_prune_add_factor(zas, local_fac->poly[i].length - 1,
                                         local_fac->exp[i]);
    zassenhaus_prune_end_add_factors(zas);

    /* check for irreducibility */
    if ((r < 2 && local_fac->exp[0] == 1) ||
         zassenhaus_prune_must_be_irreducible(zas))
    {
        fq_zech_tpoly_fit_length(F, 1, ctx);
        fq_zech_bpoly_swap(F->coeffs + 0, A, ctx);
        F->length = 1;
        success = 1;
        goto cleanup;
    }

    /* if multiple factors, get new alpha */
    for (i = 0; i < r; i++)
    {
        if (local_fac->exp[i] != 1)
            goto next_alpha;
    }

    /* done if B is constant in y */
    if (Aleny < 2)
    {
        fq_zech_tpoly_fit_length(F, r, ctx);
        F->length = r;
        for (i = 0; i < r; i++)
            fq_zech_bpoly_set_fq_zech_poly_var0(F->coeffs + i, local_fac->poly + i, ctx);
        success = 1;
        goto cleanup;
    }

    /* precision for constructing true factors */
    final_order = Aleny;

    fq_zech_bpoly_taylor_shift_var1(A, A, alpha, ctx);

    fq_zech_bpoly_lift_start(L, local_fac->poly, r, ctx);

    /* precision for lifted local factors */
    lift_order = final_order + r;
    fq_zech_bpoly_make_monic_series(monicA, A, lift_order, ctx);
    fq_zech_bpoly_lift_continue(L, monicA, lift_order, ctx);

    /* the rows of N give the combinations of local factors */
    nmod_mat_clear(N);
    nmod_mat_init(N, r, r, fq_zech_ctx_mod(ctx).n);
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
    fq_zech_neg(alpha, alpha, ctx);
    success = _zassenhaus(zas, zas_limit, F, alpha, N,
                                        L->lifted_fac, r, final_order, A, ctx);
    fq_zech_neg(alpha, alpha, ctx);
    if (success)
        goto cleanup;

    /* first attempt failed, try subsets of size 1 or 2 from now on */
    zas_limit = 2;

more:

    /* increase precision until N has fewer rows and is a nice 0-1 mat */
    old_nrows = nmod_mat_nrows(N);
    _lattice(N, L->lifted_fac, r, lift_order, CLD, A, ctx);
    if (nmod_mat_nrows(N) < old_nrows && nmod_mat_is_reduced(N))
        goto try_zas;

    lift_order += r;
    fq_zech_bpoly_make_monic_series(monicA, A, lift_order, ctx);
    fq_zech_bpoly_lift_continue(L, monicA, lift_order, ctx);
    goto more;

cleanup:

    fq_zech_bpoly_lift_clear(L, ctx);

    flint_free(CLD);

    nmod_mat_clear(N);
    fq_zech_clear(alpha, ctx);
    fq_zech_clear(Alc, ctx);
    fq_zech_poly_clear(Aeval, ctx);
    fq_zech_poly_factor_clear(local_fac, ctx);
    fq_zech_bpoly_clear(monicA, ctx);

    zassenhaus_prune_clear(zas);

    return success;
}
