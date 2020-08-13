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


void n_bpoly_fq_eval_var1(
    fq_nmod_poly_t E,
    const n_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    fq_nmod_t t;
    fq_nmod_poly_t s;

    fq_nmod_init(t, ctx);
    fq_nmod_poly_init(s, ctx);

    fq_nmod_poly_zero(E, ctx);
    for (i = A->length - 1; i >= 0; i--)
    {
        n_poly_fq_get_fq_nmod_poly(s, A->coeffs + i, ctx);
        fq_nmod_poly_evaluate_fq_nmod(t, s, alpha, ctx);
        fq_nmod_poly_set_coeff(E, i, t, ctx);
    }

    fq_nmod_clear(t, ctx);
    fq_nmod_poly_clear(s, ctx);
}

void n_bpoly_fq_make_monic_series(
    n_bpoly_t A,
    const n_bpoly_t B,
    slong order,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    n_poly_t lcinv;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(n_bpoly_fq_is_canonical(B, ctx));

    n_poly_init(lcinv);
    n_poly_fq_inv_series(lcinv, B->coeffs + B->length - 1, order, ctx);

    n_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
        n_poly_fq_mullow(A->coeffs + i, B->coeffs + i, lcinv, order, ctx);

    A->length = B->length;
    n_bpoly_normalise(A);

    n_poly_clear(lcinv);
}


static void _lattice(
    nmod_mat_t N,
    n_bpoly_struct * const * g,
    slong r,
    slong lift_order,
    slong * starts,
    const n_bpoly_t f,
    const fq_nmod_ctx_t ctx)
{
    slong deg = fq_nmod_ctx_degree(ctx);
    slong i, j, k, l;
    n_bpoly_t Q, R, dg;
    n_bpoly_struct * ld;
    nmod_mat_t M, T1, T2;
    int nlimbs;
    mp_limb_t * trow;

    nlimbs = _nmod_vec_dot_bound_limbs(r, ctx->mod);
    trow = (mp_limb_t *) flint_malloc(r*sizeof(mp_limb_t));
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(dg);
    ld = FLINT_ARRAY_ALLOC(r, n_bpoly_struct);
    for (i = 0; i < r; i++)
        n_bpoly_init(ld + i);

    /* init done */

    for (i = 0; i < r; i++)
    {
        n_bpoly_fq_divrem_series(Q, R, f, g[i], lift_order, ctx);
        FLINT_ASSERT(R->length == 0);
        n_bpoly_fq_derivative(R, g[i], ctx);
        n_bpoly_fq_mul_series(ld + i, Q, R, lift_order, ctx);
    }

    for (k = 0; k + 1 < f->length; k++)
    {
        slong d = nmod_mat_nrows(N);

        if (lift_order <= starts[k])
            continue;

        nmod_mat_init(M, deg*(lift_order - starts[k]), d, ctx->modulus->mod.n);

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
                    trow[i] = ld[i].coeffs[k].coeffs[deg*j + l];
                }
            }

            for (i = 0; i < d; i++)
                nmod_mat_entry(M, (j - starts[k])*deg + l, i) =
                 _nmod_vec_dot(trow, N->rows[i], r, ctx->modulus->mod, nlimbs);
        }

        nmod_mat_init_nullspace_tr(T1, M);

        nmod_mat_init(T2, nmod_mat_nrows(T1), nmod_mat_ncols(N),
                                                          ctx->modulus->mod.n);
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


static int _zassenhaus(
    slong limit,
    n_tpoly_t F,
    const fq_nmod_t malpha,
    const nmod_mat_t N,
    n_bpoly_struct * const * loc_fac_org,
    slong r,
    slong order,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    int success;
    n_bpoly_t Q, R, t1, t2;
    n_poly_t g;
    slong * subset;
    slong i, j, s, len, d = nmod_mat_nrows(N);
    n_bpoly_struct * loc_fac;
    n_bpoly_struct * f;
    n_bpoly_t B_copy;

    FLINT_ASSERT(nmod_mat_ncols(N) == r);

    n_poly_init(g);
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(B_copy);

    subset = FLINT_ARRAY_ALLOC(d, slong);
    loc_fac = FLINT_ARRAY_ALLOC(d, n_bpoly_struct);
    for (i = 0; i < d; i++)
    {
        subset[i] = i;
        n_bpoly_init(loc_fac + i);
    }

    for (i = 0; i < d; i++)
    {
        n_bpoly_fq_one(loc_fac + i, ctx);
        for (j = 0; j < r; j++)
        {
            if (nmod_mat_entry(N, i, j) == 0)
                continue;
            FLINT_ASSERT(nmod_mat_entry(N, i, j) == 1);
            n_bpoly_fq_mul_series(t1, loc_fac + i, loc_fac_org[j], order, ctx);
            n_bpoly_swap(t1, loc_fac + i);
        }
    }

    f = (n_bpoly_struct *) B;

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
            n_bpoly_fq_set_n_poly_fq_var1(t1, f->coeffs + f->length - 1, ctx);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    n_bpoly_fq_mul_series(t2, t1, loc_fac + subset[i], order, ctx);
                    n_bpoly_swap(t1, t2);
                }
            }

            n_bpoly_fq_make_primitive(g, t1, ctx);
            if (n_bpoly_fq_divides(Q, f, t1, ctx))
            {
                n_bpoly_fq_taylor_shift_var1_fq_nmod(t1, t1, malpha, ctx);
                n_tpoly_fit_length(F, F->length + 1);
                n_bpoly_swap(F->coeffs + F->length, t1);
                F->length++;
                f = B_copy;
                n_bpoly_swap(f, Q);
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
        n_tpoly_fit_length(F, F->length + 1);
        n_bpoly_fq_taylor_shift_var1_fq_nmod(F->coeffs + F->length, f, malpha, ctx);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(f->length == 1);
        FLINT_ASSERT(n_poly_fq_is_one(f->coeffs + 0, ctx));
    }

    success = 1;

cleanup:

    n_poly_clear(g);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_bpoly_clear(B_copy);

    for (i = 0; i < d; i++)
        n_bpoly_clear(loc_fac + i);
    flint_free(loc_fac);
    flint_free(subset);

    return success;
}

static void _hensel_build_tree(
    slong * link,
    n_bpoly_struct * v,
    n_bpoly_struct * w,
    const fq_nmod_poly_struct * local_facs,
    slong r,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    n_poly_t d;
    n_poly_struct * V;
    n_poly_struct * W;

    V = FLINT_ARRAY_ALLOC(2*r - 2, n_poly_struct);
    W = FLINT_ARRAY_ALLOC(2*r - 2, n_poly_struct);

    n_poly_init(d);
    for (i = 0; i < 2*r - 2; i++)
    {
        n_poly_init(V + i);
        n_poly_init(W + i);
    }

    for (i = 0; i < r; i++)
    {
        n_poly_fq_set_fq_nmod_poly(V + i, local_facs + i, ctx);
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

        n_poly_fq_mul(V + i, V + j, V + j + 1, ctx);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        n_poly_fq_xgcd(d, W + j, W + j + 1, V + j, V + j + 1, ctx);
        FLINT_ASSERT(n_poly_fq_is_one(d, ctx));
    }

    for (j = 0; j < 2*r - 2; j++)
    {
        n_bpoly_fq_set_n_poly_fq_var0(v + j, V + j, ctx);
        n_bpoly_fq_set_n_poly_fq_var0(w + j, W + j, ctx);
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
    const fq_nmod_ctx_t ctx)
{
    slong i;
    n_bpoly_t c, t1, t2, q, r;

    n_bpoly_init(c);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(q);
    n_bpoly_init(r);

    n_bpoly_fq_mul(t1, g, h, ctx);
    n_bpoly_fq_sub(c, f, t1, ctx);

    for (i = 0; i < c->length; i++)
    {
    #if WANT_ASSERT
        {
            slong j, d = fq_nmod_ctx_degree(ctx);
            for (j = 0; j < p0; j++)
            {
                FLINT_ASSERT(j >= c->coeffs[i].length ||
                             _n_fq_is_zero(c->coeffs[i].coeffs + d*j, d));
            }
        }
    #endif
        n_poly_fq_shift_right(c->coeffs + i, c->coeffs + i, p0, ctx);
        n_poly_fq_truncate(c->coeffs + i, p1, ctx);
    }

    n_bpoly_fq_mul_series(t1, c, b, p1, ctx);
    n_bpoly_fq_divrem_series(q, r, t1, g, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_fq_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < g->length; i++)
        n_poly_fq_truncate(g->coeffs + i, p0, ctx);
    n_bpoly_fq_add(t1, r, g, ctx);

    n_bpoly_fq_mul_series(t2, c, a, p1, ctx);
    n_bpoly_fq_divrem_series(q, r, t2, h, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_fq_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < h->length; i++)
        n_poly_fq_truncate(h->coeffs + i, p0, ctx);
    n_bpoly_fq_add(t2, r, h, ctx);

    n_bpoly_swap(G, t1);
    n_bpoly_swap(H, t2);

#if WANT_ASSERT
    {
        slong j, d = fq_nmod_ctx_degree(ctx);
        n_bpoly_fq_mul(t1, G, H, ctx);
        n_bpoly_fq_sub(c, f, t1, ctx);

        for (i = 0; i < c->length; i++)        
        for (j = 0; j < p0 + p1; j++)
        {
            FLINT_ASSERT(j >= c->coeffs[i].length ||
                         _n_fq_is_zero(c->coeffs[i].coeffs + d*j, d));
        }
    }
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
    const fq_nmod_ctx_t ctx)
{
    slong i;
    n_bpoly_t c, t1, t2, q, r;

    n_bpoly_init(c);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(q);
    n_bpoly_init(r);

    for (i = 0; i < b->length; i++)
        n_poly_fq_truncate(b->coeffs + i, p0, ctx);
    for (i = 0; i < a->length; i++)
        n_poly_fq_truncate(a->coeffs + i, p0, ctx);

    n_bpoly_fq_mul(t1, G, a, ctx);
    n_bpoly_fq_mul(t2, H, b, ctx);
    n_bpoly_fq_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        n_poly_fq_neg(c->coeffs + i, c->coeffs + i, ctx);
    n_poly_fq_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    n_bpoly_normalise(c);

    for (i = 0; i < c->length; i++)
    {
    #if WANT_ASSERT
        {
            slong j, d = fq_nmod_ctx_degree(ctx);
            for (j = 0; j < p0; j++)
            {
                FLINT_ASSERT(j >= c->coeffs[i].length ||
                             _n_fq_is_zero(c->coeffs[i].coeffs + d*j, d));
            }
        }
    #endif
        n_poly_fq_shift_right(c->coeffs + i, c->coeffs + i, p0, ctx);
        n_poly_fq_truncate(c->coeffs + i, p1, ctx);
    }

    n_bpoly_fq_mul_series(t1, c, b, p1, ctx);
    n_bpoly_fq_divrem_series(q, r, t1, G, p1, ctx);

    for (i = 0; i < r->length; i++)
        n_poly_fq_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);

    n_bpoly_fq_add(t1, r, b, ctx);

    n_bpoly_fq_mul_series(t2, c, a, p1, ctx);
    n_bpoly_fq_divrem_series(q, r, t2, H, p1, ctx);

    for (i = 0; i < r->length; i++)
        n_poly_fq_shift_left(r->coeffs + i, r->coeffs + i, p0, ctx);

    n_bpoly_fq_add(t2, r, a, ctx);

    n_bpoly_swap(t1, B);
    n_bpoly_swap(t2, A);

#if WANT_ASSERT
    n_bpoly_fq_mul(t1, G, A, ctx);
    n_bpoly_fq_mul(t2, H, B, ctx);
    n_bpoly_fq_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        n_poly_fq_neg(c->coeffs + i, c->coeffs + i, ctx);
    n_poly_fq_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    n_bpoly_normalise(c);

    {
        slong j, d = fq_nmod_ctx_degree(ctx);
        for (i = 0; i < c->length; i++)        
        for (j = 0; j < p0 + p1; j++)
        {
            FLINT_ASSERT(j >= c->coeffs[i].length ||
                         _n_fq_is_zero(c->coeffs[i].coeffs + d*j, d));
        }
    }
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
    const fq_nmod_ctx_t ctx)
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


int n_bpoly_fq_factor_smprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    int allow_shift,
    const fq_nmod_ctx_t ctx)
{
    int success;
    slong i, r;
    slong Blenx = B->length;
    slong Bleny;
    slong final_order, curr_lift_order, prev_lift_order, next_lift_order;
    slong * starts;
    fq_nmod_t alpha, Blc;
    fq_nmod_poly_t Beval;
    fq_nmod_poly_factor_t local_fac;
    n_bpoly_t monicB;
    nmod_mat_t N;
    slong * link;
    n_bpoly_struct * v, * w, ** lift_fac;
    n_tpoly_t tmp;
    slong e[FLINT_BITS];
    slong old_nrows;
    slong zas_limit;

    FLINT_ASSERT(Blenx > 1);

    fq_nmod_init(alpha, ctx);
    fq_nmod_init(Blc, ctx);
    fq_nmod_poly_init(Beval, ctx);
    fq_nmod_poly_factor_init(local_fac, ctx);

    n_bpoly_init(monicB);
    n_tpoly_init(tmp);
    nmod_mat_init(N, 0, 0, ctx->mod.n);
    starts = (slong *) flint_malloc(Blenx*sizeof(slong));
    link = (slong *) flint_malloc(sizeof(slong));
    lift_fac = FLINT_ARRAY_ALLOC(1, n_bpoly_struct *);

    /* init done */

    n_bpoly_fq_make_primitive(c, B, ctx);

    /* deg_y(B) + 1 */
    Bleny = 0;
    for (i = 0; i < B->length; i++)
        Bleny = FLINT_MAX(Bleny, (B->coeffs + i)->length);

    /* CLD bounds */
    for (i = 0; i < Blenx; i++)
        starts[i] = Bleny;

    fq_nmod_zero(alpha, ctx);
    goto got_alpha;

next_alpha:

    if (!allow_shift || !fq_nmod_next(alpha, ctx))
    {
        success = 0;
        goto cleanup;
    }

got_alpha:

    n_bpoly_fq_eval_var1(Beval, B, alpha, ctx);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blenx)
        goto next_alpha;

    local_fac->num = 0; /* stupid */
    fq_nmod_poly_factor(local_fac, Blc, Beval, ctx);

    r = local_fac->num;

    /* if multiple factors, get new alpha */
    for (i = 0; i < r; i++)
    {
        if (local_fac->exp[i] != 1)
            goto next_alpha;
    }

	/* if one factor, A is irreducible */
	if (r < 2)
	{
        n_tpoly_fit_length(F, 1);
        n_bpoly_swap(F->coeffs + 0, B);
        F->length = 1;
        success = 1;
        goto cleanup;
	}

    /* done if B is constant in y */
    if (Bleny < 2)
    {
        n_tpoly_fit_length(F, r);
        F->length = r;
        for (i = 0; i < r; i++)
            n_bpoly_fq_set_fq_nmod_poly_var0(F->coeffs + i, local_fac->poly + i, ctx);
        success = 1;
        goto cleanup;
    }

    /* precision for constructing true factors */
    final_order = Bleny;

    n_bpoly_fq_taylor_shift_var1_fq_nmod(B, B, alpha, ctx);

    nmod_mat_clear(N);
    nmod_mat_init(N, r, r, ctx->mod.n);
    for (i = 0; i < r; i++)
        nmod_mat_entry(N, i, i) = 1;

    link = (slong *) flint_realloc(link, (2*r - 2)*sizeof(slong));
    lift_fac = (n_bpoly_struct **) flint_realloc(lift_fac, r*sizeof(n_bpoly_struct *));

    n_tpoly_fit_length(tmp, 2*(2*r - 2));
    v = tmp->coeffs + 0;
    w = tmp->coeffs + (2*r - 2);

    curr_lift_order = final_order + r;

    n_bpoly_fq_make_monic_series(monicB, B, curr_lift_order, ctx);

    _hensel_build_tree(link, v, w, local_fac->poly, r, ctx);
    for (i = 0; i < 2*r - 2; i++)
        if (-link[i] - 1 >= 0)
            lift_fac[-link[i] - 1] = v + i;

    FLINT_ASSERT(curr_lift_order > 1);
    for (i = 0, e[i] = curr_lift_order; e[i] > 1; i++)
        e[i+1] = (e[i] + 1) / 2;

    for (i--; i > 0; i--)
        _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, e[i+1], e[i]-e[i+1], ctx);

    prev_lift_order = e[1];
    _hensel_lift_tree(1, link, v, w, monicB, 2*r-4, e[1], e[0]-e[1], ctx);

    zas_limit = 1;

try_zas:

    F->length = 0;
    fq_nmod_neg(alpha, alpha, ctx);
    success = _zassenhaus(zas_limit, F, alpha, N, lift_fac, r, final_order, B, ctx);
    fq_nmod_neg(alpha, alpha, ctx);
    if (success)
        goto cleanup;

    zas_limit = 2;

more:

    old_nrows = nmod_mat_nrows(N);
    _lattice(N, lift_fac, r, curr_lift_order, starts, B, ctx);
    if (nmod_mat_nrows(N) < old_nrows && nmod_mat_is_reduced(N))
        goto try_zas;

    next_lift_order = curr_lift_order + r;
    next_lift_order = FLINT_MIN(next_lift_order, 2*curr_lift_order);

    _hensel_lift_tree(-1, link, v, w, monicB, 2*r-4, prev_lift_order,
                                       curr_lift_order - prev_lift_order, ctx);

    n_bpoly_fq_make_monic_series(monicB, B, next_lift_order, ctx);

    _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, curr_lift_order,
                                       next_lift_order - curr_lift_order, ctx);

    prev_lift_order = curr_lift_order;
    curr_lift_order = next_lift_order;

    goto more;

cleanup:

    flint_free(starts);
    flint_free(link);
    flint_free(lift_fac);

    nmod_mat_clear(N);
    n_bpoly_clear(monicB);
    n_tpoly_clear(tmp);

    fq_nmod_poly_clear(Beval, ctx);
    fq_nmod_poly_factor_clear(local_fac, ctx);
    fq_nmod_clear(alpha, ctx);
    fq_nmod_clear(Blc, ctx);

    return success;
}

