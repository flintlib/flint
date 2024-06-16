/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "fmpz_poly_factor.h"
#include "fq_nmod.h"
#include "fq_nmod_poly.h"
#include "fq_nmod_poly_factor.h"
#include "n_poly.h"
#include "nmod_mpoly_factor.h"

static void n_bpoly_eval_fq_nmod_poly(
    fq_nmod_poly_t A,
    const fq_nmod_ctx_t ectx,
    const n_bpoly_t B)
{
    slong i;
    n_poly_t t;
    n_poly_t mock;
    nmod_poly_t mock2;

    n_poly_init(t);

    fq_nmod_poly_zero(A, ectx);
    for (i = B->length - 1; i >= 0; i--)
    {
        n_poly_mock(mock, ectx->modulus);
        n_poly_mod_rem(t, B->coeffs + i, mock, ectx->modulus->mod);
        nmod_poly_mock(mock2, t, ectx->modulus->mod);
        fq_nmod_poly_set_coeff(A, i, mock2, ectx);
    }

    n_poly_clear(t);
}


static void n_bpoly_mod_make_monic_mod(n_bpoly_t A, n_poly_t mk, nmod_t mod)
{
    slong i;
    n_poly_t t, lcinv;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));

    n_poly_init(t);
    n_poly_init(lcinv);
    if (!n_poly_mod_invmod(lcinv, A->coeffs + A->length - 1, mk, mod))
    {
        FLINT_ASSERT(0);
    }

    for (i = 0; i < A->length; i++)
    {
        n_poly_mod_mulmod(t, A->coeffs + i, lcinv, mk, mod);
        n_poly_swap(A->coeffs + i, t);
    }

    n_poly_clear(t);
    n_poly_clear(lcinv);
}


static void n_bpoly_set_fq_nmod_poly_gen0(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t FLINT_UNUSED(ectx))
{
    slong i;

    n_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        n_poly_set_nmod_poly(A->coeffs + i, B->coeffs + i);
}


static void n_bpoly_mod_mul_mod_poly(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const n_poly_t m,
    nmod_t ctx)
{
    slong i, j;
    n_poly_t t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_poly_init(t);

    n_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        n_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    for (j = 0; j < C->length; j++)
    {
        n_poly_mod_mul(t, B->coeffs + i, C->coeffs + j, ctx);
        n_poly_mod_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        n_poly_mod_rem(A->coeffs + i + j, A->coeffs + i + j, m, ctx);
    }

    A->length = B->length + C->length - 1;
    n_bpoly_normalise(A);

    n_poly_clear(t);
}

/* division in ((Z/nZ)[y]/m(y))[x] */
static void n_bpoly_mod_divrem_mod_poly(
    n_bpoly_t Q,
    n_bpoly_t R,
    const n_bpoly_t A,
    const n_bpoly_t B,
    const n_poly_t m,
    nmod_t ctx)
{
    slong i, qoff;
    n_poly_t q, t, Binv;

    FLINT_ASSERT(R != A);
    FLINT_ASSERT(R != B);
    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    n_poly_init(q);
    n_poly_init(t);
    n_poly_init(Binv);

    n_bpoly_set(R, A);

    Q->length = 0;

    if (!n_poly_mod_invmod(Binv, B->coeffs + B->length - 1, m, ctx))
    {
        FLINT_ASSERT(0);
    }

    while (R->length >= B->length)
    {
        n_poly_mod_mulmod(q, R->coeffs + R->length - 1, Binv, m, ctx);

        for (i = 0; i < B->length; i++)
        {
            n_poly_mod_mulmod(t, B->coeffs + i, q, m, ctx);
            n_poly_mod_sub(R->coeffs + i + R->length - B->length,
                           R->coeffs + i + R->length - B->length, t, ctx);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            n_bpoly_fit_length(Q, qoff + 1);
            for (i = Q->length; i <= qoff; i++)
                n_poly_zero(Q->coeffs + i);
            Q->length = qoff + 1;
        }

        n_poly_set(Q->coeffs + qoff, q);

        FLINT_ASSERT(n_poly_is_zero(R->coeffs + R->length - 1));

        n_bpoly_normalise(R);
    }

    n_poly_clear(q);
    n_poly_clear(t);
    n_poly_clear(Binv);
}

static int _zassenhaus(
    const zassenhaus_prune_t zas,
    slong limit,
    n_tpoly_t F,
    const n_poly_t finalmpow,
    const nmod_mat_t N,
    n_bpoly_struct * const * loc_fac_org,
    slong r,
    const n_bpoly_t B,
    nmod_t ctx)
{
    int success;
    slong total_deg;
    n_bpoly_t Q, R, t1, t2;
    n_poly_t leadf, g;
    slong i, j, k, len, d = nmod_mat_nrows(N);
    slong * subset;
    n_bpoly_struct * loc_fac;
    n_bpoly_struct * f;
    n_bpoly_t B_copy;

    FLINT_ASSERT(nmod_mat_ncols(N) == r);

    loc_fac = (n_bpoly_struct *) flint_malloc(d*sizeof(n_bpoly_struct));
    for (i = 0; i < d; i++)
        n_bpoly_init(loc_fac + i);

    n_poly_init(g);
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_poly_init(leadf);
    n_bpoly_init(B_copy);

    for (i = 0; i < d; i++)
    {
        n_bpoly_one(loc_fac + i);
        for (j = 0; j < r; j++)
        {
            if (nmod_mat_entry(N, i, j) == 0)
                continue;
            FLINT_ASSERT(nmod_mat_entry(N, i, j) == 1);
            n_bpoly_mod_mul_mod_poly(t1, loc_fac + i, loc_fac_org[j], finalmpow, ctx);
            n_bpoly_swap(t1, loc_fac + i);
        }
    }

    f = (n_bpoly_struct *) B;
    FLINT_ASSERT(f->length > 0);
    n_poly_set(leadf, f->coeffs + f->length - 1);

    subset = (slong *) flint_malloc(d * sizeof(slong));
    for (k = 0; k < d; k++)
        subset[k] = k;

    len = d;
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
                    total_deg += loc_fac[subset[i]].length - 1;
            }

            if (!zassenhaus_prune_degree_is_possible(zas, total_deg))
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
                continue;
            }

            n_bpoly_set_poly_gen1(t1, leadf);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    n_bpoly_mod_mul_mod_poly(t2, t1, loc_fac + subset[i], finalmpow, ctx);
                    n_bpoly_swap(t1, t2);
                }
            }

            n_bpoly_mod_make_primitive(g, t1, ctx);
            if (n_bpoly_mod_divides(Q, f, t1, ctx))
            {
                n_tpoly_fit_length(F, F->length + 1);
                n_bpoly_swap(F->coeffs + F->length, t1);
                F->length++;
                f = B_copy;
                n_bpoly_swap(f, Q);
                FLINT_ASSERT(f->length > 0);
                n_poly_set(leadf, f->coeffs + f->length - 1);

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
        n_bpoly_set(F->coeffs + F->length, f);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(f->length == 1);
        FLINT_ASSERT(n_poly_is_one(f->coeffs + 0));
    }

    success = 1;

cleanup:

    flint_free(subset);

    n_poly_clear(g);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_poly_clear(leadf);
    n_bpoly_clear(B_copy);

    for (i = 0; i < d; i++)
        n_bpoly_clear(loc_fac + i);
    flint_free(loc_fac);

    return success;
}


static void _hensel_build_tree(
    slong * link,
    n_bpoly_struct * v,
    n_bpoly_struct * w,
    const fq_nmod_poly_struct * local_facs,
    slong r,
    fq_nmod_ctx_t ctx)
{
    slong i, j;
    fq_nmod_poly_t d;
    fq_nmod_poly_struct * V;
    fq_nmod_poly_struct * W;

    V = (fq_nmod_poly_struct *) flint_malloc((2*r - 2)*sizeof(fq_nmod_poly_struct));
    W = (fq_nmod_poly_struct *) flint_malloc((2*r - 2)*sizeof(fq_nmod_poly_struct));

    fq_nmod_poly_init(d, ctx);
    for (i = 0; i < 2*r - 2; i++)
    {
        fq_nmod_poly_init(V + i, ctx);
        fq_nmod_poly_init(W + i, ctx);
    }

    for (i = 0; i < r; i++)
    {
        fq_nmod_poly_set(V + i, local_facs + i, ctx);
        link[i] = -i - 1;
    }

    for (i = r, j = 0; j < 2*r - 4; i++, j += 2)
    {
        slong s, minp, mind;

        minp = j;
        mind = fq_nmod_poly_degree(V + j, ctx);
        for (s = j + 1; s < i; s++)
        {
            if (fq_nmod_poly_degree(V + s, ctx) < mind)
            {
                minp = s;
                mind = fq_nmod_poly_degree(V + s, ctx);
            }
        }
        fq_nmod_poly_swap(V + j, V + minp, ctx);
        FLINT_SWAP(slong, link[j], link[minp]);

        minp = j + 1;
        mind = fq_nmod_poly_degree(V + j + 1, ctx);
        for (s = j + 2; s < i; s++)
        {
            if (fq_nmod_poly_degree(V + s, ctx) < mind)
            {
                minp = s;
                mind = fq_nmod_poly_degree(V + s, ctx);
            }
        }
        fq_nmod_poly_swap(V + j + 1, V + minp, ctx);
        FLINT_SWAP(slong, link[j + 1], link[minp]);

        fq_nmod_poly_mul(V + i, V + j, V + j + 1, ctx);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        fq_nmod_poly_xgcd(d, W + j, W + j + 1, V + j, V + j + 1, ctx);
        FLINT_ASSERT(fq_nmod_poly_is_one(d, ctx));
    }

    for (j = 0; j < 2*r - 2; j++)
    {
        n_bpoly_set_fq_nmod_poly_gen0(v + j, V + j, ctx);
        n_bpoly_set_fq_nmod_poly_gen0(w + j, W + j, ctx);
    }

    fq_nmod_poly_clear(d, ctx);
    for (i = 0; i < 2*r - 2; i++)
    {
        fq_nmod_poly_clear(V + i, ctx);
        fq_nmod_poly_clear(W + i, ctx);
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
    const n_poly_t p0,
    const n_poly_t p1,
    nmod_t ctx)
{
    slong i;
    n_bpoly_t c, t1, t2, q, r;
    n_poly_t tq, tr;

    n_bpoly_init(c);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(q);
    n_bpoly_init(r);
    n_poly_init(tq);
    n_poly_init(tr);

#if FLINT_WANT_ASSERT
    n_bpoly_mod_mul(t1, g, a, ctx);
    n_bpoly_mod_mul(t2, h, b, ctx);
    n_bpoly_mod_add(c, t1, t2, ctx);
    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        n_poly_mod_neg(c->coeffs + i, c->coeffs + i, ctx);
    n_poly_mod_add_ui(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    n_bpoly_normalise(c);
    for (i = 0; i < c->length; i++)
    {
        n_poly_mod_divrem(tq, tr, c->coeffs + i, p0, ctx);
        FLINT_ASSERT(n_poly_is_zero(tr));
    }
#endif

    n_bpoly_mod_mul(t1, g, h, ctx);
    n_bpoly_mod_sub(c, f, t1, ctx);
    for (i = 0; i < c->length; i++)
    {
        n_poly_mod_divrem(tq, tr, c->coeffs + i, p0, ctx);
        FLINT_ASSERT(n_poly_is_zero(tr));
        n_poly_mod_divrem(tr, c->coeffs + i, tq, p1, ctx);
    }

    n_bpoly_mod_mul_mod_poly(t1, c, b, p1, ctx);
    n_bpoly_mod_divrem_mod_poly(q, r, t1, g, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_mod_mul(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < g->length; i++)
        n_poly_mod_divrem(tq, g->coeffs + i, g->coeffs + i, p0, ctx);
    n_bpoly_mod_add(t1, r, g, ctx);

    n_bpoly_mod_mul_mod_poly(t2, c, a, p1, ctx);
    n_bpoly_mod_divrem_mod_poly(q, r, t2, h, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_mod_mul(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < h->length; i++)
        n_poly_mod_divrem(tq, h->coeffs + i, h->coeffs + i, p0, ctx);
    n_bpoly_mod_add(t2, r, h, ctx);

    n_bpoly_swap(G, t1);
    n_bpoly_swap(H, t2);

#if FLINT_WANT_ASSERT
    {
        n_poly_t p01;
        n_poly_init(p01);
        n_poly_mod_mul(p01, p0, p1, ctx);
        n_bpoly_mod_mul(t1, G, H, ctx);
        n_bpoly_mod_sub(c, f, t1, ctx);
        for (i = 0; i < c->length; i++)
        {
            n_poly_mod_divrem(tq, tr, c->coeffs + i, p01, ctx);
            FLINT_ASSERT(n_poly_is_zero(tr));
        }
        n_poly_clear(p01);
    }
#endif

    n_bpoly_clear(c);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_bpoly_clear(q);
    n_bpoly_clear(r);
    n_poly_clear(tq);
    n_poly_clear(tr);
}

static void _hensel_lift_inv(
    n_bpoly_t A,
    n_bpoly_t B,
    const n_bpoly_t G,
    const n_bpoly_t H,
    n_bpoly_t a,
    n_bpoly_t b,
    const n_poly_t p0,
    const n_poly_t p1,
    nmod_t ctx)
{
    slong i;
    n_bpoly_t c, t1, t2, q, r;
    n_poly_t tq, tr;

    n_bpoly_init(c);
    n_bpoly_init(t1);
    n_bpoly_init(t2);
    n_bpoly_init(q);
    n_bpoly_init(r);
    n_poly_init(tq);
    n_poly_init(tr);

    for (i = 0; i < b->length; i++)
        n_poly_mod_divrem(tq, b->coeffs + i, b->coeffs + i, p0, ctx);
    for (i = 0; i < a->length; i++)
        n_poly_mod_divrem(tq, a->coeffs + i, a->coeffs + i, p0, ctx);

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
        n_poly_mod_divrem(tq, tr, c->coeffs + i, p0, ctx);
        FLINT_ASSERT(n_poly_is_zero(tr));
        n_poly_mod_divrem(tr, c->coeffs + i, tq, p1, ctx);
    }

    n_bpoly_mod_mul_mod_poly(t1, c, b, p1, ctx);
    n_bpoly_mod_divrem_mod_poly(q, r, t1, G, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_mod_mul(r->coeffs + i, r->coeffs + i, p0, ctx);
    n_bpoly_mod_add(t1, r, b, ctx);

    n_bpoly_mod_mul_mod_poly(t2, c, a, p1, ctx);
    n_bpoly_mod_divrem_mod_poly(q, r, t2, H, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_poly_mod_mul(r->coeffs + i, r->coeffs + i, p0, ctx);
    n_bpoly_mod_add(t2, r, a, ctx);

    n_bpoly_swap(t1, B);
    n_bpoly_swap(t2, A);

#if FLINT_WANT_ASSERT
    {
        n_poly_t p01;
        n_poly_init(p01);
        n_poly_mod_mul(p01, p0, p1, ctx);
        n_bpoly_mod_mul(t1, G, A, ctx);
        n_bpoly_mod_mul(t2, H, B, ctx);
        n_bpoly_mod_add(c, t1, t2, ctx);
        FLINT_ASSERT(c->length > 0);
        for (i = 0; i < c->length; i++)
            n_poly_mod_neg(c->coeffs + i, c->coeffs + i, ctx);
        n_poly_mod_add_ui(c->coeffs + 0, c->coeffs + 0, 1, ctx);
        n_bpoly_normalise(c);
        for (i = 0; i < c->length; i++)
        {
            n_poly_mod_divrem(tq, tr, c->coeffs + i, p01, ctx);
            FLINT_ASSERT(n_poly_is_zero(tr));
        }
        n_poly_clear(p01);
    }
#endif

    n_bpoly_clear(c);
    n_bpoly_clear(t1);
    n_bpoly_clear(t2);
    n_bpoly_clear(q);
    n_bpoly_clear(r);
    n_poly_clear(tq);
    n_poly_clear(tr);
}

static void _hensel_lift_tree(
    int opt,
    slong * link,
    n_bpoly_struct * v,
    n_bpoly_struct * w,
    const n_bpoly_t f,
    slong j,
    const n_poly_t p0,
    const n_poly_t p1,
    nmod_t ctx)
{
    FLINT_ASSERT(p1->length <= p0->length);

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


static void _lattice(
    nmod_mat_t N,
    n_bpoly_struct * const * g,
    slong r,
    const n_poly_t lift_alpha_pow,
    slong * starts,
    const n_bpoly_t f,
    nmod_t ctx)
{
    slong i, j, k;
    n_bpoly_t Q, R, dg;
    n_bpoly_struct * ld;
    nmod_mat_t M, T1, T2;
    ulong * trow;
    slong lift_order = lift_alpha_pow->length - 1;

    const dot_params_t params = _nmod_vec_dot_params(r, ctx);
    trow = (ulong *) flint_malloc(r*sizeof(ulong));
    n_bpoly_init(Q);
    n_bpoly_init(R);
    n_bpoly_init(dg);
    ld = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
        n_bpoly_init(ld + i);

    /* init done */

    for (i = 0; i < r; i++)
    {
        n_bpoly_mod_divrem_mod_poly(Q, R, f, g[i], lift_alpha_pow, ctx);
        FLINT_ASSERT(R->length == 0);
        n_bpoly_mod_derivative_gen0(R, g[i], ctx);
        n_bpoly_mod_mul_mod_poly(ld + i, Q, R, lift_alpha_pow, ctx);
    }

    for (k = 0; k + 1 < f->length; k++)
    {
        slong d = nmod_mat_nrows(N);

        if (d < 2)
            break;

        if (lift_order <= starts[k])
            continue;

        nmod_mat_init(M, lift_order - starts[k], d, ctx.n);

        for (j = starts[k]; j < lift_order; j++)
        {
            for (i = 0; i < r; i++)
                trow[i] = n_bpoly_get_coeff(ld + i, k, j);

            for (i = 0; i < d; i++)
                nmod_mat_entry(M, j - starts[k], i) =
                             _nmod_vec_dot(trow, N->rows[i], r, ctx, params);
        }

        nmod_mat_init_nullspace_tr(T1, M);

        nmod_mat_init(T2, nmod_mat_nrows(T1), nmod_mat_ncols(N), ctx.n);
        nmod_mat_mul(T2, T1, N);
        nmod_mat_swap(T2, N);
        nmod_mat_rref(N);

        nmod_mat_clear(M);
        nmod_mat_clear(T1);
        nmod_mat_clear(T2);

        if (nmod_mat_is_reduced(N))
            break;
    }

    flint_free(trow);
    n_bpoly_clear(Q);
    n_bpoly_clear(R);
    n_bpoly_clear(dg);
    for (i = 0; i < r; i++)
        n_bpoly_clear(ld + i);
    flint_free(ld);
}

void n_bpoly_mod_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    nmod_t ctx)
{
    int success;
    slong i, r, deg;
    slong Blenx = B->length;
    slong Bleny;
    slong final_pow, curr_lift_pow, prev_lift_pow, next_lift_pow;
    slong * starts;
    fq_nmod_poly_t Beval;
    fq_nmod_poly_factor_t local_fac;
    fq_nmod_t Blc;
    n_bpoly_t monicB;
    nmod_mat_t N;
    slong * link;
    n_bpoly_struct * v, * w, ** lift_fac;
    n_tpoly_t tmp;
    slong e[FLINT_BITS];
    slong old_nrows;
    slong zas_limit;
    n_poly_t final_alpha_pow, curr_alpha_pow, prev_alpha_pow, next_alpha_pow;
    n_poly_t alpha, p1;
    fq_nmod_ctx_t ectx;
    zassenhaus_prune_t zas;

    FLINT_ASSERT(Blenx > 1);

    deg = 2;
    fq_nmod_ctx_init_ui(ectx, ctx.n, deg, "y");
    n_poly_init(final_alpha_pow);
    n_poly_init(curr_alpha_pow);
    n_poly_init(prev_alpha_pow);
    n_poly_init(next_alpha_pow);
    fq_nmod_poly_init(Beval, ectx);
    fq_nmod_poly_factor_init(local_fac, ectx);
    fq_nmod_init(Blc, ectx);
    n_bpoly_init(monicB);
    n_tpoly_init(tmp);
    nmod_mat_init(N, 0, 0, ctx.n);
    starts = (slong *) flint_malloc(Blenx*sizeof(slong));
    link = (slong *) flint_malloc(sizeof(slong));
    lift_fac = (n_bpoly_struct **) flint_malloc(sizeof(n_bpoly_struct *));
    n_poly_init(p1);
    zassenhaus_prune_init(zas);

    /* init done */

    n_poly_mock(alpha, ectx->modulus);

    n_bpoly_mod_make_primitive(c, B, ctx);

    Bleny = 0;
    for (i = 0; i < B->length; i++)
        Bleny = FLINT_MAX(Bleny, (B->coeffs + i)->length);

    /* CLD bounds */
    for (i = 0; i < Blenx; i++)
        starts[i] = Bleny;

    zassenhaus_prune_set_degree(zas, Blenx - 1);

    goto got_alpha;

next_alpha:

    deg++;

	fq_nmod_ctx_clear(ectx);
	fq_nmod_ctx_init_ui(ectx, ctx.n, deg, "y");

    n_poly_mock(alpha, ectx->modulus);

got_alpha:

    n_bpoly_eval_fq_nmod_poly(Beval, ectx, B);

    /* if killed leading/trailing coeff, get new alpha */
    if (Beval->length != Blenx || fq_nmod_is_zero(Beval->coeffs + 0, ectx))
        goto next_alpha;

    local_fac->num = 0;
    fq_nmod_poly_factor(local_fac, Blc, Beval, ectx);

    r = local_fac->num;

    zassenhaus_prune_start_add_factors(zas);
    for (i = 0; i < r; i++)
        zassenhaus_prune_add_factor(zas,
            fq_nmod_poly_degree(local_fac->poly + i, ectx), local_fac->exp[i]);
    zassenhaus_prune_end_add_factors(zas);

    if ((r < 2 && local_fac->exp[0] == 1) ||
        zassenhaus_prune_must_be_irreducible(zas))
    {
        n_tpoly_fit_length(F, 1);
        F->length = 1;
        n_bpoly_swap(F->coeffs + 0, B);
        goto cleanup;
    }

    /* if multiple factors, get new alpha */
    for (i = 0; i < r; i++)
    {
        if (local_fac->exp[i] != 1)
            goto next_alpha;
    }

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(local_fac->poly[i].length > 1);
        FLINT_ASSERT(fq_nmod_is_one(local_fac->poly[i].coeffs + local_fac->poly[i].length - 1, ectx));
    }

    /* precision for constructing true factors */
    final_pow = (Bleny - 1 + deg)/deg;
    n_poly_mod_pow(final_alpha_pow, alpha, final_pow, ctx);

    nmod_mat_clear(N);
    nmod_mat_init(N, r, r, ctx.n);
    for (i = 0; i < r; i++)
        nmod_mat_entry(N, i, i) = 1;

    link = (slong *) flint_realloc(link, (2*r - 2)*sizeof(slong));
    lift_fac = (n_bpoly_struct **) flint_realloc(lift_fac,
                                                   r*sizeof(n_bpoly_struct *));

    n_tpoly_fit_length(tmp, 2*(2*r - 2));
    v = tmp->coeffs + 0;
    w = tmp->coeffs + (2*r - 2);

    curr_lift_pow = final_pow + r;
    n_poly_mod_pow(curr_alpha_pow, alpha, curr_lift_pow, ctx);

    n_bpoly_set(monicB, B);
    n_bpoly_mod_make_monic_mod(monicB, curr_alpha_pow, ctx);

    _hensel_build_tree(link, v, w, local_fac->poly, r, ectx);
    for (i = 0; i < 2*r - 2; i++)
        if (-link[i] - 1 >= 0)
            lift_fac[-link[i] - 1] = v + i;

    FLINT_ASSERT(curr_lift_pow > 1);
    for (i = 0, e[i] = curr_lift_pow; e[i] > 1; i++)
        e[i+1] = (e[i] + 1) / 2;

    for (i--; i > 0; i--)
    {
        n_poly_mod_pow(prev_alpha_pow, alpha, e[i+1], ctx);
        n_poly_mod_pow(p1, alpha, e[i]-e[i+1], ctx);
        _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, prev_alpha_pow, p1, ctx);
    }

    prev_lift_pow = e[1];
    n_poly_mod_pow(prev_alpha_pow, alpha, prev_lift_pow, ctx);
    n_poly_mod_pow(p1, alpha, curr_lift_pow - prev_lift_pow, ctx);
    _hensel_lift_tree(1, link, v, w, monicB, 2*r-4, prev_alpha_pow, p1, ctx);

    zas_limit = 2;

try_zas:

    F->length = 0;
    success = _zassenhaus(zas, zas_limit, F, final_alpha_pow, N, lift_fac, r, B, ctx);
    if (success)
        goto cleanup;

    zas_limit = 3;

more:

    old_nrows = nmod_mat_nrows(N);
    _lattice(N, lift_fac, r, curr_alpha_pow, starts, B, ctx);
    if (nmod_mat_nrows(N) < old_nrows && nmod_mat_is_reduced(N))
        goto try_zas;

    next_lift_pow = curr_lift_pow + r;
    next_lift_pow = FLINT_MIN(next_lift_pow, 2*curr_lift_pow);

    n_poly_mod_pow(p1, alpha, curr_lift_pow - prev_lift_pow, ctx);
    _hensel_lift_tree(-1, link, v, w, monicB, 2*r-4, prev_alpha_pow, p1, ctx);

    n_poly_mod_pow(p1, alpha, next_lift_pow - curr_lift_pow, ctx);

    n_poly_mod_mul(next_alpha_pow, next_alpha_pow, p1, ctx);
    n_bpoly_set(monicB, B);
    n_bpoly_mod_make_monic_mod(monicB, next_alpha_pow, ctx);

    _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, curr_alpha_pow, p1, ctx);

    prev_lift_pow = curr_lift_pow;
    curr_lift_pow = next_lift_pow;
    n_poly_swap(prev_alpha_pow, curr_alpha_pow);
    n_poly_swap(curr_alpha_pow, next_alpha_pow);

    goto more;

cleanup:

    n_poly_clear(final_alpha_pow);
    n_poly_clear(curr_alpha_pow);
    n_poly_clear(prev_alpha_pow);
    n_poly_clear(next_alpha_pow);
    fq_nmod_poly_clear(Beval, ectx);
    fq_nmod_poly_factor_clear(local_fac, ectx);
    fq_nmod_clear(Blc, ectx);
    n_bpoly_clear(monicB);
    n_tpoly_clear(tmp);
    nmod_mat_clear(N);
    flint_free(starts);
    flint_free(link);
    flint_free(lift_fac);
    n_poly_clear(p1);

    fq_nmod_ctx_clear(ectx);

    zassenhaus_prune_clear(zas);

    return;
}
