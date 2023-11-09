/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "fmpz_poly_factor.h"
#include "nmod_mpoly_factor.h"
#include "fq_nmod_poly_factor.h"
#include "fq_nmod_mpoly_factor.h"

static void n_fq_bpoly_eval_sm_to_lg(
    fq_nmod_poly_t E,
    const n_bpoly_t B,
    const bad_fq_nmod_embed_t emb)
{
    slong i;
    fq_nmod_poly_fit_length(E, B->length, emb->lgctx);
    for (i = 0; i < B->length; i++)
        bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(E->coeffs + i, B->coeffs + i, emb);
    E->length = B->length;
    _fq_nmod_poly_normalise(E, emb->lgctx);
}


static void n_fq_bpoly_set_n_fq_poly_gen0_lg_to_sm(
    n_bpoly_t A,
    const n_poly_t E,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(emb->lgctx);
    slong i;
    n_bpoly_fit_length(A, E->length);
    for (i = 0; i < E->length; i++)
        bad_n_fq_embed_lg_to_sm(A->coeffs + i, E->coeffs + lgd*i, emb);
    A->length = E->length;
    n_bpoly_normalise(A);
}

static void n_fq_bpoly_make_monic_mod_poly(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_poly_t m,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    n_poly_t lcinv, t, g;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(n_fq_bpoly_is_canonical(B, ctx));

    n_poly_init(lcinv);
    n_poly_init(g);
    n_poly_init(t);

    n_fq_poly_xgcd(g, lcinv, t, B->coeffs + B->length - 1, m, ctx);
    FLINT_ASSERT(n_fq_poly_is_one(g, ctx));

    n_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
        n_fq_poly_mulmod(A->coeffs + i, B->coeffs + i, lcinv, m, ctx);

    A->length = B->length;
    n_bpoly_normalise(A);

    n_poly_clear(lcinv);
    n_poly_clear(t);
    n_poly_clear(g);
}

/* multiplication in (Fq[y]/m(y))[x], inputs need not be reduced */
static void n_fq_bpoly_mul_mod_poly(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const n_poly_t m,
    const fq_nmod_ctx_t ctx)
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
        n_fq_poly_mul(t, B->coeffs + i, C->coeffs + j, ctx);
        n_fq_poly_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        n_fq_poly_rem(A->coeffs + i + j, A->coeffs + i + j, m, ctx);
    }

    A->length = B->length + C->length - 1;
    n_bpoly_normalise(A);

    n_poly_clear(t);
}

/* division in (Fq[y]/m(y))[x], inputs need not be reduced */
static void n_fq_bpoly_divrem_mod_poly(
    n_bpoly_t Q,
    n_bpoly_t R,
    const n_bpoly_t A,
    const n_bpoly_t B,
    const n_poly_t m,
    const fq_nmod_ctx_t ctx)
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

    n_fq_bpoly_set(R, A, ctx);

    Q->length = 0;

    n_fq_poly_xgcd(q, Binv, t, B->coeffs + B->length - 1, m, ctx);
    FLINT_ASSERT(n_fq_poly_is_one(q, ctx));

    while (R->length >= B->length)
    {
        n_fq_poly_mulmod(q, R->coeffs + R->length - 1, Binv, m, ctx);

        for (i = 0; i < B->length; i++)
        {
            n_fq_poly_mulmod(t, B->coeffs + i, q, m, ctx);
            n_fq_poly_sub(R->coeffs + i + R->length - B->length,
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

        n_fq_poly_set(Q->coeffs + qoff, q, ctx);

        FLINT_ASSERT(n_poly_is_zero(R->coeffs + R->length - 1));
        n_bpoly_normalise(R);
    }

    n_poly_clear(q);
    n_poly_clear(t);
    n_poly_clear(Binv);
}


static void _lattice(
    nmod_mat_t N,
    n_bpoly_struct * const * g,
    slong r,
    const n_poly_t lift_alpha_pow,
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
    slong lift_order = lift_alpha_pow->length - 1;

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
        n_fq_bpoly_divrem_mod_poly(Q, R, f, g[i], lift_alpha_pow, ctx);
        FLINT_ASSERT(R->length == 0);
        n_fq_bpoly_derivative_gen0(R, g[i], ctx);
        n_fq_bpoly_mul_mod_poly(ld + i, Q, R, lift_alpha_pow, ctx);
    }

    for (k = 0; k + 1 < f->length; k++)
    {
        slong d = nmod_mat_nrows(N);

        if (d < 2)
            break;

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
                         _nmod_vec_dot(trow, N->rows[i], r, ctx->mod, nlimbs);
        }

        nmod_mat_init_nullspace_tr(T1, M);

        nmod_mat_init(T2, nmod_mat_nrows(T1), nmod_mat_ncols(N), ctx->mod.n);
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
    const n_poly_t final_alpha_pow,
    const nmod_mat_t N,
    n_bpoly_struct * const * loc_fac_org,
    slong r,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    int success;
    n_bpoly_t Q, R, t1, t2;
    n_poly_t g;
    slong i, j, s, len, d = nmod_mat_nrows(N);
    slong * subset;
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
        n_fq_bpoly_one(loc_fac + i, ctx);
        for (j = 0; j < r; j++)
        {
            if (nmod_mat_entry(N, i, j) == 0)
                continue;
            FLINT_ASSERT(nmod_mat_entry(N, i, j) == 1);
            n_fq_bpoly_mul_mod_poly(t1, loc_fac + i, loc_fac_org[j],
                                                         final_alpha_pow, ctx);
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
            n_fq_bpoly_set_n_fq_poly_gen1(t1, f->coeffs + f->length - 1, ctx);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    n_fq_bpoly_mul_mod_poly(t2, t1, loc_fac + subset[i],
                                                         final_alpha_pow, ctx);
                    n_bpoly_swap(t1, t2);
                }
            }

            n_fq_bpoly_make_primitive(g, t1, ctx);
            if (n_fq_bpoly_divides(Q, f, t1, ctx))
            {
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
        n_fq_bpoly_set(F->coeffs + F->length, f, ctx);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(f->length == 1);
        FLINT_ASSERT(n_fq_poly_is_one(f->coeffs + 0, ctx));
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
    const bad_fq_nmod_embed_t emb)
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
        n_fq_poly_set_fq_nmod_poly(V + i, local_facs + i, emb->lgctx);
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
        FLINT_SWAP(slong, link[j], link[minp]);

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
        FLINT_SWAP(slong, link[j + 1], link[minp]);

        n_fq_poly_mul(V + i, V + j, V + j + 1, emb->lgctx);
        link[i] = j;
    }

    for (j = 0; j < 2*r - 2; j += 2)
    {
        n_fq_poly_xgcd(d, W + j, W + j + 1, V + j, V + j + 1, emb->lgctx);
        FLINT_ASSERT(n_fq_poly_is_one(d, emb->lgctx));
    }

    for (j = 0; j < 2*r - 2; j++)
    {
        n_fq_bpoly_set_n_fq_poly_gen0_lg_to_sm(v + j, V + j, emb);
        n_fq_bpoly_set_n_fq_poly_gen0_lg_to_sm(w + j, W + j, emb);
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
    const n_poly_t p0,
    const n_poly_t p1,
    const fq_nmod_ctx_t ctx)
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

#ifdef FLINT_WANT_ASSERT
    n_fq_bpoly_mul(t1, g, a, ctx);
    n_fq_bpoly_mul(t2, h, b, ctx);
    n_fq_bpoly_add(c, t1, t2, ctx);
    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        n_fq_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
    n_fq_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    n_bpoly_normalise(c);
    for (i = 0; i < c->length; i++)
    {
        n_fq_poly_divrem(tq, tr, c->coeffs + i, p0, ctx);
        FLINT_ASSERT(n_fq_poly_is_canonical(tr, ctx));
        FLINT_ASSERT(n_poly_is_zero(tr));
    }
#endif

    n_fq_bpoly_mul(t1, g, h, ctx);
    n_fq_bpoly_sub(c, f, t1, ctx);
    for (i = 0; i < c->length; i++)
    {
        n_fq_poly_divrem(tq, tr, c->coeffs + i, p0, ctx);
        FLINT_ASSERT(n_fq_poly_is_canonical(tr, ctx));
        FLINT_ASSERT(n_poly_is_zero(tr));
        n_fq_poly_divrem(tr, c->coeffs + i, tq, p1, ctx);
    }

    n_fq_bpoly_mul_mod_poly(t1, c, b, p1, ctx);
    n_fq_bpoly_divrem_mod_poly(q, r, t1, g, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_fq_poly_mul(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < g->length; i++)
        n_fq_poly_divrem(tq, g->coeffs + i, g->coeffs + i, p0, ctx);
    n_fq_bpoly_add(t1, r, g, ctx);

    n_fq_bpoly_mul_mod_poly(t2, c, a, p1, ctx);
    n_fq_bpoly_divrem_mod_poly(q, r, t2, h, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_fq_poly_mul(r->coeffs + i, r->coeffs + i, p0, ctx);
    for (i = 0; i < h->length; i++)
        n_fq_poly_divrem(tq, h->coeffs + i, h->coeffs + i, p0, ctx);
    n_fq_bpoly_add(t2, r, h, ctx);

    n_bpoly_swap(G, t1);
    n_bpoly_swap(H, t2);

#ifdef FLINT_WANT_ASSERT
    {
        n_poly_t p01;
        n_poly_init(p01);
        n_fq_poly_mul(p01, p0, p1, ctx);
        n_fq_bpoly_mul(t1, G, H, ctx);
        n_fq_bpoly_sub(c, f, t1, ctx);
        for (i = 0; i < c->length; i++)
        {
            n_fq_poly_divrem(tq, tr, c->coeffs + i, p01, ctx);
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
    const fq_nmod_ctx_t ctx)
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

#ifdef FLINT_WANT_ASSERT
    {
        n_fq_bpoly_mul(t1, G, A, ctx);
        n_fq_bpoly_mul(t2, H, B, ctx);
        n_fq_bpoly_add(c, t1, t2, ctx);
        FLINT_ASSERT(c->length > 0);
        for (i = 0; i < c->length; i++)
            n_fq_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
        n_fq_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
        n_bpoly_normalise(c);
        for (i = 0; i < c->length; i++)
        {
            n_fq_poly_divrem(tq, tr, c->coeffs + i, p0, ctx);
            FLINT_ASSERT(n_poly_is_zero(tr));
        }
    }
#endif

    for (i = 0; i < b->length; i++)
        n_fq_poly_divrem(tq, b->coeffs + i, b->coeffs + i, p0, ctx);
    for (i = 0; i < a->length; i++)
        n_fq_poly_divrem(tq, a->coeffs + i, a->coeffs + i, p0, ctx);

    n_fq_bpoly_mul(t1, G, a, ctx);
    n_fq_bpoly_mul(t2, H, b, ctx);
    n_fq_bpoly_add(c, t1, t2, ctx);

    FLINT_ASSERT(c->length > 0);
    for (i = 0; i < c->length; i++)
        n_fq_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
    n_fq_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
    n_bpoly_normalise(c);

    for (i = 0; i < c->length; i++)
    {
        n_fq_poly_divrem(tq, tr, c->coeffs + i, p0, ctx);
        FLINT_ASSERT(n_poly_is_zero(tr));
        n_fq_poly_divrem(tr, c->coeffs + i, tq, p1, ctx);
    }

    n_fq_bpoly_mul_mod_poly(t1, c, b, p1, ctx);
    n_fq_bpoly_divrem_mod_poly(q, r, t1, G, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_fq_poly_mul(r->coeffs + i, r->coeffs + i, p0, ctx);
    n_fq_bpoly_add(t1, r, b, ctx);

    n_fq_bpoly_mul_mod_poly(t2, c, a, p1, ctx);
    n_fq_bpoly_divrem_mod_poly(q, r, t2, H, p1, ctx);
    for (i = 0; i < r->length; i++)
        n_fq_poly_mul(r->coeffs + i, r->coeffs + i, p0, ctx);
    n_fq_bpoly_add(t2, r, a, ctx);

    n_bpoly_swap(t1, B);
    n_bpoly_swap(t2, A);

#ifdef FLINT_WANT_ASSERT
    {
        n_poly_t p01;
        n_poly_init(p01);
        n_fq_poly_mul(p01, p0, p1, ctx);
        n_fq_bpoly_mul(t1, G, A, ctx);
        n_fq_bpoly_mul(t2, H, B, ctx);
        n_fq_bpoly_add(c, t1, t2, ctx);
        FLINT_ASSERT(c->length > 0);
        for (i = 0; i < c->length; i++)
            n_fq_poly_neg(c->coeffs + i, c->coeffs + i, ctx);
        n_fq_poly_add_si(c->coeffs + 0, c->coeffs + 0, 1, ctx);
        n_bpoly_normalise(c);
        for (i = 0; i < c->length; i++)
        {
            n_fq_poly_divrem(tq, tr, c->coeffs + i, p01, ctx);
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
    const fq_nmod_ctx_t ctx)
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


int n_fq_bpoly_factor_lgprime(
    n_poly_t c,
    n_tpoly_t F,
    n_bpoly_t B,
    const fq_nmod_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    slong i, r, edeg;
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
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx_mock, ctx_mock;

    FLINT_ASSERT(Blenx > 1);

    /* TODO remove mpoly_ctx_t from bad_fq_nmod_mpoly_embed_chooser interface */
    mpoly_ctx_init(ctx_mock->minfo, 2, ORD_LEX);
    *ctx_mock->fqctx = *ctx;

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, ectx_mock, ctx_mock, state);
    n_poly_init(final_alpha_pow);
    n_poly_init(curr_alpha_pow);
    n_poly_init(prev_alpha_pow);
    n_poly_init(next_alpha_pow);
    n_poly_init(p1);
    n_poly_init(alpha);
    n_bpoly_init(monicB);
    n_tpoly_init(tmp);

    fq_nmod_poly_init(Beval, ectx_mock->fqctx);
    fq_nmod_poly_factor_init(local_fac, ectx_mock->fqctx);
    fq_nmod_init(Blc, ectx_mock->fqctx);

    nmod_mat_init(N, 0, 0, ctx->mod.n);
    starts = (slong *) flint_malloc(Blenx*sizeof(slong));
    link = (slong *) flint_malloc(sizeof(slong));
    lift_fac = FLINT_ARRAY_ALLOC(1, n_bpoly_struct *);

    /* init done */

    n_fq_poly_set_fq_nmod_poly(alpha, cur_emb->h, ctx);
    edeg = cur_emb->h->length - 1;

    n_fq_bpoly_make_primitive(c, B, ctx);

    /* deg_y(B) + 1 */
    Bleny = 0;
    for (i = 0; i < B->length; i++)
        Bleny = FLINT_MAX(Bleny, (B->coeffs + i)->length);

    /* CLD bounds */
    for (i = 0; i < Blenx; i++)
        starts[i] = Bleny;

    goto got_alpha;

next_alpha:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx_mock, ctx_mock, state);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    n_fq_poly_set_fq_nmod_poly(alpha, cur_emb->h, ctx);
    edeg = cur_emb->h->length - 1;

got_alpha:

    n_fq_bpoly_eval_sm_to_lg(Beval, B, cur_emb);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blenx)
        goto next_alpha;

    fq_nmod_poly_factor(local_fac, Blc, Beval, ectx_mock->fqctx);

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

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(local_fac->poly[i].length > 1);
        FLINT_ASSERT(fq_nmod_is_one(local_fac->poly[i].coeffs + local_fac->poly[i].length - 1, ectx_mock->fqctx));
    }

    /* precision for constructing true factors */
    final_pow = (Bleny + edeg - 1)/edeg;
    n_fq_poly_pow(final_alpha_pow, alpha, final_pow, ctx);

    nmod_mat_clear(N);
    nmod_mat_init(N, r, r, ctx->mod.n);
    for (i = 0; i < r; i++)
        nmod_mat_entry(N, i, i) = 1;

    link = (slong *) flint_realloc(link, (2*r - 2)*sizeof(slong));
    lift_fac = (n_bpoly_struct **) flint_realloc(lift_fac, r*sizeof(n_bpoly_struct *));

    n_tpoly_fit_length(tmp, 2*(2*r - 2));
    v = tmp->coeffs + 0;
    w = tmp->coeffs + (2*r - 2);

    curr_lift_pow = final_pow + r;
    n_fq_poly_pow(curr_alpha_pow, alpha, curr_lift_pow, ctx);

    n_fq_bpoly_make_monic_mod_poly(monicB, B, curr_alpha_pow, ctx);

    _hensel_build_tree(link, v, w, local_fac->poly, r, cur_emb);
    for (i = 0; i < 2*r - 2; i++)
        if (-link[i] - 1 >= 0)
            lift_fac[-link[i] - 1] = v + i;

    FLINT_ASSERT(curr_lift_pow > 1);
    for (i = 0, e[i] = curr_lift_pow; e[i] > 1; i++)
        e[i+1] = (e[i] + 1) / 2;

    for (i--; i > 0; i--)
    {
        n_fq_poly_pow(prev_alpha_pow, alpha, e[i+1], ctx);
        n_fq_poly_pow(p1, alpha, e[i]-e[i+1], ctx);
        _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, prev_alpha_pow, p1, ctx);
    }

    prev_lift_pow = e[1];
    n_fq_poly_pow(prev_alpha_pow, alpha, prev_lift_pow, ctx);
    n_fq_poly_pow(p1, alpha, curr_lift_pow - prev_lift_pow, ctx);
    _hensel_lift_tree(1, link, v, w, monicB, 2*r-4, prev_alpha_pow, p1, ctx);

    zas_limit = 1;

try_zas:

    F->length = 0;
    success = _zassenhaus(zas_limit, F, final_alpha_pow, N, lift_fac, r, B, ctx);
    if (success)
        goto cleanup;

    zas_limit = 2;

more:

    old_nrows = nmod_mat_nrows(N);
    _lattice(N, lift_fac, r, curr_alpha_pow, starts, B, ctx);
    if (nmod_mat_nrows(N) < old_nrows && nmod_mat_is_reduced(N))
        goto try_zas;

    next_lift_pow = curr_lift_pow + r;
    next_lift_pow = FLINT_MIN(next_lift_pow, 2*curr_lift_pow);

    n_fq_poly_pow(p1, alpha, curr_lift_pow - prev_lift_pow, ctx);
    _hensel_lift_tree(-1, link, v, w, monicB, 2*r-4, prev_alpha_pow, p1, ctx);

    n_fq_poly_pow(p1, alpha, next_lift_pow - curr_lift_pow, ctx);

    n_fq_poly_mul(next_alpha_pow, next_alpha_pow, p1, ctx);
    n_fq_bpoly_make_monic_mod_poly(monicB, B, next_alpha_pow, ctx);

    _hensel_lift_tree(0, link, v, w, monicB, 2*r-4, curr_alpha_pow, p1, ctx);

    prev_lift_pow = curr_lift_pow;
    curr_lift_pow = next_lift_pow;
    n_poly_swap(prev_alpha_pow, curr_alpha_pow);
    n_poly_swap(curr_alpha_pow, next_alpha_pow);

    goto more;

cleanup:

    flint_free(starts);
    flint_free(link);
    flint_free(lift_fac);

    nmod_mat_clear(N);

    n_poly_clear(final_alpha_pow);
    n_poly_clear(curr_alpha_pow);
    n_poly_clear(prev_alpha_pow);
    n_poly_clear(next_alpha_pow);
    n_poly_clear(p1);
    n_poly_clear(alpha);
    n_bpoly_clear(monicB);
    n_tpoly_clear(tmp);

    fq_nmod_poly_clear(Beval, ctx);
    fq_nmod_poly_factor_clear(local_fac, ctx);
    fq_nmod_clear(Blc, ctx);

    bad_fq_nmod_mpoly_embed_chooser_clear(embc, ectx_mock, ctx_mock, state);

    return success;
}

