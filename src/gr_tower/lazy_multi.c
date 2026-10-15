/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Lazy fields: generators of functions of several arguments (pFq). The
    arguments are brought into one tower (merging their towers pairwise),
    and the generator is that of the tower with the same definition (the
    same function of equal arguments), or a new one. A generator whose
    arguments are all rational is hash-consed in a tower of its own, like
    the constants of the functions of one argument.
*/

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "fmpq_vec.h"
#include "gr_tower/lazy_impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/*
    Brings the n elements args into a common tower: sets *U and data[i]
    (initialized by this function in the context *mctx_out, the current
    context of *U at return) to the representation of args[i] in *U.
    Returns GR_UNABLE if no common tower was found; *U = NULL when all
    the elements are rational (data is not initialized then).
*/
static int
_gr_tower_lazy_common_n(gr_tower_flat_struct ** U, fmpz_mpoly_q_struct * data, fmpz_mpoly_ctx_struct ** mctx_out,
    const gr_tower_lazy_elem_struct * args, slong n, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F = NULL;
    fmpz_mpoly_ctx_struct * mctx = NULL;
    int * rational, * have;
    fmpq * q;
    slong i, j;
    int status = GR_SUCCESS;

    rational = flint_calloc(n, sizeof(int));
    have = flint_calloc(n, sizeof(int));
    q = _fmpq_vec_init(n);

    for (i = 0; i < n && status == GR_SUCCESS; i++)
    {
        gr_tower_lazy_elem_struct * x = _gr_tower_lazy_flat_view((gr_tower_lazy_elem_struct *) args + i);

        if (x->level == 0 && fmpz_mpoly_q_get_fmpq(q + i, &x->elem.flat.data, x->elem.flat.mctx))
        {
            rational[i] = 1;
            continue;
        }

        if (F == NULL)
        {
            F = x->F;
            _gr_tower_lazy_ref(F);
            gr_tower_flat_ensure(F);
            mctx = F->mctx;
            fmpz_mpoly_q_init(data + i, mctx);
            gr_tower_flat_convert(data + i, &x->elem.flat.data, x->elem.flat.mctx, F);
            have[i] = 1;
            continue;
        }

        if (x->F == F)
        {
            /* (the context of F may have moved on: all the data follows) */
            gr_tower_flat_ensure(F);
            if (F->mctx != mctx)
            {
                for (j = 0; j < i; j++)
                {
                    if (have[j])
                    {
                        fmpz_mpoly_q_t t;
                        fmpz_mpoly_q_init(t, F->mctx);
                        gr_tower_flat_convert(t, data + j, mctx, F);
                        fmpz_mpoly_q_clear(data + j, mctx);
                        fmpz_mpoly_q_init(data + j, F->mctx);
                        fmpz_mpoly_q_swap(data + j, t, F->mctx);
                        fmpz_mpoly_q_clear(t, F->mctx);
                    }
                }
                mctx = F->mctx;
            }
            fmpz_mpoly_q_init(data + i, mctx);
            gr_tower_flat_convert(data + i, &x->elem.flat.data, x->elem.flat.mctx, F);
            have[i] = 1;
            continue;
        }

        {
            /* merge F (through its element of the highest level, whose
               prefix contains the generators of all the others) with the
               tower of x */
            gr_tower_lazy_elem_struct acc;
            gr_tower_flat_struct * V = NULL;
            fmpz_mpoly_q_t xx, yy;
            slong jmax = -1, lmax = -1, lacc;

            gr_tower_flat_ensure(F);
            for (j = 0; j < i; j++)
            {
                if (have[j])
                {
                    slong l;
                    if (F->mctx != mctx)
                    {
                        fmpz_mpoly_q_t t;
                        fmpz_mpoly_q_init(t, F->mctx);
                        gr_tower_flat_convert(t, data + j, mctx, F);
                        fmpz_mpoly_q_clear(data + j, mctx);
                        fmpz_mpoly_q_init(data + j, F->mctx);
                        fmpz_mpoly_q_swap(data + j, t, F->mctx);
                        fmpz_mpoly_q_clear(t, F->mctx);
                    }
                    l = gr_tower_flat_level(data + j, F);
                    if (l > lmax)
                    {
                        lmax = l;
                        jmax = j;
                    }
                }
            }
            mctx = F->mctx;   /* (all the data was converted above) */

            _gr_tower_lazy_init(&acc, ctx);
            _gr_tower_lazy_fresh(&acc, F, F->T->num_gens, ctx);
            fmpz_mpoly_q_set(&acc.elem.flat.data, data + jmax, mctx);
            acc.reduced_version = 0;
            _gr_tower_lazy_shrink(&acc);
            lacc = acc.level;

            status = _gr_tower_lazy_common(&V, xx, yy, &acc, x, ctx);
            if (status != GR_SUCCESS || V == NULL)
            {
                if (V != NULL)
                {
                    fmpz_mpoly_q_clear(xx, V->mctx);
                    fmpz_mpoly_q_clear(yy, V->mctx);
                }
                _gr_tower_lazy_clear(&acc, ctx);
                status = GR_UNABLE;
                break;
            }

            /* the other elements of F, mapped to V */
            gr_tower_flat_ensure(V);
            for (j = 0; j < i && status == GR_SUCCESS; j++)
            {
                fmpz_mpoly_q_t t;
                if (!have[j])
                    continue;
                fmpz_mpoly_q_init(t, V->mctx);
                if (j == jmax)
                    fmpz_mpoly_q_set(t, xx, V->mctx);
                else if (V == F)
                    gr_tower_flat_convert(t, data + j, mctx, F);
                else
                    status = _gr_tower_lazy_map_element(t, data + j, mctx, F, lacc, V, ctx);
                fmpz_mpoly_q_clear(data + j, mctx);
                fmpz_mpoly_q_init(data + j, V->mctx);
                fmpz_mpoly_q_swap(data + j, t, V->mctx);
                fmpz_mpoly_q_clear(t, V->mctx);
            }

            fmpz_mpoly_q_init(data + i, V->mctx);
            fmpz_mpoly_q_swap(data + i, yy, V->mctx);
            have[i] = 1;

            fmpz_mpoly_q_clear(xx, V->mctx);
            fmpz_mpoly_q_clear(yy, V->mctx);
            _gr_tower_lazy_clear(&acc, ctx);

            _gr_tower_lazy_ref(V);
            _gr_tower_lazy_unref(LAZY(ctx), F);
            F = V;
            mctx = V->mctx;

            /* (the data may already be in a newer context of V) */
            gr_tower_flat_ensure(F);
            if (F->mctx != mctx)
            {
                for (j = 0; j <= i; j++)
                {
                    if (have[j])
                    {
                        fmpz_mpoly_q_t t;
                        fmpz_mpoly_q_init(t, F->mctx);
                        gr_tower_flat_convert(t, data + j, mctx, F);
                        fmpz_mpoly_q_clear(data + j, mctx);
                        fmpz_mpoly_q_init(data + j, F->mctx);
                        fmpz_mpoly_q_swap(data + j, t, F->mctx);
                        fmpz_mpoly_q_clear(t, F->mctx);
                    }
                }
                mctx = F->mctx;
            }
        }
    }

    if (status == GR_SUCCESS && F != NULL)
    {
        gr_tower_flat_ensure(F);
        for (i = 0; i < n; i++)
        {
            if (have[i] && F->mctx != mctx)
            {
                fmpz_mpoly_q_t t;
                fmpz_mpoly_q_init(t, F->mctx);
                gr_tower_flat_convert(t, data + i, mctx, F);
                fmpz_mpoly_q_clear(data + i, mctx);
                fmpz_mpoly_q_init(data + i, F->mctx);
                fmpz_mpoly_q_swap(data + i, t, F->mctx);
                fmpz_mpoly_q_clear(t, F->mctx);
            }
        }
        mctx = F->mctx;
        for (i = 0; i < n; i++)
        {
            if (rational[i])
            {
                fmpz_mpoly_q_init(data + i, mctx);
                fmpz_mpoly_q_set_fmpq(data + i, q + i, mctx);
                have[i] = 1;
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
            if (have[i])
                fmpz_mpoly_q_clear(data + i, mctx);
        if (F != NULL)
            _gr_tower_lazy_unref(LAZY(ctx), F);
        F = NULL;
    }

    *U = F;
    *mctx_out = mctx;

    /* (the caller holds the reference to F taken here, and releases it
       with _gr_tower_lazy_unref) */

    flint_free(rational);
    flint_free(have);
    _fmpq_vec_clear(q, n);
    return status;
}

/* the hash-consed generator kind(param) of rational arguments c[0], ..., c[n - 1] */
static int
_gr_tower_lazy_const_multi(gr_tower_lazy_elem_t res, int kind, slong param, const fmpq * c, slong n, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_flat_struct * F;
    gr_tower_lazy_const_trans_entry_struct * e;
    fmpz_mpoly_q_struct * u;
    slong i, d;
    int status;

    for (i = 0; i < L->num_const_trans; i++)
    {
        e = L->const_trans + i;
        if (e->kind == kind && e->param == param && e->nxs == n - 1 && fmpq_equal(&e->x, c) &&
            _fmpq_vec_equal(e->xs, c + 1, n - 1))
        {
            _gr_tower_lazy_set_gen_d(res, e->F, gr_tower_gid_order(e->F->T, e->gid), ctx);
            return GR_SUCCESS;
        }
    }

    F = _gr_tower_lazy_new_tower(ctx);
    gr_tower_flat_ensure(F);
    u = flint_malloc(sizeof(fmpz_mpoly_q_struct) * n);
    for (i = 0; i < n; i++)
    {
        fmpz_mpoly_q_init(u + i, F->mctx);
        fmpz_mpoly_q_set_fmpq(u + i, c + i, F->mctx);
    }
    status = gr_tower_adjoin_special_multi_flat(F->T, kind, param, u, n, F->mctx, NULL);
    for (i = 0; i < n; i++)
        fmpz_mpoly_q_clear(u + i, F->mctx);
    flint_free(u);

    if (status != GR_SUCCESS)
        return status;   /* (the new tower is empty: it is collected) */

    d = F->T->num_gens - 1;
    _gr_tower_lazy_new_def(GR_TOWER_GEN(F->T, d), F->T, ctx);

    if (L->num_const_trans == L->alloc_const_trans)
    {
        L->alloc_const_trans = FLINT_MAX(4, 2 * L->alloc_const_trans);
        L->const_trans = flint_realloc(L->const_trans, L->alloc_const_trans * sizeof(gr_tower_lazy_const_trans_entry_struct));
    }
    e = L->const_trans + L->num_const_trans;
    e->kind = kind;
    e->param = param;
    fmpq_init(&e->x);
    fmpq_init(&e->y);
    fmpq_set(&e->x, c);
    e->nxs = n - 1;
    e->xs = (n > 1) ? _fmpq_vec_init(n - 1) : NULL;
    for (i = 0; i < n - 1; i++)
        fmpq_set(e->xs + i, c + 1 + i);
    e->F = F;
    F->gc |= GR_TOWER_GC_PINNED;
    e->gid = GR_TOWER_GEN(F->T, d)->gid;
    e->def_id = GR_TOWER_GEN(F->T, d)->def_id;
    L->num_const_trans++;

    _gr_tower_lazy_set_gen_d(res, F, d, ctx);
    return GR_SUCCESS;
}

/*
    The special function generator kind(param) at the arguments args[0],
    ..., args[n - 1] (args[0] is the argument arg of the generator), with
    no canonicalization of the arguments (see lazy_hypgeom.c).
*/
static int
_gr_tower_lazy_special_gen_multi(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_struct * args, slong n,
    int kind, slong param, const gr_tower_lazy_elem_struct * W, ulong r, gr_ctx_t ctx)
{
    gr_tower_flat_struct * U;
    fmpz_mpoly_ctx_struct * mctx;
    fmpz_mpoly_q_struct * data;
    slong i, j, found_gid = -1, nd = n;
    int status;

    /* (W = args + n: the generator is a root of x^r - W, adjoined as
       such; W joins the arguments in the common tower) */
    if (W != NULL)
    {
        FLINT_ASSERT(W == args + n);
        nd = n + 1;
    }
    data = flint_malloc(sizeof(fmpz_mpoly_q_struct) * nd);
    status = _gr_tower_lazy_common_n(&U, data, &mctx, args, nd, ctx);
    if (status == GR_SUCCESS && U == NULL && W != NULL)
    {
        /* (not expected: an argument of a theta function is nonreal) */
        flint_free(data);
        return GR_UNABLE;
    }

    if (status != GR_SUCCESS)
    {
        flint_free(data);
        return status;
    }

    if (U == NULL)
    {
        /* rational arguments */
        fmpq * c = _fmpq_vec_init(n);
        for (i = 0; i < n; i++)
        {
            gr_tower_lazy_elem_struct * x = _gr_tower_lazy_flat_view((gr_tower_lazy_elem_struct *) args + i);
            (void) fmpz_mpoly_q_get_fmpq(c + i, &x->elem.flat.data, x->elem.flat.mctx);
        }
        status = _gr_tower_lazy_const_multi(res, kind, param, c, n, ctx);
        _fmpq_vec_clear(c, n);
        flint_free(data);
        return status;
    }

    /* the arguments in the prefix of a tower continuing with at least
       LAZY_FORK_SUFFIX generators (those of other computations, in a
       long session): the generator goes into a copy of the prefix, as
       for the exponentials (_gr_tower_lazy_trans_gen), keeping the
       definition of an equal generator of another tower */
    {
        slong lev = 0;
        gr_tower_flat_ensure(U);
        for (i = 0; i < nd; i++)
        {
            fmpz_mpoly_q_t t;
            fmpz_mpoly_q_init(t, U->mctx);
            gr_tower_flat_convert(t, data + i, mctx, U);
            lev = FLINT_MAX(lev, gr_tower_flat_level(t, U));
            fmpz_mpoly_q_clear(t, U->mctx);
        }
        if (lev > 0 && U->T->num_gens - lev >= LAZY_FORK_SUFFIX)
        {
            gr_tower_flat_struct * G = _gr_tower_lazy_new_tower(ctx);
            gr_tower_set_prefix(G->T, U->T, lev);
            gr_tower_flat_ensure(G);
            gr_tower_flat_ensure(U);
            for (i = 0; i < nd; i++)
            {
                fmpz_mpoly_q_t t, u;
                fmpz_mpoly_q_init(t, U->mctx);
                fmpz_mpoly_q_init(u, G->mctx);
                gr_tower_flat_convert(t, data + i, mctx, U);
                _gr_tower_flat_transport(u, t, U, G);
                fmpz_mpoly_q_clear(data + i, mctx);
                fmpz_mpoly_q_init(data + i, G->mctx);
                fmpz_mpoly_q_swap(data + i, u, G->mctx);
                fmpz_mpoly_q_clear(t, U->mctx);
                fmpz_mpoly_q_clear(u, G->mctx);
            }
            mctx = G->mctx;
            /* (the reference of U is passed on to G) */
            _gr_tower_lazy_ref(G);
            _gr_tower_lazy_unref(LAZY(ctx), U);
            U = G;
        }
    }

    /* a generator of U with this definition */
    for (j = 0; j < U->T->num_gens && found_gid < 0; j++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(U->T, j);
        slong gid = g->gid;
        truth_t eq = T_TRUE;

        if (!((g->kind == kind) || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == kind)) || g->def_param != param)
            continue;
        if (_gr_tower_gen_num_args(g) != n)
            continue;

        for (i = 0; i < n && eq == T_TRUE; i++)
        {
            const gr_tower_flat_elem_struct * a = _gr_tower_gen_arg_ptr(g, i);
            fmpz_mpoly_q_t d, t;
            gr_tower_flat_ensure(U);
            fmpz_mpoly_q_init(d, U->mctx);
            fmpz_mpoly_q_init(t, U->mctx);
            gr_tower_flat_convert(d, &a->data, a->mctx, U);
            gr_tower_flat_convert(t, data + i, mctx, U);
            fmpz_mpoly_q_sub(d, d, t, U->mctx);
            eq = gr_tower_flat_num_is_zero(d, U);
            fmpz_mpoly_q_clear(d, U->mctx);
            fmpz_mpoly_q_clear(t, U->mctx);
            /* (g may be stale after a zero test which restructured U; if
               it was eliminated, it is not this generator) */
            if (eq == T_TRUE)
            {
                slong d0 = gr_tower_gid_order(U->T, gid);
                if (d0 < 0)
                    eq = T_UNKNOWN;
                else
                    g = GR_TOWER_GEN(U->T, d0);
            }
        }

        if (eq == T_TRUE)
            found_gid = gid;
    }

    if (found_gid >= 0)
    {
        _gr_tower_lazy_set_gen_d(res, U, gr_tower_gid_order(U->T, found_gid), ctx);
    }
    else
    {
        fmpz_mpoly_q_struct * w;
        fmpz_mpoly_ctx_struct * wctx;
        char * name;
        ulong def_id;

        /* the same definition in another tower */
        def_id = _gr_tower_lazy_find_def(&name, kind, param, args, n, U->T, ctx);

        /* the data, in the current context of U */
        gr_tower_flat_ensure(U);
        wctx = U->mctx;
        w = flint_malloc(sizeof(fmpz_mpoly_q_struct) * nd);
        for (i = 0; i < nd; i++)
        {
            fmpz_mpoly_q_init(w + i, U->mctx);
            gr_tower_flat_convert(w + i, data + i, mctx, U);
        }

        status = gr_tower_adjoin_special_multi_flat(U->T, kind, param, w, n, U->mctx, NULL);

        if (status == GR_SUCCESS && W != NULL)
        {
            /* the modulus x^r - W (W involves only earlier generators) */
            fmpz_mpoly_q_struct * m = flint_malloc(sizeof(fmpz_mpoly_q_struct) * (r + 1));
            slong d = U->T->num_gens - 1;
            gr_tower_flat_ensure(U);
            for (i = 0; i <= (slong) r; i++)
                fmpz_mpoly_q_init(m + i, U->mctx);
            gr_tower_flat_convert(m + 0, w + n, wctx, U);
            fmpz_mpoly_q_neg(m + 0, m + 0, U->mctx);
            fmpz_mpoly_q_one(m + r, U->mctx);
            if (!_gr_tower_make_algebraic(U->T, d, m, r + 1, U->mctx, GR_TOWER_STATUS_DYNAMIC))
                status = GR_UNABLE;
            for (i = 0; i <= (slong) r; i++)
                fmpz_mpoly_q_clear(m + i, U->mctx);
            flint_free(m);
        }

        for (i = 0; i < nd; i++)
            fmpz_mpoly_q_clear(w + i, wctx);
        flint_free(w);

        if (status == GR_SUCCESS)
        {
            gr_tower_gen_struct * ng = GR_TOWER_GEN(U->T, U->T->num_gens - 1);
            slong gid = ng->gid;
            _gr_tower_lazy_set_def(ng, U->T, def_id, name, ctx);
            if (W == NULL)
                _gr_tower_search_relations(U, GR_TOWER_DEFAULT_PREC);
            _gr_tower_lazy_set_gen_d(res, U, gr_tower_gid_order(U->T, gid), ctx);
        }
        flint_free(name);
    }

    for (i = 0; i < nd; i++)
        fmpz_mpoly_q_clear(data + i, mctx);
    flint_free(data);
    _gr_tower_lazy_unref(LAZY(ctx), U);

    return status;
}

int
_gr_tower_lazy_special_gen_multi_locked(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_struct * args, slong n,
    int kind, slong param, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * t;
    slong i;
    int status;

    _gr_tower_lazy_lock(ctx);

    /* (owned copies: shallow elements, and res aliasing an argument) */
    t = flint_malloc(sizeof(gr_tower_lazy_elem_struct) * n);
    status = GR_SUCCESS;
    for (i = 0; i < n; i++)
    {
        _gr_tower_lazy_init(t + i, ctx);
        status |= _gr_tower_lazy_set(t + i, args + i, ctx);
    }

    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_special_gen_multi(res, t, n, kind, param, NULL, 0, ctx);

    for (i = 0; i < n; i++)
        _gr_tower_lazy_clear(t + i, ctx);
    flint_free(t);

    _gr_tower_lazy_unlock(ctx);
    return status;
}

/*
    The generator kind(param) at args[0], ..., args[n - 1] which is
    known to be a root of x^r - W (the root selected by the numerical
    value of the function): adjoined as an algebraic generator with this
    modulus and the definition of the function, so that equal arguments
    give the same generator.
*/
int
_gr_tower_lazy_special_alg_gen_multi_locked(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_struct * args, slong n,
    int kind, slong param, const gr_tower_lazy_elem_t W, ulong r, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * t;
    slong i;
    int status;

    _gr_tower_lazy_lock(ctx);

    t = flint_malloc(sizeof(gr_tower_lazy_elem_struct) * (n + 1));
    status = GR_SUCCESS;
    for (i = 0; i <= n; i++)
    {
        _gr_tower_lazy_init(t + i, ctx);
        status |= _gr_tower_lazy_set(t + i, (i < n) ? args + i : W, ctx);
    }

    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_special_gen_multi(res, t, n, kind, param, t + n, r, ctx);

    for (i = 0; i <= n; i++)
        _gr_tower_lazy_clear(t + i, ctx);
    flint_free(t);

    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
_gr_tower_lazy_special_multi(gr_tower_lazy_elem_t res, int kind, slong param, const gr_tower_lazy_elem_struct * args, slong nargs, gr_ctx_t ctx)
{
    int status;

    if (nargs == 1)
        return _gr_tower_lazy_special(res, kind, param, args, ctx);

    _gr_tower_lazy_lock(ctx);
    if (kind == GR_TOWER_HYPGEOM)
        status = _gr_tower_lazy_hypgeom_args(res, param, args, nargs, ctx);
    else if (kind == GR_TOWER_JACOBI_THETA && nargs == 2)
        status = _gr_tower_lazy_jacobi_theta_j(res, param, args, args + 1, ctx);
    else if (kind == GR_TOWER_HURWITZ_ZETA && nargs == 2)
        status = _gr_tower_lazy_hurwitz_general(res, args, args + 1, ctx);
    else
        status = GR_UNABLE;
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* anchors: values linked to the generators of the context               */
/* -------------------------------------------------------------------- */

/* the nesting of the linking steps (Landen's transformation of K and E,
   the quadratic transformations of 2F1), changed by delta */
int
_gr_tower_lazy_hyp_anchored(gr_ctx_t ctx, int delta)
{
    LAZY(ctx)->hyp_anchored += delta;
    return LAZY(ctx)->hyp_anchored;
}

/*
    Whether a generator K(x), E(x) (or, unless elliptic_only, 2F1(...; x))
    of the context has x numerically equal to t or to t/(t - 1) (the
    canonical arguments are those of Pfaff's transformation, or of the
    imaginary-modulus transformation of K and E). A test for the
    candidates of the linking steps: the values are then computed
    exactly, and lead to the generator when it is the same.
*/
int
_gr_tower_lazy_hyp_anchor_present(const acb_t t, int elliptic_only, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_ptr e;
    acb_t u, x;
    slong i, d;
    int found = 0;

    if (!acb_is_finite(t))
        return 0;

    GR_TMP_INIT(e, ctx);
    acb_init(u);
    acb_init(x);
    acb_sub_ui(u, t, 1, 128);
    acb_div(u, t, u, 128);

    for (i = 0; i < L->num_towers && !found; i++)
    {
        gr_tower_flat_struct * F = L->towers[i];
        gr_tower_struct * T = F->T;

        for (d = 0; d < T->num_gens && !found; d++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
            int kind = (g->kind == GR_TOWER_ALGEBRAIC) ? g->def_kind : g->kind;

            if (!(kind == GR_TOWER_ELLIPTIC_K || kind == GR_TOWER_ELLIPTIC_E ||
                  (!elliptic_only && kind == GR_TOWER_HYPGEOM && g->def_param == GR_TOWER_HYPGEOM_PARAM(2, 1))))
                continue;
            if (g->arg.mctx == NULL)
                continue;
            _gr_tower_lazy_set_flat(e, F, &g->arg.data, g->arg.mctx, ctx);
            if (gr_tower_lazy_get_acb(x, e, 64, ctx) == GR_SUCCESS &&
                (acb_overlaps(x, t) || acb_overlaps(x, u)))
                found = 1;
        }
    }

    acb_clear(u);
    acb_clear(x);
    GR_TMP_CLEAR(e, ctx);
    return found;
}

POP_OPTIONS
