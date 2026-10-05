/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: common towers (merging, primitive elements, aliases). */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

/* -------------------------------------------------------------------- */
/* common towers                                                         */
/* -------------------------------------------------------------------- */

static void
_add_alias(gr_tower_flat_struct * F, ulong def_id, const fmpz_mpoly_q_t data, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_lazy_alias_struct * a;

    if (L->num_aliases == L->alloc_aliases)
    {
        L->alloc_aliases = FLINT_MAX(4, 2 * L->alloc_aliases);
        L->aliases = flint_realloc(L->aliases, L->alloc_aliases * sizeof(gr_tower_lazy_alias_struct));
    }

    a = L->aliases + L->num_aliases;
    a->F = F;
    a->def_id = def_id;
    a->mctx = F->mctx;
    fmpz_mpoly_q_init(&a->data, a->mctx);
    fmpz_mpoly_q_set(&a->data, data, a->mctx);
    L->num_aliases++;
}

/*
    Writes to img (a polynomial in F's current context) the image in F
    of the generator with definition order d of S: a generator with the
    same definition, an alias, or a power of a root of unity (or of an
    integer) generator. Returns 0 if absent. Alias data is converted to
    the current context if needed.
*/
static int
_def_image(fmpz_mpoly_q_t img, gr_tower_flat_struct * F, gr_tower_flat_struct * S, slong d, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    const gr_tower_gen_struct * g = GR_TOWER_GEN(S->T, d);
    ulong def_id = g->def_id;
    slong j, i;

    gr_tower_flat_ensure(F);

    j = gr_tower_find_def(F->T, def_id);
    if (j > 0)
    {
        fmpz_mpoly_q_gen(img, GR_TOWER_FLAT_VAR(F, j), F->mctx);
        return 1;
    }

    j = gr_tower_find_trans_def(F->T, def_id);
    if (j > 0)
    {
        fmpz_mpoly_q_gen(img, GR_TOWER_FLAT_TVAR(F, j), F->mctx);
        return 1;
    }

    for (i = 0; i < L->num_aliases; i++)
    {
        if (L->aliases[i].F == F && L->aliases[i].def_id == def_id)
        {
            gr_tower_lazy_alias_struct * a = L->aliases + i;
            if (a->mctx != F->mctx)
            {
                fmpz_mpoly_q_t t;
                fmpz_mpoly_q_init(t, F->mctx);
                gr_tower_flat_convert(t, &a->data, a->mctx, F);
                fmpz_mpoly_q_clear(&a->data, a->mctx);
                a->mctx = F->mctx;
                fmpz_mpoly_q_init(&a->data, a->mctx);
                fmpz_mpoly_q_swap(&a->data, t, a->mctx);
                fmpz_mpoly_q_clear(t, a->mctx);
            }
            fmpz_mpoly_q_set(img, &a->data, F->mctx);
            return 1;
        }
    }

    /* a root of unity or of an integer, as a product of powers of the
       roots of prime power orders present in F */
    return _gr_tower_structured_image(img, F->T, g, !(LAZY(ctx)->merge_flags & GR_TOWER_MERGE_REAL_FIRST));
}

/* Whether F has an image of the generator with definition order d of S. */
static int
_tower_has_def(gr_tower_flat_struct * F, gr_tower_flat_struct * S, slong d, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    const gr_tower_gen_struct * g = GR_TOWER_GEN(S->T, d);
    ulong def_id = g->def_id;
    slong i;

    if (gr_tower_find_def(F->T, def_id) > 0 || gr_tower_find_trans_def(F->T, def_id) > 0)
        return 1;

    for (i = 0; i < L->num_aliases; i++)
        if (L->aliases[i].F == F && L->aliases[i].def_id == def_id)
            return 1;

    if (g->kind == GR_TOWER_ALGEBRAIC && (g->def_kind == GR_TOWER_ROOT_OF_UNITY || g->def_kind == GR_TOWER_ROOT))
        return _gr_tower_structured_present(F->T, g, !(L->merge_flags & GR_TOWER_MERGE_REAL_FIRST));

    return 0;
}

/* Definition id of the generator with definition order d in T. */
static ulong
_gen_def_id(const gr_tower_struct * T, slong d)
{
    return GR_TOWER_GEN(T, d)->def_id;
}

static int
_tower_contains_prefix(gr_tower_flat_struct * F, gr_tower_flat_struct * S, slong p, gr_ctx_t ctx)
{
    slong d;

    if (F == S)
        return 1;

    for (d = 0; d < p; d++)
        if (!_tower_has_def(F, S, d, ctx))
            return 0;

    return 1;
}

/*
    Maps the flat element x (in the context x_mctx of S, using only the
    first p generators of S in definition order) to the flat
    representation in U, substituting the images of the definitions.
    res must be initialized in U->mctx.
*/
int
_gr_tower_lazy_map_element(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, gr_tower_flat_struct * S, slong p, gr_tower_flat_struct * U, gr_ctx_t ctx)
{
    slong nvars, d;
    fmpz_mpoly_q_struct ** imgs;
    fmpz_mpoly_q_t xs;
    int status = GR_SUCCESS;

    gr_tower_flat_ensure(U);
    gr_tower_flat_ensure(S);

    fmpz_mpoly_q_init(xs, S->mctx);
    gr_tower_flat_convert(xs, x, x_mctx, S);

    nvars = S->cap;
    imgs = flint_calloc(nvars, sizeof(fmpz_mpoly_q_struct *));

    for (d = 0; d < p && status == GR_SUCCESS; d++)
    {
        slong var = GR_TOWER_FLAT_VAR_D(S, d);

        imgs[var] = flint_malloc(sizeof(fmpz_mpoly_q_struct));
        fmpz_mpoly_q_init(imgs[var], U->mctx);
        if (!_def_image(imgs[var], U, S, d, ctx))
        {
            status = GR_DOMAIN;
        }
    }

    if (status == GR_SUCCESS)
        status = gr_tower_flat_compose(res, xs, S->mctx, imgs, U);

    for (d = 0; d < nvars; d++)
    {
        if (imgs[d] != NULL)
        {
            fmpz_mpoly_q_clear(imgs[d], U->mctx);
            flint_free(imgs[d]);
        }
    }
    flint_free(imgs);
    fmpz_mpoly_q_clear(xs, S->mctx);

    return status;
}

/* Whether the first p generators of S include a conjugate (same origin
   polynomial) of a generator of F: absorbing them extends a splitting
   tower, which is done in place whatever the size of the towers. */
static int
_prefix_extends_splitting(gr_tower_flat_struct * F, gr_tower_flat_struct * S, slong p)
{
    slong d, e;
    for (d = 0; d < p; d++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(S->T, d);
        if (g->kind != GR_TOWER_ALGEBRAIC || g->origin == NULL)
            continue;
        for (e = 0; e < F->T->num_gens; e++)
            if (GR_TOWER_GEN(F->T, e)->origin != NULL && fmpz_poly_equal(GR_TOWER_GEN(F->T, e)->origin, g->origin))
                return 1;
    }
    return 0;
}

/* Whether all algebraic generators among the first p of S are present
   in F (by definition id). */
static int
_prefix_alg_covered(gr_tower_flat_struct * F, gr_tower_flat_struct * S, slong p, gr_ctx_t ctx)
{
    slong d;
    for (d = 0; d < p; d++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(S->T, d);

        if (g->kind != GR_TOWER_ALGEBRAIC || _tower_has_def(F, S, d, ctx))
            continue;

        /* a conjugate of a generator of F (same origin polynomial)
           extends the splitting tower of that polynomial in place */
        if (g->origin != NULL)
        {
            slong e;
            int found = 0;
            for (e = 0; e < F->T->num_gens && !found; e++)
                if (GR_TOWER_GEN(F->T, e)->origin != NULL && fmpz_poly_equal(GR_TOWER_GEN(F->T, e)->origin, g->origin))
                    found = 1;
            if (found)
                continue;
        }

        return 0;
    }
    return 1;
}

/* Generators created inside the towers (by the relation search: a
   primitive exponential generator inserted before its powers, say) get
   definition ids, so that they are found again in copies of the tower. */
void
_gr_tower_lazy_assign_def_ids(gr_tower_flat_struct * F, gr_ctx_t ctx)
{
    slong k;
    for (k = 0; k < F->T->num_gens; k++)
        if (GR_TOWER_GEN(F->T, k)->def_id == 0)
            _gr_tower_lazy_new_def(GR_TOWER_GEN(F->T, k), F->T, ctx);
}

/*
    Whether the transcendental generators of U are known to be linearly
    independent over Q together with pi i without a search: pi, and the
    constant logarithms log(p) of rational primes p and log(a + b i) of
    Gaussian primes (a, b > 0, a^2 + b^2 prime). A relation
    sum c_j log(q_j) = c pi i would give prod q_j^{c_j} = a root of
    unity, impossible for distinct primes up to units by unique
    factorization in Z and Z[i] (the generators of one value are
    hash-consed, and the conjugate a - b i is never a generator: see
    _gr_tower_lazy_log of Gaussian rationals).
*/
static int
_gr_tower_lazy_trans_gens_independent(gr_tower_flat_struct * U, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_struct * T = U->T;
    slong d, i;

    for (d = 0; d < T->num_gens; d++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
        const gr_tower_lazy_const_trans_entry_struct * e = NULL;

        if (g->kind == GR_TOWER_ALGEBRAIC)
            continue;
        if (g->kind == GR_TOWER_PI)
            continue;
        if (g->kind != GR_TOWER_LOG || g->def_id == 0)
            return 0;

        for (i = 0; i < L->num_const_trans; i++)
            if (L->const_trans[i].def_id == g->def_id)
                e = L->const_trans + i;
        if (e == NULL || e->kind != GR_TOWER_LOG)
            return 0;
        if (!fmpz_is_one(fmpq_denref(&e->x)) || !fmpz_is_one(fmpq_denref(&e->y)))
            return 0;
        if (fmpz_sgn(fmpq_numref(&e->x)) <= 0 || fmpz_sgn(fmpq_numref(&e->y)) < 0)
            return 0;
        if (!fmpz_abs_fits_ui(fmpq_numref(&e->x)) || !fmpz_abs_fits_ui(fmpq_numref(&e->y)))
            return 0;

        if (fmpz_is_zero(fmpq_numref(&e->y)))
        {
            if (!n_is_prime(fmpz_get_ui(fmpq_numref(&e->x))))
                return 0;
        }
        else
        {
            ulong a = fmpz_get_ui(fmpq_numref(&e->x)), b = fmpz_get_ui(fmpq_numref(&e->y)), N;
            if (a >= (UWORD(1) << (FLINT_BITS / 2 - 1)) || b >= (UWORD(1) << (FLINT_BITS / 2 - 1)))
                return 0;
            N = a * a + b * b;
            if (!n_is_prime(N))
                return 0;
        }
    }

    return 1;
}


/*
    Finds a tower containing the first lx definitions of Fx and the first
    ly definitions of Fy, creating one if necessary.
*/
gr_tower_flat_struct *
_gr_tower_lazy_common_tower(gr_tower_flat_struct * Fx, slong lx, gr_tower_flat_struct * Fy, slong ly, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_flat_struct * U = NULL;
    slong i;

    _gr_tower_lazy_assign_def_ids(Fx, ctx);
    _gr_tower_lazy_assign_def_ids(Fy, ctx);

    /* one of the towers already contains the other prefix */
    if (_tower_contains_prefix(Fx, Fy, ly, ctx))
        return Fx;
    if (_tower_contains_prefix(Fy, Fx, lx, ctx))
        return Fy;

    for (i = 0; i < L->num_towers; i++)
    {
        gr_tower_flat_struct * V = L->towers[i];

        if (U != NULL && (gr_tower_degree(V->T) > gr_tower_degree(U->T) ||
                (gr_tower_degree(V->T) == gr_tower_degree(U->T) && V->T->num_gens >= U->T->num_gens)))
            continue;

        if (_tower_contains_prefix(V, Fx, lx, ctx) && _tower_contains_prefix(V, Fy, ly, ctx))
            U = V;
    }

    if (U != NULL)
        return U;

    /*
        Small towers are combined into a new tower (the prefix of Fx
        extended by the prefix of Fy), so that small fields stay small.
        Beyond a few generators, the prefix of the smaller tower is
        absorbed into the larger tower in place: flat elements survive
        extensions, and this avoids a combinatorial explosion of towers
        (a sum of many logarithms would otherwise create a tower for
        each partial sum).
    */
    if ((lx + ly > LAZY(ctx)->options[GR_TOWER_OPT_GROW_THRESHOLD] || _prefix_extends_splitting(Fx, Fy, ly) || _prefix_extends_splitting(Fy, Fx, lx)) &&
        Fx != L->trivial && Fy != L->trivial &&
        (_prefix_alg_covered(Fx, Fy, ly, ctx) || _prefix_alg_covered(Fy, Fx, lx, ctx)))
    {
        /* the algebraic generators of the absorbed prefix are already in
           the growing tower, so that its degree does not change */
        if (!_prefix_alg_covered(Fx, Fy, ly, ctx) ||
            (_prefix_alg_covered(Fy, Fx, lx, ctx) && Fy->T->num_gens > Fx->T->num_gens))
        {
            gr_tower_flat_struct * t = Fx; Fx = Fy; Fy = t;
            i = lx; lx = ly; ly = i;
        }
        U = Fx;
    }
    else
    {
        U = _gr_tower_lazy_new_tower(ctx);
        gr_tower_set_prefix(U->T, Fx->T, lx);
    }

    {
        gr_tower_map_t my;
        slong k, ntrans0 = U->T->num_trans;

        gr_tower_map_init(my, Fy->T, U->T);
        if (gr_tower_absorb_prefix(U->T, my, Fy->T, ly, L->merge_flags) != GR_SUCCESS)
        {
            gr_tower_map_clear(my);
            return NULL;
        }

        gr_tower_map_sync(my);

        /* generators created by the absorption itself (a root of unity of
           a larger order, say) are new definitions */
        for (k = 0; k < U->T->num_gens; k++)
            if (GR_TOWER_GEN(U->T, k)->def_id == 0)
                _gr_tower_lazy_new_def(GR_TOWER_GEN(U->T, k), U->T, ctx);

        /* record definitions of Fy which were expressed rather than adjoined */
        for (k = 0; k < ly; k++)
        {
            ulong id = _gen_def_id(Fy->T, k);
            if (gr_tower_find_def(U->T, id) == 0 && gr_tower_find_trans_def(U->T, id) == 0)
            {
                _add_alias(U, id, gr_tower_map_image(my, k), ctx);
            }
        }

        gr_tower_map_clear(my);

        /* transcendental generators from both sides may be related (a
           tangent of a copied prefix and those of the growing tower):
           a cheap search now, as at adjunction, rather than in the next
           zero test, which would first convert its operands to the
           nested representation over unreduced generators (not for
           logarithms of distinct primes, which are independent) */
        if (ntrans0 > 0 && U->T->num_trans > ntrans0 && !_gr_tower_lazy_trans_gens_independent(U, ctx))
            _gr_tower_search_relations(U, GR_TOWER_DEFAULT_PREC);
    }

    return U;
}

/*
    Brings x and y into a common tower: sets *U and writes copies xx, yy
    in U's current context (the caller clears them with U->mctx).
*/
/*
    Primitive elements (GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT): a tower U
    over QQ of several proven algebraic steps (no transcendental
    generators, no roots of unity or tangents of rational multiples of pi,
    whose structure the other operations use), of degree D within the
    limit, is replaced by the fresh tower V = QQ(theta) for theta = l sum
    c_i a_i over the steps a_i of degree > 1, with small integers c_i such
    that 1, theta, ..., theta^(D-1) are linearly independent (then the
    minimal polynomial of theta has degree D), scaled by an integer l
    making it an algebraic integer. V is the tower of the algebraic number
    theta (hash-consed); the definitions of U are recorded as aliases in V
    (polynomials in theta), so that the merges which follow find V
    containing both operands. Elements of V have the dense form.
*/

/* the coordinates of the flat element x (reduced, integer denominator)
   on the monomials prod a_k^(e_k), index sum e_k C_k */
int
_gr_tower_lazy_prim_coords(fmpq * v, const fmpz_mpoly_q_t x, const slong * var, const slong * C, slong m, slong D, const fmpz_mpoly_ctx_t mctx)
{
    slong nvars = mctx->minfo->nvars, t, k;
    ulong * e;
    fmpz_t den;
    int ok = 1;

    if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), mctx))
        return 0;

    e = flint_malloc(sizeof(ulong) * nvars);
    fmpz_init(den);
    fmpz_mpoly_get_fmpz(den, fmpz_mpoly_q_denref(x), mctx);
    for (k = 0; k < D; k++)
        fmpq_zero(v + k);

    for (t = 0; t < fmpz_mpoly_q_numref(x)->length && ok; t++)
    {
        slong idx = 0, j;
        fmpz_mpoly_get_term_exp_ui(e, fmpz_mpoly_q_numref(x), t, mctx);
        for (k = 0; k < m; k++)
        {
            idx += (slong) e[var[k]] * C[k];
            e[var[k]] = 0;
        }
        for (j = 0; j < nvars; j++)
            if (e[j] != 0)
                ok = 0;
        if (ok && idx < D)
        {
            fmpz_mpoly_get_term_coeff_fmpz(fmpq_numref(v + idx), fmpz_mpoly_q_numref(x), t, mctx);
            fmpz_set(fmpq_denref(v + idx), den);
            fmpq_canonicalise(v + idx);
        }
        else
            ok = 0;
    }

    fmpz_clear(den);
    flint_free(e);
    return ok;
}

/* whether V has the definition def_id (a generator or an alias) */
static int
_tower_has_def_id(gr_tower_flat_struct * V, ulong def_id, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    slong i;
    if (gr_tower_find_def(V->T, def_id) > 0 || gr_tower_find_trans_def(V->T, def_id) > 0)
        return 1;
    for (i = 0; i < L->num_aliases; i++)
        if (L->aliases[i].F == V && L->aliases[i].def_id == def_id)
            return 1;
    return 0;
}

static gr_tower_flat_struct *
_gr_tower_lazy_primitive_tower(gr_tower_flat_struct * U, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_struct * T = U->T;
    gr_tower_flat_struct * V = NULL;
    slong limit = L->options[GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT];
    slong m = 0, D = 1, k, i, j, attempt;
    slong * var, * deg, * C;
    fmpq_mat_t X;
    fmpz_mpoly_q_t theta;
    int found = 0;

    if (limit <= 0 || !GR_TOWER_BASE_IS_CONSTS(T) || T->length < 2 || !GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL))
        return NULL;

    /* (sought once per state of the tower) */
    if (U->primitive_tried[0] == T->version + 1 && U->primitive_tried[1] == (ulong) T->length)
        return NULL;
    U->primitive_tried[0] = T->version + 1;
    U->primitive_tried[1] = T->length;

    gr_tower_flat_ensure(U);
    for (k = 1; k <= T->length; k++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_STEP(T, k - 1);
        slong d = gr_tower_step_degree(T, k);
        if (d > 1)
        {
            /* (linearly disjoint simple extensions of QQ: proven steps
               with monic integral univariate moduli, not roots of unity or
               tangents) */
            if (g->status != GR_TOWER_STATUS_PROVEN || g->def_kind == GR_TOWER_ROOT_OF_UNITY || g->def_kind == GR_TOWER_TAN_PI ||
                U->ideal_univar == NULL || U->ideal_univar[k - 1] == NULL)
                return NULL;
            m++;
            D *= d;
            if (D > limit)
                return NULL;
        }
    }
    if (m < 2)
        return NULL;
    for (k = 0; k < T->num_gens; k++)
        if (GR_TOWER_GEN(T, k)->def_id == 0)
            return NULL;

    gr_tower_flat_ensure(U);
    var = flint_malloc(sizeof(slong) * 3 * (m + 1));
    deg = var + m + 1;
    C = deg + m + 1;
    for (k = 1, i = 0; k <= T->length; k++)
    {
        slong d = gr_tower_step_degree(T, k);
        if (d > 1)
        {
            var[i] = GR_TOWER_FLAT_VAR(U, k);
            deg[i] = d;
            i++;
        }
    }
    C[0] = 1;
    for (i = 0; i < m; i++)
        C[i + 1] = C[i] * deg[i];

    fmpq_mat_init(X, D, m + 1);
    fmpz_mpoly_q_init(theta, U->mctx);

    for (attempt = 0; attempt < 8 && !found; attempt++)
    {
        fmpz_mpoly_q_t pw, t;
        fmpq_mat_t P, R;
        fmpq * v;
        int ok = 1;

        fmpz_mpoly_q_init(pw, U->mctx);
        fmpz_mpoly_q_init(t, U->mctx);
        fmpq_mat_init(P, D, D);
        fmpq_mat_init(R, D, m + 1);
        v = _fmpq_vec_init(D);

        /* theta = sum c_i a_i, c_i = 1 + (attempt i) mod (2 attempt + 3) */
        fmpz_mpoly_q_zero(theta, U->mctx);
        for (i = 0; i < m; i++)
        {
            fmpz_mpoly_q_gen(t, var[i], U->mctx);
            fmpz_mpoly_q_mul_si(t, t, 1 + (attempt * i) % (2 * attempt + 3), U->mctx);
            fmpz_mpoly_q_add(theta, theta, t, U->mctx);
        }

        /* the coordinates of theta^j (j < D) as the columns of P; theta^D
           and the a_i as the right-hand sides */
        fmpz_mpoly_q_one(pw, U->mctx);
        for (j = 0; j <= D && ok; j++)
        {
            ok = _gr_tower_lazy_prim_coords(v, pw, var, C, m, D, U->mctx);
            for (i = 0; i < D && ok; i++)
                fmpq_set((j < D) ? fmpq_mat_entry(P, i, j) : fmpq_mat_entry(R, i, 0), v + i);
            if (j < D && ok)
            {
                fmpz_mpoly_q_mul(pw, pw, theta, U->mctx);
                ok = (gr_tower_flat_reduce(pw, U) == GR_SUCCESS);
            }
        }
        for (i = 0; i < m && ok; i++)
            fmpq_set_si(fmpq_mat_entry(R, C[i], i + 1), 1, 1);

        if (ok && fmpq_mat_solve(X, P, R))
            found = 1;

        fmpz_mpoly_q_clear(pw, U->mctx);
        fmpz_mpoly_q_clear(t, U->mctx);
        fmpq_mat_clear(P);
        fmpq_mat_clear(R);
        _fmpq_vec_clear(v, D);
    }

    if (found)
    {
        /* theta^D = sum X[j][0] theta^j; theta' = l theta has the integral
           minimal polynomial x^D - sum X[j][0] l^(D-j) x^j */
        fmpz_t l, c;
        fmpq_t q;
        fmpz_poly_t mp;
        qqbar_t th;
        acb_t z;
        slong prec, gid = -1;
        int valid = 0;

        fmpz_init(l);
        fmpz_init(c);
        fmpq_init(q);
        fmpz_poly_init(mp);
        qqbar_init(th);
        acb_init(z);

        /* the least l making the coefficients integral: for each prime p
           of the denominators, v_p(l) = max_j ceil(v_p(den_j) / (D - j))
           (the lcm of the denominators when they do not factor quickly) */
        fmpz_one(l);
        for (j = 0; j < D; j++)
            fmpz_lcm(l, l, fmpq_denref(fmpq_mat_entry(X, j, 0)));
        if (!fmpz_is_one(l) && fmpz_bits(l) <= 128)
        {
            fmpz_factor_t fac;
            fmpz_t pe, r;
            slong pi;
            fmpz_factor_init(fac);
            fmpz_init(pe);
            fmpz_init(r);
            fmpz_factor(fac, l);
            fmpz_one(l);
            for (pi = 0; pi < fac->num; pi++)
            {
                slong e = 0;
                for (j = 0; j < D; j++)
                {
                    slong v = fmpz_remove(r, fmpq_denref(fmpq_mat_entry(X, j, 0)), fac->p + pi);
                    e = FLINT_MAX(e, (v + (D - j) - 1) / (D - j));
                }
                fmpz_pow_ui(pe, fac->p + pi, e);
                fmpz_mul(l, l, pe);
            }
            fmpz_factor_clear(fac);
            fmpz_clear(pe);
            fmpz_clear(r);
        }
        fmpz_poly_set_coeff_si(mp, D, 1);
        for (j = 0; j < D; j++)
        {
            fmpz_pow_ui(c, l, D - j);
            fmpq_mul_fmpz(q, fmpq_mat_entry(X, j, 0), c);
            fmpz_neg(c, fmpq_numref(q));        /* (an integer) */
            fmpz_poly_set_coeff_fmpz(mp, j, c);
        }

        /* an enclosure of theta' isolating it among the roots */
        for (prec = 2 * GR_TOWER_DEFAULT_PREC; prec <= GR_TOWER_OPTION(T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT) && !valid; prec *= 2)
        {
            if (gr_tower_flat_get_acb(z, theta, prec, U) == GR_SUCCESS)
            {
                acb_mul_fmpz(z, z, l, prec);
                valid = qqbar_set_fmpz_poly_root(th, mp, z, 2 * prec);
            }
        }

        if (valid)
        {
            if (_gr_tower_lazy_qqbar_tower(&V, &gid, th, GR_TOWER_ALGEBRAIC, 0, NULL, ctx) != GR_SUCCESS)
                V = NULL;
        }

        if (V != NULL)
        {
            /* the images of the steps (by variable of U), then the aliases */
            slong nv = U->mctx->minfo->nvars, vt;
            fmpz_mpoly_q_struct ** imgs = flint_calloc(nv, sizeof(fmpz_mpoly_q_struct *));
            fmpz_mpoly_q_t g, img;

            gr_tower_flat_ensure(V);
            vt = GR_TOWER_FLAT_VAR_D(V, gr_tower_gid_order(V->T, gid));
            fmpz_mpoly_q_init(g, V->mctx);
            fmpz_mpoly_q_init(img, V->mctx);

            for (i = 0; i < m; i++)
            {
                imgs[var[i]] = flint_malloc(sizeof(fmpz_mpoly_q_struct));
                fmpz_mpoly_q_init(imgs[var[i]], V->mctx);
                /* a_i = sum X[j][i+1] l^(-j) theta'^j */
                for (j = D - 1; j >= 0; j--)
                {
                    fmpz_mpoly_q_gen(g, vt, V->mctx);
                    fmpz_mpoly_q_mul(imgs[var[i]], imgs[var[i]], g, V->mctx);
                    fmpz_pow_ui(c, l, j);
                    fmpq_div_fmpz(q, fmpq_mat_entry(X, j, i + 1), c);
                    fmpz_mpoly_q_add_fmpq(imgs[var[i]], imgs[var[i]], q, V->mctx);
                }
                GR_MUST_SUCCEED(gr_tower_flat_reduce(imgs[var[i]], V));
            }

            /* every generator of U (the steps of degree one through their
               values), and the aliases recorded in U */
            for (k = 0; k < T->num_gens; k++)
            {
                const gr_tower_gen_struct * h = GR_TOWER_GEN(T, k);
                if (_tower_has_def_id(V, h->def_id, ctx))
                    continue;
                {
                    fmpz_mpoly_q_t u;
                    fmpz_mpoly_q_init(u, U->mctx);
                    fmpz_mpoly_q_gen(u, GR_TOWER_FLAT_VAR_D(U, k), U->mctx);
                    if (gr_tower_flat_reduce(u, U) == GR_SUCCESS &&
                        gr_tower_flat_compose(img, u, U->mctx, imgs, V) == GR_SUCCESS &&
                        gr_tower_flat_reduce(img, V) == GR_SUCCESS)
                        _add_alias(V, h->def_id, img, ctx);
                    fmpz_mpoly_q_clear(u, U->mctx);
                }
            }
            for (k = 0; k < L->num_aliases; k++)
            {
                gr_tower_lazy_alias_struct * a = L->aliases + k;
                if (a->F != U || _tower_has_def_id(V, a->def_id, ctx))
                    continue;
                {
                    fmpz_mpoly_q_t u;
                    ulong id = a->def_id;
                    fmpz_mpoly_q_init(u, U->mctx);
                    gr_tower_flat_convert(u, &a->data, a->mctx, U);
                    if (gr_tower_flat_reduce(u, U) == GR_SUCCESS &&
                        gr_tower_flat_compose(img, u, U->mctx, imgs, V) == GR_SUCCESS &&
                        gr_tower_flat_reduce(img, V) == GR_SUCCESS)
                        _add_alias(V, id, img, ctx);   /* (may move L->aliases) */
                    fmpz_mpoly_q_clear(u, U->mctx);
                }
            }

            for (i = 0; i < nv; i++)
            {
                if (imgs[i] != NULL)
                {
                    fmpz_mpoly_q_clear(imgs[i], V->mctx);
                    flint_free(imgs[i]);
                }
            }
            flint_free(imgs);
            fmpz_mpoly_q_clear(g, V->mctx);
            fmpz_mpoly_q_clear(img, V->mctx);

            /* (every definition of U must be present) */
            for (k = 0; k < T->num_gens && V != NULL; k++)
                if (!_tower_has_def_id(V, GR_TOWER_GEN(T, k)->def_id, ctx))
                    V = NULL;
        }

        fmpz_clear(l);
        fmpz_clear(c);
        fmpq_clear(q);
        fmpz_poly_clear(mp);
        qqbar_clear(th);
        acb_clear(z);
    }

    fmpz_mpoly_q_clear(theta, U->mctx);
    fmpq_mat_clear(X);
    flint_free(var);
    return V;
}

int
_gr_tower_lazy_common(gr_tower_flat_struct ** U, fmpz_mpoly_q_t xx, fmpz_mpoly_q_t yy,
    gr_tower_lazy_elem_t x, gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    x = _gr_tower_lazy_flat_view(x);
    y = _gr_tower_lazy_flat_view(y);

    if (x->F == y->F)
    {
        *U = x->F;
        fmpz_mpoly_q_init(xx, (*U)->mctx);
        fmpz_mpoly_q_init(yy, (*U)->mctx);
        fmpz_mpoly_q_set(xx, &x->elem.flat.data, (*U)->mctx);
        fmpz_mpoly_q_set(yy, &y->elem.flat.data, (*U)->mctx);
    }
    else if (x->level == 0 || y->level == 0)
    {
        /* a rational operand: set directly in the tower of the other */
        gr_tower_lazy_elem_struct * r = (x->level == 0) ? x : y;
        gr_tower_lazy_elem_struct * o = (x->level == 0) ? y : x;
        fmpz_t cn, cd;

        *U = o->F;
        fmpz_mpoly_q_init(xx, (*U)->mctx);
        fmpz_mpoly_q_init(yy, (*U)->mctx);
        fmpz_init(cn);
        fmpz_init(cd);
        fmpz_mpoly_get_fmpz(cn, fmpz_mpoly_q_numref(&r->elem.flat.data), r->elem.flat.mctx);
        fmpz_mpoly_get_fmpz(cd, fmpz_mpoly_q_denref(&r->elem.flat.data), r->elem.flat.mctx);
        if (r == x)
        {
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_numref(xx), cn, (*U)->mctx);
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(xx), cd, (*U)->mctx);
            fmpz_mpoly_q_set(yy, &y->elem.flat.data, (*U)->mctx);
        }
        else
        {
            fmpz_mpoly_q_set(xx, &x->elem.flat.data, (*U)->mctx);
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_numref(yy), cn, (*U)->mctx);
            fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(yy), cd, (*U)->mctx);
        }
        fmpz_clear(cn);
        fmpz_clear(cd);
    }
    else
    {
        *U = _gr_tower_lazy_common_tower(x->F, x->level, y->F, y->level, ctx);
        if (*U == NULL)
            return GR_UNABLE;

        fmpz_mpoly_q_init(xx, (*U)->mctx);
        fmpz_mpoly_q_init(yy, (*U)->mctx);
        if (x->F == *U)
        {
            x = _gr_tower_lazy_flat_view(x);
            fmpz_mpoly_q_set(xx, &x->elem.flat.data, (*U)->mctx);
        }
        else
            status |= _gr_tower_lazy_map_element(xx, &x->elem.flat.data, x->elem.flat.mctx, x->F, x->level, *U, ctx);
        if (y->F == *U)
        {
            y = _gr_tower_lazy_flat_view(y);
            fmpz_mpoly_q_set(yy, &y->elem.flat.data, (*U)->mctx);
        }
        else
            status |= _gr_tower_lazy_map_element(yy, &y->elem.flat.data, y->elem.flat.mctx, y->F, y->level, *U, ctx);

        /* a number field of several steps: its primitive element tower */
        if (status == GR_SUCCESS && LAZY(ctx)->options[GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT] > 0)
        {
            fmpz_mpoly_ctx_struct * old_mctx = (*U)->mctx;
            gr_tower_flat_struct * V = _gr_tower_lazy_primitive_tower(*U, ctx);
            if ((*U)->mctx != old_mctx)
            {
                /* (the search brought U's context up to date) */
                fmpz_mpoly_q_t a;
                fmpz_mpoly_q_init(a, (*U)->mctx);
                gr_tower_flat_convert(a, xx, old_mctx, *U);
                fmpz_mpoly_q_clear(xx, old_mctx);
                fmpz_mpoly_q_init(xx, (*U)->mctx);
                fmpz_mpoly_q_swap(xx, a, (*U)->mctx);
                gr_tower_flat_convert(a, yy, old_mctx, *U);
                fmpz_mpoly_q_clear(yy, old_mctx);
                fmpz_mpoly_q_init(yy, (*U)->mctx);
                fmpz_mpoly_q_swap(yy, a, (*U)->mctx);
                fmpz_mpoly_q_clear(a, (*U)->mctx);
            }
            if (V != NULL)
            {
                fmpz_mpoly_q_t a, b;
                fmpz_mpoly_q_init(a, V->mctx);
                fmpz_mpoly_q_init(b, V->mctx);
                if (_gr_tower_lazy_map_element(a, xx, (*U)->mctx, *U, (*U)->T->num_gens, V, ctx) == GR_SUCCESS &&
                    _gr_tower_lazy_map_element(b, yy, (*U)->mctx, *U, (*U)->T->num_gens, V, ctx) == GR_SUCCESS)
                {
                    fmpz_mpoly_q_clear(xx, (*U)->mctx);
                    fmpz_mpoly_q_clear(yy, (*U)->mctx);
                    *U = V;
                    fmpz_mpoly_q_init(xx, V->mctx);
                    fmpz_mpoly_q_init(yy, V->mctx);
                    fmpz_mpoly_q_swap(xx, a, V->mctx);
                    fmpz_mpoly_q_swap(yy, b, V->mctx);
                }
                fmpz_mpoly_q_clear(a, V->mctx);
                fmpz_mpoly_q_clear(b, V->mctx);
            }
        }
    }

    return status;
}
