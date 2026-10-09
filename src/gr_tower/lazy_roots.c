/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: roots of polynomials. */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* -------------------------------------------------------------------- */
/* polynomial roots                                                      */
/* -------------------------------------------------------------------- */

/* Sets res to the element x of the top field of U's tower. */
static int
_gr_tower_lazy_set_nested_top(gr_tower_lazy_elem_t res, gr_tower_flat_struct * U, gr_srcptr x, gr_ctx_t ctx)
{
    fmpz_mpoly_q_t t;
    int status;

    gr_tower_flat_ensure(U);
    fmpz_mpoly_q_init(t, U->mctx);
    status = gr_tower_flat_set_nested_at(t, x, U->T->length, U);
    if (status == GR_SUCCESS)
    {
        _gr_tower_lazy_install(res, U, t, ctx);
        res->reduced_version = 0;
    }
    fmpz_mpoly_q_clear(t, U->mctx);
    return status;
}

/*
    Roots of a squarefree polynomial f over the lazy field, appended to
    roots. All coefficients are brought to a common tower U; each root is
    expressed in U if possible and adjoined to U otherwise, and f is
    divided by the corresponding linear factor before the next root is
    looked for.
*/
/* the roots of a quadratic a x^2 + b x + c by the formula, the square
   root of the discriminant being structured (2^(1/4) for x^2 - sqrt 2) */
static int
_gr_tower_lazy_roots_quadratic(gr_vec_t roots, const gr_poly_t f, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct d, s, r, a2;
    int status = GR_SUCCESS;

    _gr_tower_lazy_init(&d, ctx);
    _gr_tower_lazy_init(&s, ctx);
    _gr_tower_lazy_init(&r, ctx);
    _gr_tower_lazy_init(&a2, ctx);

    /* d = b^2 - 4 a c, a2 = 2 a */
    status |= _gr_tower_lazy_mul(&d, gr_poly_coeff_srcptr(f, 1, ctx), gr_poly_coeff_srcptr(f, 1, ctx), ctx);
    status |= _gr_tower_lazy_mul(&s, gr_poly_coeff_srcptr(f, 2, ctx), gr_poly_coeff_srcptr(f, 0, ctx), ctx);
    status |= _gr_tower_lazy_set_si(&a2, 4, ctx);
    status |= _gr_tower_lazy_mul(&s, &s, &a2, ctx);
    status |= _gr_tower_lazy_sub(&d, &d, &s, ctx);
    status |= _gr_tower_lazy_set_si(&a2, 2, ctx);
    status |= _gr_tower_lazy_mul(&a2, &a2, gr_poly_coeff_srcptr(f, 2, ctx), ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_root_ui(&s, &d, 2, ctx);
    if (status == GR_SUCCESS)
    {
        /* (-b + s) / (2a), (-b - s) / (2a) */
        status |= _gr_tower_lazy_sub(&r, &s, gr_poly_coeff_srcptr(f, 1, ctx), ctx);
        status |= _gr_tower_lazy_div(&r, &r, &a2, ctx);
        status |= gr_vec_append(roots, &r, ctx);
        status |= _gr_tower_lazy_add(&r, &s, gr_poly_coeff_srcptr(f, 1, ctx), ctx);
        status |= _gr_tower_lazy_neg(&r, &r, ctx);
        status |= _gr_tower_lazy_div(&r, &r, &a2, ctx);
        status |= gr_vec_append(roots, &r, ctx);
    }

    _gr_tower_lazy_clear(&d, ctx);
    _gr_tower_lazy_clear(&s, ctx);
    _gr_tower_lazy_clear(&r, ctx);
    _gr_tower_lazy_clear(&a2, ctx);
    return status;
}

/*
    Roots of f which are algebraic numbers of small degree and height
    (i, roots of unity, small radicals): each numerical root is guessed
    as such and the guess verified exactly; the roots found are appended
    and divided out of f.
*/
#define LAZY_ROOTS_GUESS_DEGREE 6
#define LAZY_ROOTS_GUESS_BITS 24

static int
_gr_tower_lazy_roots_guess_algebraic(gr_vec_t roots, gr_poly_t f, gr_ctx_t ctx)
{
    slong prec = 200, n = f->length - 1, i, num, D;
    acb_poly_t fz;
    acb_ptr zs;
    int status = GR_SUCCESS;

    if (n < 3)
        return GR_SUCCESS;

    acb_poly_init(fz);
    acb_poly_fit_length(fz, n + 1);
    for (i = 0; i <= n && status == GR_SUCCESS; i++)
        status = gr_tower_lazy_get_acb(fz->coeffs + i, gr_poly_coeff_srcptr(f, i, ctx), prec, ctx);
    _acb_poly_set_length(fz, n + 1);
    if (status != GR_SUCCESS)
    {
        acb_poly_clear(fz);
        return GR_SUCCESS;   /* (no guesses) */
    }

    /* the degree of the field of the coefficients, at least (the
       largest degree of the algebraic prefixes of their towers): a root
       of degree n D over QQ is a generic root, which the exact
       factorization over the tower finds more cheaply than the
       verification of the guess (a merge with a factorization of the
       guessed minimal polynomial) */
    {
        slong j;
        D = 1;
        for (i = 0; i <= n; i++)
        {
            gr_tower_lazy_elem_struct * c = _gr_tower_lazy_flat_view(((gr_tower_lazy_elem_struct *) f->coeffs) + i);
            slong k, Dc = 1;
            k = _gr_tower_lazy_alg_level(c);
            for (j = 1; j <= k && Dc < WORD(1) << 20; j++)
                Dc *= gr_tower_step_degree(c->F->T, j);
            D = FLINT_MAX(D, Dc);
        }
    }

    zs = _acb_vec_init(n);
    num = acb_poly_find_roots(zs, fz, NULL, 0, prec);

    if (num == n)
    {
        for (i = 0; i < n && f->length >= 3; i++)
        {
            qqbar_t q;
            qqbar_init(q);
            if (acb_rel_accuracy_bits(zs + i) > 100 &&
                qqbar_guess(q, zs + i, (D == 1) ? LAZY_ROOTS_GUESS_DEGREE : FLINT_MIN(LAZY_ROOTS_GUESS_DEGREE, n * D - 1),
                    LAZY_ROOTS_GUESS_BITS, 0, prec))
            {
                gr_tower_lazy_elem_struct r, v;
                _gr_tower_lazy_init(&r, ctx);
                _gr_tower_lazy_init(&v, ctx);
                if (_gr_tower_lazy_set_qqbar_structured(&r, q, ctx) == GR_SUCCESS &&
                    gr_poly_evaluate(&v, f, &r, ctx) == GR_SUCCESS &&
                    _gr_tower_lazy_is_zero(&v, ctx) == T_TRUE)
                {
                    /* a root: divide it out */
                    gr_poly_t quo;
                    gr_poly_init(quo, ctx);
                    status |= gr_poly_div_root(quo, &v, f, &r, ctx);
                    if (status == GR_SUCCESS)
                    {
                        status |= gr_vec_append(roots, &r, ctx);
                        gr_poly_swap(f, quo, ctx);
                    }
                    gr_poly_clear(quo, ctx);
                }
                _gr_tower_lazy_clear(&r, ctx);
                _gr_tower_lazy_clear(&v, ctx);
            }
            qqbar_clear(q);
        }
    }

    _acb_vec_clear(zs, n);
    acb_poly_clear(fz);
    return status;
}

static int _gr_tower_lazy_roots_squarefree_tower(gr_vec_t roots, const gr_poly_t f, int * factor_only, gr_ctx_t ctx);

/* a x^n + b with a, b nonzero (and n >= 2) */
static int
_gr_tower_lazy_poly_is_binomial(const gr_poly_t f, gr_ctx_t ctx)
{
    slong i, n = f->length - 1;
    if (n < 2 || _gr_tower_lazy_is_zero(gr_poly_coeff_srcptr(f, 0, ctx), ctx) != T_FALSE)
        return 0;
    for (i = 1; i < n; i++)
        if (_gr_tower_lazy_is_zero(gr_poly_coeff_srcptr(f, i, ctx), ctx) != T_TRUE)
            return 0;
    return 1;
}

/* the roots of a x^n + b: root_n(-b/a) times the n-th roots of unity */
static int
_gr_tower_lazy_roots_binomial(gr_vec_t roots, const gr_poly_t f, gr_ctx_t ctx)
{
    slong k, n = f->length - 1;
    gr_tower_lazy_elem_struct r0, r, z;
    int status;

    _gr_tower_lazy_init(&r0, ctx);
    _gr_tower_lazy_init(&r, ctx);
    _gr_tower_lazy_init(&z, ctx);

    status = _gr_tower_lazy_div(&r0, gr_poly_coeff_srcptr(f, 0, ctx), gr_poly_coeff_srcptr(f, n, ctx), ctx);
    status |= _gr_tower_lazy_neg(&r0, &r0, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_root_ui(&r0, &r0, n, ctx);

    for (k = 0; k < n && status == GR_SUCCESS; k++)
    {
        if (k == 0)
            status = _gr_tower_lazy_set(&r, &r0, ctx);
        else
        {
            status = _gr_tower_lazy_root_of_unity(&z, k, n, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_mul(&r, &r0, &z, ctx);
        }
        if (status == GR_SUCCESS)
            status = gr_vec_append(roots, &r, ctx);
    }

    _gr_tower_lazy_clear(&r0, ctx);
    _gr_tower_lazy_clear(&r, ctx);
    _gr_tower_lazy_clear(&z, ctx);
    return status;
}

static int
_gr_tower_lazy_roots_squarefree(gr_vec_t roots, const gr_poly_t f_in, gr_ctx_t ctx)
{
    gr_poly_t f;
    int status;

    if (f_in->length - 1 < 1)
        return GR_SUCCESS;

    gr_poly_init(f, ctx);
    status = gr_poly_set(f, f_in, ctx);

    /* algebraic roots of small degree (i, say) are found without
       extending the tower */
    if (status == GR_SUCCESS && f->length - 1 >= 3)
        status = _gr_tower_lazy_roots_guess_algebraic(roots, f, ctx);

    /* exact factorization over the field of the coefficients, when it
       applies, splits f into factors whose roots are found by the
       structured methods */
    if (status == GR_SUCCESS && f->length - 1 >= 3 && !_gr_tower_lazy_poly_is_binomial(f, ctx))
    {
        int factor_only = 1;
        status = _gr_tower_lazy_roots_squarefree_tower(roots, f, &factor_only, ctx);
        if (status == GR_SUCCESS && factor_only)
        {
            gr_poly_clear(f, ctx);
            return GR_SUCCESS;
        }
    }

    if (status == GR_SUCCESS)
    {
        if (f->length - 1 == 1)
        {
            gr_tower_lazy_elem_struct r;
            _gr_tower_lazy_init(&r, ctx);
            status = _gr_tower_lazy_div(&r, gr_poly_coeff_srcptr(f, 0, ctx), gr_poly_coeff_srcptr(f, 1, ctx), ctx);
            status |= _gr_tower_lazy_neg(&r, &r, ctx);
            status |= gr_vec_append(roots, &r, ctx);
            _gr_tower_lazy_clear(&r, ctx);
        }
        else if (f->length - 1 == 2)
            status = _gr_tower_lazy_roots_quadratic(roots, f, ctx);
        else if (f->length - 1 >= 3 && _gr_tower_lazy_poly_is_binomial(f, ctx))
            status = _gr_tower_lazy_roots_binomial(roots, f, ctx);
        else if (f->length - 1 >= 3)
            status = _gr_tower_lazy_roots_squarefree_tower(roots, f, NULL, ctx);
    }

    gr_poly_clear(f, ctx);
    return status;
}

static int _gr_tower_lazy_roots_squarefree(gr_vec_t roots, const gr_poly_t f_in, gr_ctx_t ctx);

/* degree limit ([F : F_0] deg g) for exact factorization in root finding
   (smaller over rational function fields, where the norms are factored
   as multivariate polynomials) */

/* After a new algebraic step on top of T (top: the previous top field):
   its definition, and g (over top) promoted to the new top field */
static int
_promote_to_new_top(gr_poly_t g, gr_ctx_struct * top, gr_tower_t T, gr_ctx_t ctx)
{
    gr_poly_t g2;
    gr_ctx_struct * newtop = gr_tower_field(T);
    slong i;
    int status = GR_SUCCESS;

    _gr_tower_lazy_new_def(GR_TOWER_STEP(T, T->length - 1), T, ctx);

    gr_poly_init(g2, newtop);
    gr_poly_fit_length(g2, g->length, newtop);
    for (i = 0; i < g->length && status == GR_SUCCESS; i++)
        status |= gr_tower_promote(gr_poly_coeff_ptr(g2, i, newtop), gr_poly_coeff_srcptr(g, i, top), T->length - 1, T->length, T);
    _gr_poly_set_length(g2, g->length, newtop);
    gr_poly_clear(g, top);
    gr_poly_init(g, newtop);
    gr_poly_swap(g, g2, newtop);
    gr_poly_clear(g2, newtop);
    return status;
}

/*
    One step of root finding by exact factorization of the monic
    squarefree g over the top field of U->T (whose steps are proven, or
    provable): the roots of the linear factors are appended to roots,
    and g is replaced by the product of the other factors. If that is
    not constant, a root of its irreducible factor of least degree is
    adjoined as a proven step, g is promoted to the new top field and
    *r (a heap element of the old top field, reallocated) is set to the
    new generator. Returns 1 in that case, 2 if g became constant,
    0 if the factorization does not apply (nothing changed) and -1 on
    error.
*/
static int
_gr_tower_lazy_roots_factor_step(gr_vec_t roots, gr_poly_t g, gr_ptr * r, gr_tower_flat_struct * U, gr_ctx_t ctx)
{
    gr_tower_struct * T = U->T;
    gr_ctx_struct * top = gr_tower_field(T);
    gr_ctx_t pctx;
    gr_vec_t fac;
    fmpz_vec_t fexp;
    gr_ptr c;
    slong i, D = gr_tower_degree(T), best = -1;
    int status, result = 0;

    if (!GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_FACTOR))
        return 0;
    /* (with transcendental generators, the factorization over the function
       field is costlier: a smaller fraction of the degree limit) */
    if (D * (g->length - 1) > ((T->num_trans == 0) ? GR_TOWER_OPTION(T, GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT) : GR_TOWER_OPTION(T, GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT) * 3 / 8))
        return 0;
    /* (the proofs are attempted once per version of the moduli; they
       may refine the steps in place, which leaves g valid) */
    if (!gr_tower_prove_modular(T, GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES)))
        return 0;
    top = gr_tower_field(T);

    gr_ctx_init_gr_poly(pctx, top);
    gr_vec_init(fac, 0, pctx);
    fmpz_vec_init(fexp, 0);
    GR_TMP_INIT(c, top);

    status = gr_tower_poly_factor_limit(c, fac, fexp, g, T->length, GR_TOWER_OPTION(T, GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT), T);

    if (status == GR_SUCCESS && gr_tower_field(T) == top)
    {
        gr_poly_t rest;
        gr_poly_init(rest, top);
        status |= gr_poly_one(rest, top);

        for (i = 0; i < fac->length && status == GR_SUCCESS; i++)
        {
            gr_poly_struct * h = gr_vec_entry_ptr(fac, i, pctx);

            if (h->length == 2)
            {
                gr_tower_lazy_elem_struct e;
                _gr_tower_lazy_init(&e, ctx);
                status |= gr_neg(*r, gr_poly_coeff_srcptr(h, 0, top), top);
                status |= _gr_tower_lazy_set_nested_top(&e, U, *r, ctx);
                if (status == GR_SUCCESS)
                    status |= gr_vec_append(roots, &e, ctx);
                _gr_tower_lazy_clear(&e, ctx);
            }
            else
            {
                status |= gr_poly_mul(rest, rest, h, top);
                if (best == -1 || h->length < ((gr_poly_struct *) gr_vec_entry_ptr(fac, best, pctx))->length)
                    best = i;
            }
        }

        /* several nonlinear factors, or a single one of degree 2 or a
           binomial: the roots of each factor are found separately, by
           the structured methods where they apply (a new step for each
           irreducible factor of degree >= 3, through this function) */
        if (status == GR_SUCCESS && best != -1)
        {
            slong num = 0, nz, j;
            gr_poly_struct * h = gr_vec_entry_ptr(fac, best, pctx);

            for (i = 0; i < fac->length; i++)
                num += (((gr_poly_struct *) gr_vec_entry_ptr(fac, i, pctx))->length > 2);
            for (nz = 0, j = 0; j < h->length; j++)
                nz += (gr_is_zero(gr_poly_coeff_srcptr(h, j, top), top) != T_TRUE);

            if (num >= 2 || h->length == 3 || nz == 2)
            {
                /* (all the factors are converted before the recursive
                   calls, which may extend the tower) */
                gr_poly_struct * hl = flint_malloc(sizeof(gr_poly_struct) * fac->length);

                for (i = 0; i < fac->length; i++)
                {
                    h = gr_vec_entry_ptr(fac, i, pctx);
                    gr_poly_init2(hl + i, h->length, ctx);
                    if (h->length <= 2)
                        continue;
                    for (j = 0; j < h->length && status == GR_SUCCESS; j++)
                        status |= _gr_tower_lazy_set_nested_top(gr_poly_coeff_ptr(hl + i, j, ctx), U, gr_poly_coeff_srcptr(h, j, top), ctx);
                    _gr_poly_set_length(hl + i, h->length, ctx);
                }

                for (i = 0; i < fac->length && status == GR_SUCCESS; i++)
                    if (hl[i].length > 2)
                        status = _gr_tower_lazy_roots_squarefree(roots, hl + i, ctx);

                for (i = 0; i < fac->length; i++)
                    gr_poly_clear(hl + i, ctx);
                flint_free(hl);
                best = -1;
                status |= gr_poly_one(rest, top);
            }
        }

        if (status == GR_SUCCESS && best == -1)
        {
            gr_poly_swap(g, rest, top);
            result = 2;

            /* (the recursive calls may have extended the tower: g and *r
               are moved to the current top field) */
            if (gr_tower_field(T) != top)
            {
                gr_poly_clear(g, top);
                gr_poly_init(g, gr_tower_field(T));
                status |= gr_poly_one(g, gr_tower_field(T));
                gr_heap_clear(*r, top);
                *r = gr_heap_init(gr_tower_field(T));
            }
        }
        else if (status == GR_SUCCESS)
        {
            gr_poly_struct * h = gr_vec_entry_ptr(fac, best, pctx);
            slong m = h->length - 1, prec;
            acb_ptr zs = _acb_vec_init(m);
            int isolated = 0;

            gr_poly_swap(g, rest, top);

            for (prec = GR_TOWER_DEFAULT_PREC; prec <= 4096 && !isolated; prec *= 2)
            {
                acb_poly_t hz;
                acb_poly_init(hz);
                if (_gr_tower_poly_get_acb_poly(hz, h, T->length, prec, T) == GR_SUCCESS)
                    isolated = (acb_poly_find_roots(zs, hz, NULL, 0, prec) == m);
                acb_poly_clear(hz);
            }

            if (!isolated)
                status = GR_UNABLE;
            else
                status = gr_tower_adjoin_algebraic(T, h, zs, GR_TOWER_STATUS_PROVEN, NULL);

            if (status == GR_SUCCESS)
            {
                gr_ctx_struct * newtop = gr_tower_field(T);
                status = _promote_to_new_top(g, top, T, ctx);

                gr_heap_clear(*r, top);
                *r = gr_heap_init(newtop);
                status |= gr_gen(*r, newtop);
                result = 1;
            }

            _acb_vec_clear(zs, m);
        }

        gr_poly_clear(rest, top);
        if (status != GR_SUCCESS)
            result = -1;
    }
    else if (status == GR_SUCCESS)
    {
        /* (the proofs refined the tower: g is still valid, the numeric
           path takes over) */
    }

    GR_TMP_CLEAR(c, top);
    fmpz_vec_clear(fexp);
    gr_vec_clear(fac, pctx);
    gr_ctx_clear(pctx);
    return result;
}

/*
    The polynomial f (of degree n >= 1) over the top field of a tower U
    of its coefficients, made monic: *Uout = U and g initialized over the
    top field (on success; *Uout = NULL when the coefficients have no
    common tower, g is then not initialized). When the coefficients
    involve only some of the generators of U (a tower shared with earlier
    computations: the roots of other polynomials, unrelated radicals), U
    is a tower of those generators (with the generators their definitions
    involve), so that neither the factorizations nor the new steps
    involve the others.
*/
static int
_gr_tower_lazy_roots_prepare(gr_tower_flat_struct ** Uout, gr_poly_t g, const gr_poly_t f, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_flat_struct * U;
    gr_tower_struct * T;
    gr_ctx_struct * top;
    slong i, n = f->length - 1;
    int status = GR_SUCCESS;

    *Uout = NULL;

    /* common tower of the coefficients */
    {
        gr_tower_lazy_elem_struct * c = _gr_tower_lazy_flat_view(f->coeffs);
        slong lu;
        U = c->F;
        lu = c->level;
        for (i = 1; i <= n; i++)
        {
            c = _gr_tower_lazy_flat_view(((gr_tower_lazy_elem_struct *) f->coeffs) + i);
            if (c->F == U)
                lu = FLINT_MAX(lu, c->level);
            else
            {
                U = _gr_tower_lazy_common_tower(U, lu, c->F, c->level, ctx);
                if (U == NULL)
                    return GR_UNABLE;
                lu = U->T->num_gens;
            }
        }
        if (U == L->trivial)
            U = _gr_tower_lazy_new_tower(ctx);
    }

    /* the coefficients in U; when they involve only some of the
       generators of U (a tower shared with earlier computations: the
       roots of other polynomials, unrelated radicals), the roots are
       sought over a tower of those generators (with the generators their
       definitions involve), so that neither the factorizations nor the
       new steps involve the others */
    {
        fmpz_mpoly_q_struct * t = flint_malloc(sizeof(fmpz_mpoly_q_struct) * (n + 1));
        int * mark = flint_calloc(FLINT_MAX(U->T->num_gens, 1), sizeof(int));
        slong * order_map = flint_malloc(sizeof(slong) * FLINT_MAX(U->T->num_gens, 1));
        slong need = 0;

        gr_tower_flat_ensure(U);
        for (i = 0; i <= n; i++)
        {
            gr_tower_lazy_elem_struct * c = _gr_tower_lazy_flat_view(((gr_tower_lazy_elem_struct *) f->coeffs) + i);
            fmpz_mpoly_q_init(t + i, U->mctx);
            if (status == GR_SUCCESS)
            {
                if (c->F == U)
                    fmpz_mpoly_q_set(t + i, &c->elem.flat.data, U->mctx);
                else
                    status |= _gr_tower_lazy_map_element(t + i, &c->elem.flat.data, c->elem.flat.mctx, c->F, c->level, U, ctx);
            }
            if (status == GR_SUCCESS)
                _gr_tower_lazy_mark_used_gens(mark, t + i, U->mctx, U->T);
        }

        _gr_tower_involved_gens_closure(mark, U->T);
        for (i = 0; i < U->T->num_gens; i++)
        {
            order_map[i] = mark[i] ? need : -1;
            need += (mark[i] != 0);
        }

        if (status == GR_SUCCESS && need < U->T->num_gens)
        {
            gr_tower_flat_struct * G = _gr_tower_lazy_new_tower(ctx);
            gr_tower_set_subset(G->T, U->T, mark);
            gr_tower_flat_ensure(G);
            for (i = 0; i <= n; i++)
            {
                fmpz_mpoly_q_t u;
                fmpz_mpoly_q_init(u, G->mctx);
                _gr_tower_flat_transport_map(u, t + i, U, G, order_map);
                fmpz_mpoly_q_clear(t + i, U->mctx);
                fmpz_mpoly_q_init(t + i, G->mctx);
                fmpz_mpoly_q_swap(t + i, u, G->mctx);
                fmpz_mpoly_q_clear(u, G->mctx);
            }
            U = G;
        }

        T = U->T;
        top = gr_tower_field(T);

        /* the polynomial over the top field, made monic */
        gr_poly_init(g, top);
        gr_poly_fit_length(g, n + 1, top);
        for (i = 0; i <= n; i++)
        {
            if (status == GR_SUCCESS)
                status |= gr_tower_flat_get_nested_at(gr_poly_coeff_ptr(g, i, top), t + i, T->length, U);
            fmpz_mpoly_q_clear(t + i, U->mctx);
        }
        _gr_poly_set_length(g, n + 1, top);
        flint_free(t);
        flint_free(mark);
        flint_free(order_map);
    }

    if (status == GR_SUCCESS)
    {
        gr_ptr lc;
        GR_TMP_INIT(lc, top);
        status |= gr_tower_inv(lc, gr_poly_coeff_srcptr(g, n, top), T);
        if (status == GR_SUCCESS)
            status |= _gr_vec_mul_scalar(g->coeffs, g->coeffs, n + 1, lc, top);
        GR_TMP_CLEAR(lc, top);
    }

    if (status != GR_SUCCESS)
    {
        gr_poly_clear(g, top);
        return status;
    }

    *Uout = U;
    return GR_SUCCESS;
}

static int
_gr_tower_lazy_roots_squarefree_tower(gr_vec_t roots, const gr_poly_t f, int * factor_only, gr_ctx_t ctx)
{
    gr_tower_flat_struct * U;
    gr_tower_struct * T;
    gr_ctx_struct * top;
    gr_poly_t g;
    slong i, n = f->length - 1;
    slong roots_len0 = roots->length;
    int status = GR_SUCCESS;

    if (n < 1)
        return GR_SUCCESS;

    status = _gr_tower_lazy_roots_prepare(&U, g, f, ctx);
    if (U == NULL)
    {
        if (status == GR_UNABLE && factor_only != NULL)
        {
            *factor_only = 0;
            return GR_SUCCESS;
        }
        return status;
    }
    T = U->T;
    top = gr_tower_field(T);

    while (status == GR_SUCCESS && g->length >= 2)
    {
        int first = (roots_len0 == roots->length && g->length == n + 1);
        slong m = g->length - 1;
        gr_ptr r;
        acb_t z;
        int found = 0;

        /* (g lives in the top field of T as it was when g was made or
           promoted; a rebuild of the nested chain by one of the
           operations below -- a reordering of the generators by a
           relation found during a zero test, say -- would invalidate
           it: the zero tests run on copies of the tower, so this is
           not expected, but checked) */
        if (gr_tower_field(T) != top)
        {
            status = GR_UNABLE;
            _gr_poly_set_length(g, 0, top);   /* (the coefficients are lost with their contexts) */
            break;
        }

        r = gr_heap_init(top);
        acb_init(z);

        if (m == 1)
        {
            status |= gr_neg(r, gr_poly_coeff_srcptr(g, 0, top), top);
            found = 1;
        }
        else if (m == 2 && !first)
        {
            /* a quadratic cofactor: by the formula (a structured square
               root of the discriminant, rather than a lattice search in
               the tower followed by a new step) */
            gr_poly_t h;
            gr_poly_init2(h, 3, ctx);
            for (i = 0; i < 3 && status == GR_SUCCESS; i++)
                status |= _gr_tower_lazy_set_nested_top(gr_poly_coeff_ptr(h, i, ctx), U, gr_poly_coeff_srcptr(g, i, top), ctx);
            _gr_poly_set_length(h, 3, ctx);
            /* (g and r are released first: the square root may extend
               the tower in place) */
            gr_heap_clear(r, top);
            acb_clear(z);
            gr_poly_clear(g, top);
            gr_poly_init(g, top);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_roots_quadratic(roots, h, ctx);
            gr_poly_clear(h, ctx);
            break;
        }
        else if ((found = _gr_tower_lazy_roots_factor_step(roots, g, &r, U, ctx)) == 0 && factor_only != NULL && first)
        {
            /* (factorization does not apply: nothing done) */
            *factor_only = 0;
            gr_heap_clear(r, top);
            acb_clear(z);
            break;
        }
        else if (found != 0)
        {
            /* (roots of the linear factors appended, g replaced by the
               product of the other factors; a root of one of them
               adjoined as a proven step, g promoted, r its generator) */
            if (found < 0)
            {
                status = GR_UNABLE;
                found = 0;
            }
            else if (found == 2)
                found = 0;   /* (only linear factors: g is now constant) */
            top = gr_tower_field(T);
        }
        else
        {
            /* isolate one root numerically */
            slong prec;
            acb_ptr zs = _acb_vec_init(m);
            int isolated = 0;

            for (prec = GR_TOWER_DEFAULT_PREC; prec <= 4096 && !isolated; prec *= 2)
            {
                acb_poly_t gz;
                acb_poly_init(gz);
                if (_gr_tower_poly_get_acb_poly(gz, g, T->length, prec, T) == GR_SUCCESS)
                {
                    slong num = acb_poly_find_roots(zs, gz, NULL, 0, prec);
                    if (num == m)
                        isolated = 1;
                }
                acb_poly_clear(gz);
            }

            if (isolated)
            {
                slong j;

                /* a root lying in the tower, if any, first (unless the
                   reduction of g at a place of the tower has no root,
                   which excludes roots in the tower and saves the
                   lattice searches) */
                int none = gr_tower_poly_no_roots_modular(g, T, GR_TOWER_OPTION(T, GR_TOWER_OPT_NO_ROOTS_TRIES));

                for (j = 0; j < m && !found && !none; j++)
                {
                    acb_set(z, zs + j);
                    if (gr_tower_express_limit(r, g, z, GR_TOWER_MERGE_EXPRESS_PREC(T, gr_tower_degree(T)), T) == GR_SUCCESS)
                        found = 1;
                }
                if (!found)
                    acb_set(z, zs);

                if (found)
                {
                }
                else
                {
                    status = gr_tower_adjoin_algebraic(T, g, z, GR_TOWER_STATUS_DYNAMIC, NULL);
                    if (status == GR_SUCCESS)
                    {
                        gr_ctx_struct * newtop = gr_tower_field(T);
                        /* (the root is the generator) */
                        status = _promote_to_new_top(g, top, T, ctx);

                        gr_heap_clear(r, top);
                        top = newtop;
                        r = gr_heap_init(top);
                        status |= gr_gen(r, top);
                        found = 1;
                    }
                }
            }
            else
            {
                status = GR_UNABLE;
            }

            _acb_vec_clear(zs, m);
        }

        if (status == GR_SUCCESS && found)
        {
            gr_tower_lazy_elem_struct e;

            _gr_tower_lazy_init(&e, ctx);
            status |= _gr_tower_lazy_set_nested_top(&e, U, r, ctx);
            if (status == GR_SUCCESS)
                status |= gr_vec_append(roots, &e, ctx);
            _gr_tower_lazy_clear(&e, ctx);

            /* g = g / (x - r) */
            if (status == GR_SUCCESS)
            {
                gr_ptr rem;
                GR_TMP_INIT(rem, top);
                status |= gr_poly_div_root(g, rem, g, r, top);
                GR_TMP_CLEAR(rem, top);
            }
        }

        gr_heap_clear(r, top);
        acb_clear(z);
    }

    gr_poly_clear(g, top);
    return status;
}

/*
    res = the root of the squarefree polynomial f (of degree >= 1, with
    coefficients in the field) which overlaps ref, or (with pm set) one
    of ref and -ref; exactly one root must match, numerically. The root
    is adjoined as an algebraic step of the tower of the coefficients
    (with dynamic status: the polynomial may be reducible over it),
    without the factorizations and lattice searches of the general root
    finding: for callers which know that the root is generically new
    (the values of modular functions at a new point, from a modular
    equation). Copies of the step made by separate calls are identified
    when the towers are merged. GR_UNABLE when the roots cannot be
    isolated or the match is not unique.
*/
int
_gr_tower_lazy_poly_root_near(gr_tower_lazy_elem_t res, const gr_poly_t f, const acb_t ref, int pm, gr_ctx_t ctx)
{
    gr_tower_flat_struct * U;
    gr_tower_struct * T;
    gr_ctx_struct * top;
    gr_poly_t g;
    acb_ptr zs;
    acb_t mref;
    slong j, m, prec, found = -1;
    int status, isolated = 0;

    if (f->length < 2)
        return GR_DOMAIN;

    status = _gr_tower_lazy_roots_prepare(&U, g, f, ctx);
    if (U == NULL)
        return (status == GR_SUCCESS) ? GR_UNABLE : status;

    T = U->T;
    top = gr_tower_field(T);
    m = g->length - 1;

    if (m == 1)
    {
        gr_ptr r;
        GR_TMP_INIT(r, top);
        status = gr_neg(r, gr_poly_coeff_srcptr(g, 0, top), top);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_set_nested_top(res, U, r, ctx);
        GR_TMP_CLEAR(r, top);
        gr_poly_clear(g, top);
        return status;
    }

    zs = _acb_vec_init(m);
    acb_init(mref);
    acb_neg(mref, ref);

    for (prec = GR_TOWER_DEFAULT_PREC; prec <= 4096 && !isolated; prec *= 2)
    {
        acb_poly_t gz;
        acb_poly_init(gz);
        if (_gr_tower_poly_get_acb_poly(gz, g, T->length, prec, T) == GR_SUCCESS)
            isolated = (acb_poly_find_roots(zs, gz, NULL, 0, prec) == m);
        acb_poly_clear(gz);
    }

    if (isolated)
    {
        for (j = 0; j < m; j++)
        {
            if (acb_overlaps(zs + j, ref) || (pm && acb_overlaps(zs + j, mref)))
            {
                if (found != -1)
                {
                    found = -2;
                    break;
                }
                found = j;
            }
        }
    }

    if (found < 0)
        status = GR_UNABLE;
    else
        status = gr_tower_adjoin_algebraic(T, g, zs + found, GR_TOWER_STATUS_DYNAMIC, NULL);

    if (status == GR_SUCCESS)
    {
        _gr_tower_lazy_new_def(GR_TOWER_STEP(T, T->length - 1), T, ctx);
        _gr_tower_lazy_set_gen(res, U, T->length, ctx);
    }

    _acb_vec_clear(zs, m);
    acb_clear(mref);
    gr_poly_clear(g, top);
    return status;
}

/*
    Appends the roots of the irreducible primitive integer polynomial p
    to roots. Linear and quadratic polynomials are solved by radicals
    (through the structured roots); for degree d >= 3 the roots are
    algebraic numbers.
*/
static int
_gr_tower_lazy_roots_irreducible(gr_vec_t roots, const fmpz_poly_t p, gr_ctx_t ctx)
{
    slong d = fmpz_poly_degree(p), i;
    gr_tower_lazy_elem_struct e, t;
    int status = GR_SUCCESS;

    _gr_tower_lazy_init(&e, ctx);
    _gr_tower_lazy_init(&t, ctx);

    if (d == 1)
    {
        fmpq_t q;
        fmpq_init(q);
        fmpz_neg(fmpq_numref(q), p->coeffs + 0);
        fmpz_set(fmpq_denref(q), p->coeffs + 1);
        fmpq_canonicalise(q);
        status = _gr_tower_lazy_set_fmpq(&e, q, ctx);
        status |= gr_vec_append(roots, &e, ctx);
        fmpq_clear(q);
    }
    else if (d == 2)
    {
        /* (-b +/- sqrt(b^2 - 4ac)) / (2a) */
        fmpq_t disc, q;
        fmpq_init(disc);
        fmpq_init(q);
        fmpz_mul(fmpq_numref(disc), p->coeffs + 1, p->coeffs + 1);
        fmpz_submul(fmpq_numref(disc), p->coeffs + 2, p->coeffs + 0);
        fmpz_submul(fmpq_numref(disc), p->coeffs + 2, p->coeffs + 0);
        fmpz_submul(fmpq_numref(disc), p->coeffs + 2, p->coeffs + 0);
        fmpz_submul(fmpq_numref(disc), p->coeffs + 2, p->coeffs + 0);
        status = _gr_tower_lazy_root_fmpq(&t, disc, 2, ctx);
        for (i = 0; i < 2 && status == GR_SUCCESS; i++)
        {
            fmpz_neg(fmpq_numref(q), p->coeffs + 1);
            fmpz_mul_ui(fmpq_denref(q), p->coeffs + 2, 2);
            fmpq_canonicalise(q);
            status |= _gr_tower_lazy_set_fmpq(&e, q, ctx);
            fmpq_one(q);
            fmpz_mul_ui(fmpq_denref(q), p->coeffs + 2, 2);
            fmpq_canonicalise(q);
            if (i == 0)
                status |= gr_addmul_fmpq(&e, &t, q, ctx);
            else
                status |= gr_submul_fmpq(&e, &t, q, ctx);
            status |= gr_vec_append(roots, &e, ctx);
        }
        fmpq_clear(disc);
        fmpq_clear(q);
    }
    else
    {
        /* each root is an algebraic number of its own (hash-consed by
           value); when conjugate roots meet in one tower, the merge
           divides the conjugates already present out of the minimal
           polynomial (see _divide_out_conjugates in map.c), so that the
           tower of all the roots has the Cauchy moduli and symmetric
           functions of the roots reduce to the coefficients */
        qqbar_ptr rts = _qqbar_vec_init(d);
        qqbar_roots_fmpz_poly(rts, p, QQBAR_ROOTS_IRREDUCIBLE);
        int binomial = 1;
        for (i = 1; i < d && binomial; i++)
            binomial = fmpz_is_zero(p->coeffs + i);
        for (i = 0; i < d && status == GR_SUCCESS; i++)
        {
            /* (binomials: principal radicals times roots of unity) */
            status = binomial ? _gr_tower_lazy_set_qqbar_structured(&e, rts + i, ctx) : _gr_tower_lazy_set_qqbar(&e, rts + i, ctx);
            status |= gr_vec_append(roots, &e, ctx);
        }
        _qqbar_vec_clear(rts, d);
    }

    _gr_tower_lazy_clear(&e, ctx);
    _gr_tower_lazy_clear(&t, ctx);
    return status;
}

/*
    The roots of a squarefree f whose coefficients are rational functions
    of transcendental generators only (no algebraic generator, all in one
    tower): the linear factors of f as a polynomial over Q[t_1, ..., t_r]
    (multivariate factorization). The roots found are appended; *rest is
    set to f divided by them (the factors of higher degree, for the
    general method). Returns 0 if this does not apply (rest untouched).
    The characteristic polynomial of log(M) for M with eigenvalues a, b, c
    is (x - log a)(x - log b)(x - log c).
*/
static int
_gr_tower_lazy_roots_trans_linear(gr_vec_t roots, gr_poly_t rest, const gr_poly_t f, int * status, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_lazy_elem_struct ** cs;   /* (flat views of the coefficients) */
    gr_tower_flat_struct * F = NULL;
    gr_tower_struct * T;
    const fmpz_mpoly_ctx_struct * mctx;
    fmpz_mpoly_ctx_t mctx2;
    fmpz_mpoly_t D, g, t, P, xi;
    fmpz_mpoly_factor_t pf;
    slong * cto, * cback, r, i, v, found = 0;
    int * used, ok = 1, trans = 0;

    if (f->length < 3)
        return 0;
    cs = flint_malloc(sizeof(gr_tower_lazy_elem_struct *) * f->length);
    for (i = 0; i < f->length; i++)
    {
        cs[i] = _gr_tower_lazy_flat_view(((gr_tower_lazy_elem_struct *) f->coeffs) + i);
        if (cs[i]->F == L->trivial)
            continue;
        if (F == NULL)
            F = cs[i]->F;
        else if (cs[i]->F != F)
        {
            flint_free(cs);
            return 0;
        }
    }
    if (F == NULL)
    {
        flint_free(cs);
        return 0;
    }
    T = F->T;
    gr_tower_flat_ensure(F);
    mctx = F->mctx;
    r = mctx->minfo->nvars;

    /* (rational coefficients are in the trivial tower: their data are
       constants, valid in any context) */
    used = flint_calloc(r, sizeof(int));
    for (i = 0; i < f->length && ok; i++)
    {
        int * u = flint_malloc(sizeof(int) * r);
        if (cs[i]->F == F)
        {
            if (cs[i]->elem.flat.mctx != mctx)
                ok = 0;
            else
            {
                fmpz_mpoly_q_used_vars(u, &cs[i]->elem.flat.data, mctx);
                for (v = 0; v < r; v++)
                    used[v] |= u[v];
            }
        }
        flint_free(u);
    }
    for (v = 0; v < r && ok; v++)
    {
        slong d;
        if (!used[v])
            continue;
        for (d = 0; d < T->num_gens; d++)
            if (GR_TOWER_FLAT_VAR_D(F, d) == v)
                break;
        if (d == T->num_gens || GR_TOWER_GEN(T, d)->kind == GR_TOWER_ALGEBRAIC)
            ok = 0;
        else
            trans = 1;
    }
    flint_free(used);
    if (!ok || !trans)
    {
        flint_free(cs);
        return 0;
    }

    fmpz_mpoly_ctx_init(mctx2, r + 1, ORD_LEX);
    fmpz_mpoly_init(D, mctx);
    fmpz_mpoly_init(g, mctx);
    fmpz_mpoly_init(t, mctx);
    fmpz_mpoly_init(P, mctx2);
    fmpz_mpoly_init(xi, mctx2);
    fmpz_mpoly_factor_init(pf, mctx2);
    cto = flint_malloc(sizeof(slong) * (r + 1));
    cback = flint_malloc(sizeof(slong) * (r + 1));
    for (v = 0; v < r; v++)
        cto[v] = cback[v] = v;
    cback[r] = -1;

    /* the coefficients as fractions in mctx (constants of the trivial
       tower converted) */
    {
        fmpz_mpoly_q_struct * q = flint_malloc(sizeof(fmpz_mpoly_q_struct) * f->length);
        for (i = 0; i < f->length; i++)
        {
            fmpz_mpoly_q_init(q + i, mctx);
            if (cs[i]->F == F)
                fmpz_mpoly_q_set(q + i, &cs[i]->elem.flat.data, mctx);
            else
            {
                fmpz_t a;
                fmpz_init(a);
                fmpz_mpoly_get_fmpz(a, fmpz_mpoly_q_numref(&cs[i]->elem.flat.data), cs[i]->elem.flat.mctx);
                fmpz_mpoly_set_fmpz(fmpz_mpoly_q_numref(q + i), a, mctx);
                fmpz_mpoly_get_fmpz(a, fmpz_mpoly_q_denref(&cs[i]->elem.flat.data), cs[i]->elem.flat.mctx);
                fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(q + i), a, mctx);
                fmpz_clear(a);
            }
        }

        /* D = lcm of the denominators; P = D f in Q[t, x] */
        fmpz_mpoly_one(D, mctx);
        for (i = 0; i < f->length; i++)
        {
            if (fmpz_mpoly_q_is_zero(q + i, mctx))
                continue;
            fmpz_mpoly_gcd(g, D, fmpz_mpoly_q_denref(q + i), mctx);
            fmpz_mpoly_mul(t, D, fmpz_mpoly_q_denref(q + i), mctx);
            fmpz_mpoly_divides(D, t, g, mctx);
        }
        fmpz_mpoly_zero(P, mctx2);
        for (i = 0; i < f->length; i++)
        {
            fmpz_mpoly_t u;
            if (fmpz_mpoly_q_is_zero(q + i, mctx))
                continue;
            fmpz_mpoly_init(u, mctx2);
            fmpz_mpoly_divides(t, D, fmpz_mpoly_q_denref(q + i), mctx);
            fmpz_mpoly_mul(t, t, fmpz_mpoly_q_numref(q + i), mctx);
            fmpz_mpoly_compose_fmpz_mpoly_gen(u, t, cto, mctx, mctx2);
            fmpz_mpoly_gen(xi, r, mctx2);
            fmpz_mpoly_pow_ui(xi, xi, i, mctx2);
            fmpz_mpoly_mul(u, u, xi, mctx2);
            fmpz_mpoly_add(P, P, u, mctx2);
            fmpz_mpoly_clear(u, mctx2);
        }
        for (i = 0; i < f->length; i++)
            fmpz_mpoly_q_clear(q + i, mctx);
        flint_free(q);
    }

    if (!fmpz_mpoly_factor(pf, P, mctx2))
        ok = 0;

    if (ok)
    {
        gr_poly_t lin;
        gr_poly_init(lin, ctx);
        *status = gr_poly_set(rest, f, ctx);
        for (i = 0; i < pf->num && *status == GR_SUCCESS; i++)
        {
            slong var = r;
            ulong e0 = 0, e1 = 1;
            fmpz_mpoly_t c0, c1;
            fmpz_mpoly_q_t root;
            gr_tower_lazy_elem_struct z;

            if (fmpz_mpoly_degree_si(pf->poly + i, r, mctx2) != 1)
                continue;
            fmpz_mpoly_init(c0, mctx2);
            fmpz_mpoly_init(c1, mctx2);
            fmpz_mpoly_q_init(root, mctx);
            fmpz_mpoly_get_coeff_vars_ui(c0, pf->poly + i, &var, &e0, 1, mctx2);
            fmpz_mpoly_get_coeff_vars_ui(c1, pf->poly + i, &var, &e1, 1, mctx2);
            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_numref(root), c0, cback, mctx2, mctx);
            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_denref(root), c1, cback, mctx2, mctx);
            fmpz_mpoly_neg(fmpz_mpoly_q_numref(root), fmpz_mpoly_q_numref(root), mctx);
            fmpz_mpoly_q_canonicalise(root, mctx);
            _gr_tower_lazy_init(&z, ctx);
            _gr_tower_lazy_set_flat(&z, F, root, mctx, ctx);
            *status |= gr_vec_append(roots, &z, ctx);
            /* rest /= (x - z) */
            *status |= gr_poly_set_coeff_si(lin, 1, 1, ctx);
            *status |= gr_neg(&z, &z, ctx);
            *status |= gr_poly_set_coeff_scalar(lin, 0, &z, ctx);
            *status |= gr_poly_divexact(rest, rest, lin, ctx);
            found++;
            _gr_tower_lazy_clear(&z, ctx);
            fmpz_mpoly_q_clear(root, mctx);
            fmpz_mpoly_clear(c0, mctx2);
            fmpz_mpoly_clear(c1, mctx2);
        }
        gr_poly_clear(lin, ctx);
    }

    fmpz_mpoly_factor_clear(pf, mctx2);
    fmpz_mpoly_clear(D, mctx);
    fmpz_mpoly_clear(g, mctx);
    fmpz_mpoly_clear(t, mctx);
    fmpz_mpoly_clear(P, mctx2);
    fmpz_mpoly_clear(xi, mctx2);
    fmpz_mpoly_ctx_clear(mctx2);
    flint_free(cto);
    flint_free(cback);
    flint_free(cs);
    return ok && found > 0;
}

int
_gr_tower_lazy_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    gr_ptr c;
    slong i;
    int status;

    if (poly->length == 0)
        return GR_DOMAIN;

    gr_vec_set_length(roots, 0, ctx);
    fmpz_vec_set_length(mult, 0);

    /* rational coefficients: use qqbar (irreducible factors, hash-consed) */
    {
        int rational = 1;
        gr_tower_lazy_elem_struct * cs = poly->coeffs;
        for (i = 0; i < poly->length; i++)
            if (cs[i].F != L->trivial)
                rational = 0;

        if (rational)
        {
            /* factor over QQ; each irreducible factor is handled by
               _gr_tower_lazy_roots_irreducible */
            fmpq_poly_t fq;
            fmpz_poly_t fz;
            fmpz_poly_factor_t fac;

            fmpq_poly_init(fq);
            fmpz_poly_init(fz);
            fmpz_poly_factor_init(fac);

            fmpq_poly_fit_length(fq, poly->length);
            for (i = 0; i < poly->length; i++)
            {
                fmpq_t q;
                fmpq_init(q);
                GR_MUST_SUCCEED(_gr_tower_lazy_get_fmpq(q, cs + i, ctx));
                fmpq_poly_set_coeff_fmpq(fq, i, q);
                fmpq_clear(q);
            }
            fmpq_poly_get_numerator(fz, fq);

            status = GR_SUCCESS;
            if (fmpz_poly_is_zero(fz))
                status = GR_DOMAIN;
            else
            {
                fmpz_poly_factor(fac, fz);
                for (i = 0; i < fac->num && status == GR_SUCCESS; i++)
                {
                    slong before = roots->length, j;
                    status = _gr_tower_lazy_roots_irreducible(roots, fac->p + i, ctx);
                    for (j = before; j < roots->length; j++)
                    {
                        fmpz_t e;
                        fmpz_init_set_ui(e, fac->exp[i]);
                        fmpz_vec_append(mult, e);
                        fmpz_clear(e);
                    }
                }
            }

            fmpz_poly_factor_clear(fac);
            fmpz_poly_clear(fz);
            fmpq_poly_clear(fq);
            return status;
        }
    }

    c = gr_heap_init(ctx);
    gr_poly_vec_init(fac, 0, ctx);
    fmpz_vec_init(exp, 0);

    status = gr_poly_factor_squarefree(c, fac, exp, poly, ctx);

    for (i = 0; i < fac->length && status == GR_SUCCESS; i++)
    {
        gr_poly_struct * fac_i = gr_poly_vec_entry_ptr(fac, i, ctx);
        slong before = roots->length, j;

        {
            gr_poly_t rest;
            gr_poly_init(rest, ctx);
            if (_gr_tower_lazy_roots_trans_linear(roots, rest, fac_i, &status, ctx))
            {
                if (status == GR_SUCCESS && rest->length > 1)
                    status |= _gr_tower_lazy_roots_squarefree(roots, rest, ctx);
            }
            else
                status |= _gr_tower_lazy_roots_squarefree(roots, fac_i, ctx);
            gr_poly_clear(rest, ctx);
        }
        for (j = before; j < roots->length; j++)
            fmpz_vec_append(mult, exp->entries + i);
    }

    gr_poly_vec_clear(fac, ctx);
    fmpz_vec_clear(exp);
    gr_heap_clear(c, ctx);

    return status;
}

POP_OPTIONS
