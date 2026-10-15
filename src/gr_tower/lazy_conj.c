/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: conversions to rational numbers, complex conjugation and transfer between contexts. */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* -------------------------------------------------------------------- */
/* conversions to rationals                                              */
/* -------------------------------------------------------------------- */

/*
    Reduced representations are canonical modulo the triangular set, so a
    rational number is represented by a constant unless a modulus is
    reducible (dynamic step which has not been refined yet): then a
    non-constant representation gives GR_UNABLE rather than GR_DOMAIN.
*/
static int _gr_tower_lazy_get_fmpq_inner(fmpq_t res, const gr_tower_lazy_elem_t x_in, gr_ctx_t ctx);

int
_gr_tower_lazy_get_fmpq(fmpq_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status;
    if (x->shallow)
    {
        gr_tower_lazy_elem_struct t;
        _gr_tower_lazy_init(&t, ctx);
        status = _gr_tower_lazy_set(&t, x, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_get_fmpq_inner(res, &t, ctx);
        _gr_tower_lazy_clear(&t, ctx);
        return status;
    }
    return _gr_tower_lazy_get_fmpq_inner(res, x, ctx);
}

static int
_gr_tower_lazy_get_fmpq_inner(fmpq_t res, const gr_tower_lazy_elem_t x_in, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    gr_tower_struct * T;

    x = _gr_tower_lazy_flat_view(x);
    GR_MUST_SUCCEED(gr_tower_flat_reduce(&x->elem.flat.data, x->F));

    if (fmpz_mpoly_q_get_fmpq(res, &x->elem.flat.data, x->elem.flat.mctx))
    {
        return GR_SUCCESS;
    }

    T = x->F->T;
    if (gr_tower_prove_modular(T, GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES)))
        return GR_DOMAIN;

    /* (only the generators x involves matter: a conjecturally
       transcendental generator elsewhere in the tower does not) */
    {
        int * mark = flint_calloc(FLINT_MAX(T->num_gens, 1), sizeof(int));
        slong d;
        int proven = 1;

        gr_tower_flat_ensure(x->F);
        _gr_tower_flat_involved(mark, &x->elem.flat.data, x->F);
        for (d = 0; d < T->num_gens && proven; d++)
            if (mark[d] && T->gens[d].status != GR_TOWER_STATUS_PROVEN)
                proven = 0;
        flint_free(mark);

        return proven ? GR_DOMAIN : GR_UNABLE;
    }
}

int
_gr_tower_lazy_get_fmpz(fmpz_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    fmpq_t q;
    int status;

    fmpq_init(q);
    status = _gr_tower_lazy_get_fmpq(q, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_is_one(fmpq_denref(q)))
            fmpz_set(res, fmpq_numref(q));
        else
            status = GR_DOMAIN;
    }
    fmpq_clear(q);
    return status;
}

int
_gr_tower_lazy_get_si(slong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    status = _gr_tower_lazy_get_fmpz(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_fits_si(t))
            *res = fmpz_get_si(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

int
_gr_tower_lazy_get_ui(ulong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    status = _gr_tower_lazy_get_fmpz(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(t) >= 0 && fmpz_abs_fits_ui(t))
            *res = fmpz_get_ui(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

/* -------------------------------------------------------------------- */
/* complex conjugation                                                   */
/* -------------------------------------------------------------------- */

int
_gr_tower_lazy_i(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    return _gr_tower_lazy_root_of_unity(res, 1, 4, ctx);
}

/*
    Whether the generator with definition order d is known to be real
    (enclosure with an exactly zero imaginary part, and a definition which
    is conjugation-invariant: pi, exp of a real, log of a positive real,
    a root of a polynomial over real generators).
*/
static int
_gen_is_real(gr_tower_struct * T, slong d)
{
    gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
    acb_t z;
    int real;

    acb_init(z);
    if (g->kind == GR_TOWER_ALGEBRAIC)
        real = (gr_tower_step_get_acb(z, T, g->index, GR_TOWER_DEFAULT_PREC) == GR_SUCCESS);
    else
        real = (gr_tower_trans_get_acb(z, T, g->index, GR_TOWER_DEFAULT_PREC) == GR_SUCCESS);
    real = real && arb_is_zero(acb_imagref(z));
    acb_clear(z);

    return real;
}

/* Whether the generator is i (root of x^2 + 1 with positive imaginary part). */
int
_gr_tower_lazy_gen_is_i(gr_tower_struct * T, slong d)
{
    gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
    const gr_poly_struct * m;
    gr_ctx_struct * below;

    if (g->kind != GR_TOWER_ALGEBRAIC)
        return 0;
    m = gr_tower_step_minpoly(T, g->index);
    below = gr_tower_field_at(T, g->index - 1);
    return (m->length == 3 && gr_is_one(gr_poly_coeff_srcptr(m, 0, below), below) == T_TRUE
                           && gr_is_zero(gr_poly_coeff_srcptr(m, 1, below), below) == T_TRUE
                           && arb_is_positive(acb_imagref(&g->enclosure)));
}

/* Sets res to the flat element x of F (as a lazy element). */
void
_gr_tower_lazy_set_flat(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, gr_ctx_t ctx)
{
    fmpz_mpoly_q_t t;
    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(t, F->mctx);
    gr_tower_flat_convert(t, x, x_mctx, F);
    _gr_tower_lazy_install(res, F, t, ctx);
    res->reduced_version = 0;   /* (converted, not reduced) */
    fmpz_mpoly_q_clear(t, F->mctx);
}

/*
    Where the argument u of a generator lies with respect to the negative
    real axis (the branch cut of log and of the principal roots):
    1 = on it, -1 = away from it, 0 = undecided.
*/
static int
_flat_on_cut(const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t u_mctx, gr_tower_flat_struct * F)
{
    acb_t z;
    fmpz_mpoly_q_t t;
    int res = 0;

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(t, F->mctx);
    gr_tower_flat_convert(t, u, u_mctx, F);
    acb_init(z);
    if (gr_tower_flat_get_acb(z, t, GR_TOWER_DEFAULT_PREC, F) == GR_SUCCESS)
    {
        if (arb_is_zero(acb_imagref(z)) && arb_is_negative(acb_realref(z)))
            res = 1;
        else if (!arb_contains_zero(acb_imagref(z)) || arb_is_positive(acb_realref(z)))
            res = -1;
    }
    acb_clear(z);
    fmpz_mpoly_q_clear(t, F->mctx);
    return res;
}

/* Evaluates the polynomial f (in F's context) at the lazy elements imgs
   (indexed by variable), accumulating into res. */
static int
_gr_tower_lazy_eval_poly(gr_tower_lazy_elem_t res, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_struct * fctx, slong nvars, gr_tower_lazy_elem_struct ** imgs, gr_ctx_t ctx)
{
    /* (f is read in its own context fctx, with nvars variables and images:
       the operations on the images may change the layout of the tower
       meanwhile, but the images are elements which follow) */
    slong i, v, n = fmpz_mpoly_length(f, fctx);
    ulong * exp;
    fmpz_t c;
    gr_tower_lazy_elem_struct term, t;
    int status = GR_SUCCESS;

    exp = flint_malloc(sizeof(ulong) * nvars);
    fmpz_init(c);
    _gr_tower_lazy_init(&term, ctx);
    _gr_tower_lazy_init(&t, ctx);

    status = _gr_tower_lazy_zero(res, ctx);
    for (i = 0; i < n && status == GR_SUCCESS; i++)
    {
        fmpz_mpoly_get_term_coeff_fmpz(c, f, i, fctx);
        fmpz_mpoly_get_term_exp_ui(exp, f, i, fctx);
        status = _gr_tower_lazy_set_fmpz(&term, c, ctx);
        for (v = 0; v < nvars && status == GR_SUCCESS; v++)
        {
            if (exp[v] != 0)
            {
                status = gr_pow_ui(&t, imgs[v], exp[v], ctx);
                status |= _gr_tower_lazy_mul(&term, &term, &t, ctx);
            }
        }
        status |= _gr_tower_lazy_add(res, res, &term, ctx);
    }

    _gr_tower_lazy_clear(&term, ctx);
    _gr_tower_lazy_clear(&t, ctx);
    fmpz_clear(c);
    flint_free(exp);
    return status;
}

/*
    Complex conjugation. Real generators are fixed, i is sent to -i and
    roots of unity to their inverses; when only such generators occur,
    the conjugate is a substitution within the tower. Otherwise the
    conjugate of each generator is computed as an element of the field:
    conj(exp(u)) = exp(conj(u)), conj(log(u)) = log(conj(u)) and
    conj(root_n(u)) = root_n(conj(u)) away from the branch cut, while
    on the negative real axis log(u) - 2 pi i and root_n(u) e^{-2 pi i/n};
    the element is then evaluated at the images. Other non-real
    algebraic generators are not handled (GR_UNABLE).
*/
/* res = x with the exponent e of each variable v with ord[v] != 0
   replaced by (-e) mod ord[v] (not reduced) */
static void
_flat_conj_exponents(fmpz_mpoly_t res, const fmpz_mpoly_t x, const ulong * ord, const fmpz_mpoly_ctx_t mctx)
{
    slong nvars = mctx->minfo->nvars, t, v;
    ulong * exp = flint_malloc(sizeof(ulong) * FLINT_MAX(nvars, 1));
    fmpz_t c;
    fmpz_mpoly_t r;

    fmpz_init(c);
    fmpz_mpoly_init(r, mctx);
    for (t = 0; t < x->length; t++)
    {
        fmpz_mpoly_get_term_exp_ui(exp, x, t, mctx);
        for (v = 0; v < nvars; v++)
            if (ord[v] != 0)
                exp[v] = (ord[v] - exp[v] % ord[v]) % ord[v];
        fmpz_mpoly_get_term_coeff_fmpz(c, x, t, mctx);
        fmpz_mpoly_push_term_fmpz_ui(r, c, exp, mctx);
    }
    fmpz_mpoly_sort_terms(r, mctx);
    fmpz_mpoly_combine_like_terms(r, mctx);
    fmpz_mpoly_swap(res, r, mctx);
    fmpz_mpoly_clear(r, mctx);
    fmpz_clear(c);
    flint_free(exp);
}

/*
    The conjugate of the algebraic generator d of F which is not an
    absolute algebraic number (a root of a polynomial over transcendental
    generators: the roots of x^5 - pi x - 1): a root of the modulus with
    conjugated coefficients, the one whose enclosure meets the conjugate
    of the generator's enclosure. Root finding may extend the tower; the
    caller detects that and starts over.
*/
static int _gr_tower_lazy_eval_poly(gr_tower_lazy_elem_t res, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_struct * fctx, slong nvars, gr_tower_lazy_elem_struct ** imgs, gr_ctx_t ctx);

static int
_gr_tower_lazy_conj_alg_gen(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, slong d, slong nvars, gr_tower_lazy_elem_struct ** limgs, gr_ctx_t ctx)
{
    gr_tower_struct * T = F->T;
    const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
    slong k = g->index, v = GR_TOWER_FLAT_VAR_D(F, d), i, prec, found = -1;
    const fmpz_mpoly_ctx_struct * mctx = F->mctx;
    fmpz_mpoly_univar_t u;
    fmpz_mpoly_q_t c;
    gr_poly_t P;
    gr_vec_t roots;
    fmpz_vec_t mult;
    gr_tower_lazy_elem_struct a, ca, z;
    acb_t w, e;
    int status = GR_SUCCESS;
    ulong layout = F->layout_version, version = T->version;

    fmpz_mpoly_univar_init(u, mctx);
    fmpz_mpoly_q_init(c, mctx);
    gr_poly_init(P, ctx);
    gr_vec_init(roots, 0, ctx);
    fmpz_vec_init(mult, 0);
    _gr_tower_lazy_init(&a, ctx);
    _gr_tower_lazy_init(&ca, ctx);
    _gr_tower_lazy_init(&z, ctx);
    acb_init(w);
    acb_init(e);

    /* the generator (as an element, before the tower can change) */
    _gr_tower_lazy_set_gen_d(&z, F, d, ctx);

    /* (the univariate form and c stay in this context: F may move to a
       new one while the coefficients are conjugated) */
    fmpz_mpoly_to_univar(u, _gr_tower_flat_ideal_elem(F, k), v, mctx);
    for (i = 0; i < u->length && status == GR_SUCCESS; i++)
    {
        /* through the images of the generators (computed once in the
           caller: conjugating each coefficient on its own would compute
           the conjugates of the generators again, in new towers), or
           else recursively */
        int have = 1, * uv = flint_malloc(sizeof(int) * nvars);
        slong w;
        fmpz_mpoly_used_vars(uv, u->coeffs + i, mctx);
        for (w = 0; w < nvars; w++)
            if (uv[w] && (limgs == NULL || limgs[w] == NULL))
                have = 0;
        flint_free(uv);
        if (have)
            status = _gr_tower_lazy_eval_poly(&ca, u->coeffs + i, mctx, nvars, limgs, ctx);
        else
        {
            fmpz_mpoly_set(fmpz_mpoly_q_numref(c), u->coeffs + i, mctx);
            fmpz_mpoly_one(fmpz_mpoly_q_denref(c), mctx);
            _gr_tower_lazy_set_flat(&a, F, c, mctx, ctx);
            status = _gr_tower_lazy_conj(&ca, &a, ctx);
        }
        if (status == GR_SUCCESS)
            status = gr_poly_set_coeff_scalar(P, fmpz_get_si(u->exps + i), &ca, ctx);
    }

    /* conjugating the coefficients extended the tower: the caller starts
       over (it detects the change after a successful image), when the
       conjugates of the coefficients are already in the tower */
    if (status == GR_SUCCESS && (F->layout_version != layout || F->mctx != mctx || T->version != version))
    {
        status = gr_zero(res, ctx);
        goto cleanup;
    }

    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_poly_roots(roots, mult, P, 0, ctx);

    /* the unique root whose enclosure meets it (the roots are distinct:
       the modulus is irreducible, hence squarefree) */
    for (prec = GR_TOWER_DEFAULT_PREC; status == GR_SUCCESS && found < 0 && prec <= 16 * GR_TOWER_DEFAULT_PREC; prec *= 2)
    {
        slong hits = 0, last = -1;
        if (_gr_tower_lazy_get_acb_impl(w, &z, prec, ctx) != GR_SUCCESS)
            break;
        acb_conj(w, w);
        for (i = 0; i < roots->length; i++)
        {
            if (_gr_tower_lazy_get_acb_impl(e, gr_vec_entry_ptr(roots, i, ctx), prec, ctx) != GR_SUCCESS)
            {
                hits = 2;
                break;
            }
            if (acb_overlaps(e, w))
            {
                hits++;
                last = i;
            }
        }
        if (hits == 1)
            found = last;
        else if (hits == 0)
            break;
    }

    if (status == GR_SUCCESS)
    {
        if (found >= 0)
            status = gr_set(res, gr_vec_entry_ptr(roots, found, ctx), ctx);
        else
            status = GR_UNABLE;
    }

cleanup:
    fmpz_mpoly_univar_clear(u, mctx);
    fmpz_mpoly_q_clear(c, mctx);
    gr_poly_clear(P, ctx);
    gr_vec_clear(roots, ctx);
    fmpz_vec_clear(mult);
    _gr_tower_lazy_clear(&a, ctx);
    _gr_tower_lazy_clear(&ca, ctx);
    _gr_tower_lazy_clear(&z, ctx);
    acb_clear(w);
    acb_clear(e);
    return status;
}

static int
_gr_tower_lazy_conj_attempt(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, int * layout_changed, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    fmpz_mpoly_q_struct ** imgs;
    gr_tower_lazy_elem_struct ** limgs = NULL;
    fmpz_mpoly_q_t r;
    int * used;
    slong d, nvars;
    int simple = 1;
    int status = GR_SUCCESS;

    *layout_changed = 0;
    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;
    gr_tower_flat_ensure(F);

    nvars = F->cap;
    used = flint_malloc(sizeof(int) * nvars);
    fmpz_mpoly_q_used_vars(used, &x->elem.flat.data, F->mctx);
    imgs = flint_calloc(nvars, sizeof(fmpz_mpoly_q_struct *));

    /* substitutions within the tower */
    for (d = 0; d < T->num_gens; d++)
    {
        slong v = GR_TOWER_FLAT_VAR_D(F, d);
        const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);

        if (!used[v])
            continue;

        imgs[v] = flint_malloc(sizeof(fmpz_mpoly_q_struct));
        fmpz_mpoly_q_init(imgs[v], F->mctx);
        fmpz_mpoly_q_gen(imgs[v], v, F->mctx);

        if (_gr_tower_lazy_gen_is_i(T, d))
            fmpz_mpoly_q_neg(imgs[v], imgs[v], F->mctx);
        else if (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_ROOT_OF_UNITY)
        {
            /* the conjugate of a root of unity is its inverse */
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(imgs[v]), fmpz_mpoly_q_numref(imgs[v]), g->def_param - 1, F->mctx);
            status |= gr_tower_flat_reduce(imgs[v], F);
        }
        else if (!_gen_is_real(T, d))
            simple = 0;
    }

    if (status == GR_SUCCESS && simple)
    {
        /* (the only nonreal generators are roots of unity: their inverses
           as exponent maps, e -> -e mod the order, reduced once -- rather
           than a composition with the dense polynomials zeta^(n-1)) */
        ulong * ord = flint_calloc(nvars, sizeof(ulong));
        for (d = 0; d < T->num_gens; d++)
        {
            slong v = GR_TOWER_FLAT_VAR_D(F, d);
            const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
            if (!used[v])
                continue;
            if (_gr_tower_lazy_gen_is_i(T, d))
                ord[v] = 4;
            else if (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_ROOT_OF_UNITY)
                ord[v] = g->def_param;
        }
        fmpz_mpoly_q_init(r, F->mctx);
        _flat_conj_exponents(fmpz_mpoly_q_numref(r), fmpz_mpoly_q_numref(&x->elem.flat.data), ord, F->mctx);
        _flat_conj_exponents(fmpz_mpoly_q_denref(r), fmpz_mpoly_q_denref(&x->elem.flat.data), ord, F->mctx);
        fmpz_mpoly_q_canonicalise(r, F->mctx);
        status = gr_tower_flat_reduce(r, F);
        if (status == GR_SUCCESS)
            _gr_tower_lazy_install(res, F, r, ctx);
        fmpz_mpoly_q_clear(r, F->mctx);
        flint_free(ord);
    }
    else if (status == GR_SUCCESS)
    {
        /* images as elements of the field */
        gr_tower_lazy_elem_struct num, den, u, cu;
        ulong layout = F->layout_version;
        ulong version = T->version;
        const fmpz_mpoly_ctx_struct * mctx0 = F->mctx;   /* (the capacity may grow without a new layout version) */

        limgs = flint_calloc(nvars, sizeof(gr_tower_lazy_elem_struct *));

        /* the generators in the moduli of algebraic generators over
           transcendental ones get images too (in definition order, before
           the generators whose conjugates need them) */
        {
            int * tdep = flint_calloc(T->num_gens, sizeof(int));
            slong * var2d = flint_malloc(sizeof(slong) * nvars);
            int * uv = flint_malloc(sizeof(int) * nvars);
            slong w;
            for (w = 0; w < nvars; w++)
                var2d[w] = -1;
            for (d = 0; d < T->num_gens; d++)
                var2d[GR_TOWER_FLAT_VAR_D(F, d)] = d;
            for (d = 0; d < T->num_gens; d++)
            {
                const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
                if (g->kind != GR_TOWER_ALGEBRAIC)
                {
                    tdep[d] = 1;
                    continue;
                }
                fmpz_mpoly_used_vars(uv, _gr_tower_flat_ideal_elem(F, g->index), F->mctx);
                for (w = 0; w < nvars; w++)
                    if (uv[w] && var2d[w] >= 0 && var2d[w] != d && tdep[var2d[w]])
                        tdep[d] = 1;
            }
            for (d = T->num_gens - 1; d >= 0; d--)
            {
                const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
                if (!used[GR_TOWER_FLAT_VAR_D(F, d)] || g->kind != GR_TOWER_ALGEBRAIC || !tdep[d] ||
                    g->def_kind == GR_TOWER_ROOT_OF_UNITY || _gen_is_real(T, d))
                    continue;
                fmpz_mpoly_used_vars(uv, _gr_tower_flat_ideal_elem(F, g->index), F->mctx);
                for (w = 0; w < nvars; w++)
                    if (uv[w] && var2d[w] >= 0)
                        used[w] = 1;
            }
            flint_free(tdep);
            flint_free(var2d);
            flint_free(uv);

            /* (substitution images of the added generators, as above) */
            for (d = 0; d < T->num_gens; d++)
            {
                slong v = GR_TOWER_FLAT_VAR_D(F, d);
                const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
                if (!used[v] || imgs[v] != NULL)
                    continue;
                imgs[v] = flint_malloc(sizeof(fmpz_mpoly_q_struct));
                fmpz_mpoly_q_init(imgs[v], F->mctx);
                fmpz_mpoly_q_gen(imgs[v], v, F->mctx);
                if (_gr_tower_lazy_gen_is_i(T, d))
                    fmpz_mpoly_q_neg(imgs[v], imgs[v], F->mctx);
                else if (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_ROOT_OF_UNITY)
                {
                    fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(imgs[v]), fmpz_mpoly_q_numref(imgs[v]), g->def_param - 1, F->mctx);
                    status |= gr_tower_flat_reduce(imgs[v], F);
                }
            }
        }

        _gr_tower_lazy_init(&num, ctx);
        _gr_tower_lazy_init(&den, ctx);
        _gr_tower_lazy_init(&u, ctx);
        _gr_tower_lazy_init(&cu, ctx);

        for (d = 0; d < T->num_gens && status == GR_SUCCESS; d++)
        {
            slong v = GR_TOWER_FLAT_VAR_D(F, d);
            const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);

            if (!used[v])
                continue;

            limgs[v] = flint_malloc(sizeof(gr_tower_lazy_elem_struct));
            _gr_tower_lazy_init(limgs[v], ctx);

            if (_gr_tower_lazy_gen_is_i(T, d) || _gen_is_real(T, d) ||
                (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_ROOT_OF_UNITY))
            {
                _gr_tower_lazy_set_flat(limgs[v], F, imgs[v], F->mctx, ctx);
            }
            else if (g->kind == GR_TOWER_PI)
            {
                _gr_tower_lazy_set_gen_d(limgs[v], F, d, ctx);
            }
            else if ((g->kind == GR_TOWER_EXP || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_EXP)) && g->arg.mctx != NULL)
            {
                /* (also an exponential which became algebraic, a radical
                   over other generators: conjugated by its definition) */
                _gr_tower_lazy_set_flat(&u, F, &g->arg.data, g->arg.mctx, ctx);
                status = _gr_tower_lazy_conj(&cu, &u, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_exp_log(limgs[v], &cu, GR_TOWER_EXP, ctx);
            }
            else if (g->kind == GR_TOWER_LOG || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_LOG) ||
                     (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_ROOT && g->arg.mctx != NULL))
            {
                int cut = _flat_on_cut(&g->arg.data, g->arg.mctx, F);
                int is_log = (g->kind == GR_TOWER_LOG || g->def_kind == GR_TOWER_LOG);

                if (cut == 1)
                {
                    /* on the negative real axis */
                    _gr_tower_lazy_set_gen_d(limgs[v], F, d, ctx);
                    if (is_log)
                    {
                        status = _gr_tower_lazy_pi(&u, ctx);
                        status |= _gr_tower_lazy_i(&cu, ctx);
                        status |= _gr_tower_lazy_mul(&u, &u, &cu, ctx);
                        status |= gr_mul_ui(&u, &u, 2, ctx);
                        status |= _gr_tower_lazy_sub(limgs[v], limgs[v], &u, ctx);
                    }
                    else
                    {
                        status = _gr_tower_lazy_root_of_unity(&u, -1, g->def_param, ctx);
                        status |= _gr_tower_lazy_mul(limgs[v], limgs[v], &u, ctx);
                    }
                }
                else if (cut == -1)
                {
                    /* (g may be a dangling pointer after the recursive
                       conjugation, which may extend the tower) */
                    ulong def_param = g->def_param;
                    _gr_tower_lazy_set_flat(&u, F, &g->arg.data, g->arg.mctx, ctx);
                    status = _gr_tower_lazy_conj(&cu, &u, ctx);
                    if (status == GR_SUCCESS)
                    {
                        if (is_log)
                            status = _gr_tower_lazy_exp_log(limgs[v], &cu, GR_TOWER_LOG, ctx);
                        else
                            status = _gr_tower_lazy_root_ui(limgs[v], &cu, def_param, ctx);
                    }
                }
                else
                    status = GR_UNABLE;
            }
            else if ((g->kind == GR_TOWER_JACOBI_THETA || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_JACOBI_THETA)) &&
                     g->arg.mctx != NULL && g->num_xargs == 1)
            {
                /* conj(theta_j(z, tau)) = theta_j(conj(z), -conj(tau)) */
                slong param = g->def_param;
                gr_tower_lazy_elem_struct cargs[2];
                _gr_tower_lazy_init(cargs, ctx);
                _gr_tower_lazy_init(cargs + 1, ctx);
                _gr_tower_lazy_set_flat(cargs, F, &g->arg.data, g->arg.mctx, ctx);
                _gr_tower_lazy_set_flat(cargs + 1, F, &g->xargs[0].data, g->xargs[0].mctx, ctx);
                status = _gr_tower_lazy_conj(cargs, cargs, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_conj(cargs + 1, cargs + 1, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_neg(cargs + 1, cargs + 1, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_special_multi(limgs[v], GR_TOWER_JACOBI_THETA, param, cargs, 2, ctx);
                _gr_tower_lazy_clear(cargs, ctx);
                _gr_tower_lazy_clear(cargs + 1, ctx);
            }
            else if ((g->kind == GR_TOWER_HURWITZ_ZETA || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_HURWITZ_ZETA)) &&
                     g->arg.mctx != NULL && g->num_xargs == 1)
            {
                /* conj(zeta(s, a)) = zeta(conj(s), conj(a)) (Re(a) > 0) */
                gr_tower_lazy_elem_struct cargs[2];
                _gr_tower_lazy_init(cargs, ctx);
                _gr_tower_lazy_init(cargs + 1, ctx);
                _gr_tower_lazy_set_flat(cargs, F, &g->arg.data, g->arg.mctx, ctx);
                _gr_tower_lazy_set_flat(cargs + 1, F, &g->xargs[0].data, g->xargs[0].mctx, ctx);
                status = _gr_tower_lazy_conj(cargs, cargs, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_conj(cargs + 1, cargs + 1, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_special_multi(limgs[v], GR_TOWER_HURWITZ_ZETA, 0, cargs, 2, ctx);
                _gr_tower_lazy_clear(cargs, ctx);
                _gr_tower_lazy_clear(cargs + 1, ctx);
            }
            else if ((g->kind == GR_TOWER_HYPGEOM || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_HYPGEOM)) &&
                     g->arg.mctx != NULL)
            {
                /* conj(pFq(a; b; z)) = pFq(conj(a); conj(b); conj(z)) off
                   the branch cut z >= 1 (p = q + 1) */
                slong i, n = _gr_tower_gen_num_args(g);
                slong param = g->def_param;
                gr_tower_lazy_elem_struct * cargs;
                int ok = 1;

                cargs = flint_malloc(sizeof(gr_tower_lazy_elem_struct) * n);
                for (i = 0; i < n; i++)
                    _gr_tower_lazy_init(cargs + i, ctx);
                /* (g may be a dangling pointer after a recursive
                   conjugation: the arguments are read first) */
                for (i = 0; i < n; i++)
                {
                    const gr_tower_flat_elem_struct * a = _gr_tower_gen_arg_ptr(g, i);
                    _gr_tower_lazy_set_flat(cargs + i, F, &a->data, a->mctx, ctx);
                }

                if (GR_TOWER_HYPGEOM_P(param) == GR_TOWER_HYPGEOM_Q(param) + 1)
                {
                    acb_t z;
                    acb_init(z);
                    if (_gr_tower_lazy_get_acb_impl(z, cargs, GR_TOWER_DEFAULT_PREC, ctx) != GR_SUCCESS)
                        ok = 0;
                    else if (arb_contains_zero(acb_imagref(z)))
                    {
                        arb_t t;
                        arb_init(t);
                        arb_sub_ui(t, acb_realref(z), 1, GR_TOWER_DEFAULT_PREC);
                        ok = (_gr_tower_lazy_is_real_exact(cargs, ctx) == T_TRUE) && arb_is_negative(t);
                        arb_clear(t);
                    }
                    acb_clear(z);
                }

                status = ok ? GR_SUCCESS : GR_UNABLE;
                for (i = 0; i < n && status == GR_SUCCESS; i++)
                    status = _gr_tower_lazy_conj(cargs + i, cargs + i, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_special_multi(limgs[v], GR_TOWER_HYPGEOM, param, cargs, n, ctx);

                for (i = 0; i < n; i++)
                    _gr_tower_lazy_clear(cargs + i, ctx);
                flint_free(cargs);
            }
            else if ((GR_TOWER_KIND_IS_SPECIAL(g->kind) || (g->kind == GR_TOWER_ALGEBRAIC && GR_TOWER_KIND_IS_SPECIAL(g->def_kind))) &&
                     (g->arg.mctx != NULL || g->def_kind == GR_TOWER_CONSTANT))
            {
                /* conj(F(u)) = F(conj(u)) (W_k: W_{-k}) away from the
                   branch cuts (on a cut, the value is not real while
                   F(conj(u)) = F(u)); the canonical form of F(conj(u))
                   is computed by the special function itself */
                int kind = g->def_kind;
                slong param = g->def_param;
                int ok = 1;

                if (kind != GR_TOWER_CONSTANT)
                    _gr_tower_lazy_set_flat(&u, F, &g->arg.data, g->arg.mctx, ctx);

                if (kind == GR_TOWER_POLYLOG || kind == GR_TOWER_LAMBERTW ||
                    kind == GR_TOWER_ELLIPTIC_K || kind == GR_TOWER_ELLIPTIC_E)
                {
                    acb_t z;
                    acb_init(z);
                    if (_gr_tower_lazy_get_acb_impl(z, &u, GR_TOWER_DEFAULT_PREC, ctx) != GR_SUCCESS)
                        ok = 0;
                    else if (arb_contains_zero(acb_imagref(z)))
                    {
                        /* a real point: off the cut only if the
                           function is real there */
                        acb_t w;
                        acb_init(w);
                        acb_set(w, z);
                        arb_zero(acb_imagref(w));
                        ok = (_gr_tower_lazy_is_real_exact(&u, ctx) == T_TRUE) &&
                             _gr_tower_special_real_at(kind, param, w, GR_TOWER_DEFAULT_PREC);
                        acb_clear(w);
                    }
                    acb_clear(z);
                }

                if (!ok)
                    status = GR_UNABLE;
                else if (kind == GR_TOWER_CONSTANT)
                    status = _gr_tower_lazy_special(limgs[v], kind, param, NULL, ctx);
                else if (kind == GR_TOWER_MODULAR_LAMBDA)
                {
                    /* conj(lambda(tau)) = lambda(-conj(tau)) */
                    status = _gr_tower_lazy_conj(&cu, &u, ctx);
                    if (status == GR_SUCCESS)
                        status = _gr_tower_lazy_neg(&cu, &cu, ctx);
                    if (status == GR_SUCCESS)
                        status = _gr_tower_lazy_special(limgs[v], kind, param, &cu, ctx);
                }
                else
                {
                    status = _gr_tower_lazy_conj(&cu, &u, ctx);
                    /* (L(s, chi): conj = L(conj(s), conj(chi)), the
                       conjugate character chi_q(k^(-1))) */
                    if (kind == GR_TOWER_DIRICHLET_L)
                        param = GR_TOWER_DIRICHLET_PARAM(GR_TOWER_DIRICHLET_Q(param),
                            n_invmod(GR_TOWER_DIRICHLET_K(param), GR_TOWER_DIRICHLET_Q(param)));
                    if (status == GR_SUCCESS)
                        status = _gr_tower_lazy_special(limgs[v], kind, (kind == GR_TOWER_LAMBERTW) ? -param : param, &cu, ctx);
                }
            }
            else if (g->kind == GR_TOWER_ALGEBRAIC)
            {
                /* an algebraic number: conjugate as a qqbar; the conjugate
                   is a root of the same polynomial and is absorbed into F
                   next to the generator */
                qqbar_t q;
                qqbar_init(q);
                _gr_tower_lazy_set_gen_d(&u, F, d, ctx);
                status = _gr_tower_lazy_get_qqbar(q, &u, ctx);
                if (status == GR_SUCCESS)
                {
                    qqbar_conj(q, q);
                    status = _gr_tower_lazy_set_qqbar(limgs[v], q, ctx);
                }
                else
                {
                    /* over transcendental generators: a root of the
                       conjugated modulus */
                    status = _gr_tower_lazy_conj_alg_gen(limgs[v], F, d, nvars, limgs, ctx);
                }
                qqbar_clear(q);
            }
            else
                status = GR_UNABLE;

            /* the image computations may extend or refine the tower: the
               variables (and the arrays indexed by them) are then stale;
               start over */
            if (status == GR_SUCCESS && (F->layout_version != layout || F->mctx != mctx0 || T->version != version))
            {
                *layout_changed = 1;
                status = GR_UNABLE;
            }
        }

        /* the computations above may have changed the layout of F, or
           refined its moduli (so that x, reduced, may involve generators
           it did not involve before, without images): the caller starts
           over in that case */
        if (status == GR_SUCCESS && (F->layout_version != layout || F->mctx != mctx0 || T->version != version))
        {
            *layout_changed = 1;
            status = GR_UNABLE;
        }

        if (status == GR_SUCCESS)
        {
            x = _gr_tower_lazy_flat_view(x);
            /* every variable of x must have an image */
            fmpz_mpoly_q_used_vars(used, &x->elem.flat.data, F->mctx);
            for (d = 0; d < nvars; d++)
                if (used[d] && limgs[d] == NULL)
                    status = GR_UNABLE;
        }

        if (status == GR_SUCCESS)
        {
            /* (a copy of the data: x itself may be updated by the
               operations below, which would move it to another context) */
            fmpz_mpoly_q_t xd;
            const fmpz_mpoly_ctx_struct * xctx = x->elem.flat.mctx;
            fmpz_mpoly_q_init(xd, xctx);
            fmpz_mpoly_q_set(xd, &x->elem.flat.data, xctx);
            status = _gr_tower_lazy_eval_poly(&num, fmpz_mpoly_q_numref(xd), xctx, nvars, limgs, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_eval_poly(&den, fmpz_mpoly_q_denref(xd), xctx, nvars, limgs, ctx);
            fmpz_mpoly_q_clear(xd, xctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_div(res, &num, &den, ctx);
        }

        for (d = 0; d < nvars; d++)
        {
            if (limgs[d] != NULL)
            {
                _gr_tower_lazy_clear(limgs[d], ctx);
                flint_free(limgs[d]);
            }
        }
        flint_free(limgs);
        _gr_tower_lazy_clear(&num, ctx);
        _gr_tower_lazy_clear(&den, ctx);
        _gr_tower_lazy_clear(&u, ctx);
        _gr_tower_lazy_clear(&cu, ctx);
    }

    for (d = 0; d < nvars; d++)
    {
        if (imgs[d] != NULL)
        {
            fmpz_mpoly_q_clear(imgs[d], F->mctx);
            flint_free(imgs[d]);
        }
    }
    flint_free(imgs);
    flint_free(used);

    return status;
}

int
_gr_tower_lazy_conj(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int layout_changed, iter;
    int status = GR_UNABLE;

    /* conjugating the generators may adjoin new generators to the field
       (which changes the layout of the flat representation): retry with
       the extended field, where the images are already present */
    LAZY(ctx)->conj_depth++;
    for (iter = 0; iter < 8; iter++)
    {
        status = _gr_tower_lazy_conj_attempt(res, x, &layout_changed, ctx);
        if (!(status == GR_UNABLE && layout_changed))
            break;
    }
    LAZY(ctx)->conj_depth--;

    return status;
}

/* -------------------------------------------------------------------- */
/* transfer between lazy contexts                                        */
/* -------------------------------------------------------------------- */

/* the data of one generator of the source, read under the source lock */
typedef struct
{
    int kind;                   /* def_kind, or GR_TOWER_ALGEBRAIC for a root given by its modulus */
    slong param;
    int has_arg;
    fmpz_mpoly_q_struct arg;    /* in the snapshot context */
    slong num_xargs;            /* functions of several arguments: the others (in the snapshot context) */
    fmpz_mpoly_q_struct * xargs;
    int absolute;               /* an absolute algebraic number: q */
    qqbar_struct q;
    fmpz_mpoly_struct modulus;  /* a root of this polynomial (in the snapshot context, in variable var) */
    int has_modulus;
    acb_struct z;               /* its enclosure (256 bits) */
    slong var;
}

transfer_gen_struct;

/*
    res (in ctx) = x (in x_ctx, another lazy context): the generators x
    involves (with their definitions, in definition order) are recreated
    in ctx from their definitions, and x is evaluated at the images.
    Every definition kind is transferred exactly (the same function of
    the image of the argument, with the same branch conventions, the same
    algebraic number); a root given only by its modulus over
    transcendental generators is identified among the roots of the mapped
    modulus by its enclosure (a unique overlap at increasing precision).
    The source is read under its lock, then released before anything is
    built in ctx (ctx's lock is held by the caller). Two threads
    transferring in opposite directions at once would be a deadlock; the
    contexts are therefore locked in address order and the caller's lock
    on ctx is taken recursively.
*/
int
_gr_tower_lazy_transfer(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    fmpz_mpoly_ctx_struct * mctx;
    fmpz_mpoly_q_t data;
    transfer_gen_struct * gens = NULL;
    gr_tower_lazy_elem_struct ** imgs = NULL;
    slong * var_of = NULL;
    int * need = NULL, * clos = NULL, * used = NULL;
    slong G, nvars, d, e, v;
    int status = GR_SUCCESS;

    /* ---- the source, under its lock ---- */
    _gr_tower_lazy_lock(x_ctx);
    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;
    gr_tower_flat_ensure(F);
    mctx = F->mctx;
    G = T->num_gens;
    nvars = F->cap;

    fmpz_mpoly_q_init(data, mctx);
    if (x->elem.flat.mctx == mctx)
        fmpz_mpoly_q_set(data, &x->elem.flat.data, mctx);
    else
        gr_tower_flat_convert(data, &x->elem.flat.data, x->elem.flat.mctx, F);

    used = flint_calloc(nvars, sizeof(int));
    need = flint_calloc(FLINT_MAX(G, 1), sizeof(int));
    clos = flint_calloc(FLINT_MAX(G * G, 1), sizeof(int));
    var_of = flint_malloc(sizeof(slong) * FLINT_MAX(G, 1));
    fmpz_mpoly_q_used_vars(used, data, mctx);
    for (d = 0; d < G; d++)
    {
        var_of[d] = GR_TOWER_FLAT_VAR_D(F, d);
        clos[d * G + d] = 1;
        _gr_tower_involved_gens_closure(clos + d * G, T);
    }
    for (d = 0; d < G; d++)
        if (used[var_of[d]])
            for (e = 0; e < G; e++)
                if (clos[d * G + e])
                    need[e] = 1;

    gens = flint_calloc(FLINT_MAX(G, 1), sizeof(transfer_gen_struct));
    for (d = 0; d < G && status == GR_SUCCESS; d++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
        transfer_gen_struct * r = gens + d;

        if (!need[d])
            continue;
        r->kind = g->def_kind;
        r->param = g->def_param;
        r->var = var_of[d];
        if (g->arg.mctx != NULL && GR_TOWER_KIND_HAS_ARG(g->def_kind))
        {
            fmpz_mpoly_q_init(&r->arg, mctx);
            gr_tower_flat_convert(&r->arg, &g->arg.data, g->arg.mctx, F);
            r->has_arg = 1;
            if (g->num_xargs > 0)
            {
                slong i;
                r->num_xargs = g->num_xargs;
                r->xargs = flint_malloc(sizeof(fmpz_mpoly_q_struct) * g->num_xargs);
                for (i = 0; i < g->num_xargs; i++)
                {
                    fmpz_mpoly_q_init(r->xargs + i, mctx);
                    gr_tower_flat_convert(r->xargs + i, &g->xargs[i].data, g->xargs[i].mctx, F);
                }
            }
        }
        else if (g->def_kind == GR_TOWER_ROOT && g->arg.mctx != NULL)
        {
            fmpz_mpoly_q_init(&r->arg, mctx);
            gr_tower_flat_convert(&r->arg, &g->arg.data, g->arg.mctx, F);
            r->has_arg = 1;
        }
        else if (g->def_kind == GR_TOWER_ALGEBRAIC && g->kind == GR_TOWER_ALGEBRAIC)
        {
            /* a root: an absolute algebraic number, or a root of its
               modulus over transcendental generators */
            int absolute = 1;
            for (e = 0; e < G; e++)
                if (clos[d * G + e] && GR_TOWER_GEN(T, e)->kind != GR_TOWER_ALGEBRAIC)
                    absolute = 0;
            if (absolute)
            {
                gr_tower_lazy_elem_struct u;
                _gr_tower_lazy_init(&u, x_ctx);
                _gr_tower_lazy_set_gen_d(&u, F, d, x_ctx);
                qqbar_init(&r->q);
                status = _gr_tower_lazy_get_qqbar(&r->q, &u, x_ctx);
                r->absolute = 1;
                _gr_tower_lazy_clear(&u, x_ctx);
                /* (the tower may have moved to a new context) */
                gr_tower_flat_ensure(F);
                if (F->mctx != mctx)
                    status = GR_UNABLE;
            }
            else
            {
                fmpz_mpoly_init(&r->modulus, mctx);
                fmpz_mpoly_set(&r->modulus, _gr_tower_flat_ideal_elem(F, g->index), mctx);
                r->has_modulus = 1;
                acb_init(&r->z);
                status = gr_tower_step_get_acb(&r->z, T, g->index, 256);
            }
        }
        else if (g->def_kind == GR_TOWER_FREE)
            status = GR_UNABLE;
    }
    _gr_tower_lazy_unlock(x_ctx);

    /* ---- the images, in ctx ---- */
    if (status == GR_SUCCESS)
    {
        gr_tower_lazy_elem_struct u, num, den;
        _gr_tower_lazy_init(&u, ctx);
        _gr_tower_lazy_init(&num, ctx);
        _gr_tower_lazy_init(&den, ctx);
        imgs = flint_calloc(nvars, sizeof(gr_tower_lazy_elem_struct *));
        for (d = 0; d < G && status == GR_SUCCESS; d++)
        {
            transfer_gen_struct * r = gens + d;
            gr_tower_lazy_elem_struct * img;

            if (!need[d])
                continue;
            img = flint_malloc(sizeof(gr_tower_lazy_elem_struct));
            _gr_tower_lazy_init(img, ctx);
            imgs[r->var] = img;

            if (r->has_arg)
            {
                status = _gr_tower_lazy_eval_poly(&num, fmpz_mpoly_q_numref(&r->arg), mctx, nvars, imgs, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_eval_poly(&den, fmpz_mpoly_q_denref(&r->arg), mctx, nvars, imgs, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_div(&u, &num, &den, ctx);
                if (status != GR_SUCCESS)
                    break;
            }

            switch (r->kind)
            {
                case GR_TOWER_PI:
                    status = _gr_tower_lazy_pi(img, ctx);
                    break;
                case GR_TOWER_ROOT_OF_UNITY:
                    status = _gr_tower_lazy_root_of_unity(img, 1, (ulong) r->param, ctx);
                    break;
                case GR_TOWER_TAN_PI:
                    {
                        qqbar_t z;
                        qqbar_init(z);
                        status = _gr_tower_lazy_tan_pi_gen(img, (ulong) r->param, z, ctx);
                        qqbar_clear(z);
                    }
                    break;
                case GR_TOWER_ROOT:
                    status = r->has_arg ? _gr_tower_lazy_root_ui(img, &u, (ulong) r->param, ctx) : GR_UNABLE;
                    break;
                case GR_TOWER_EXP:
                case GR_TOWER_LOG:
                    status = r->has_arg ? _gr_tower_lazy_exp_log(img, &u, r->kind, ctx) : GR_UNABLE;
                    break;
                case GR_TOWER_TAN:
                case GR_TOWER_ATAN:
                    status = r->has_arg ? _gr_tower_lazy_trig_locked(img, &u, r->kind, ctx) : GR_UNABLE;
                    break;
                case GR_TOWER_CONSTANT:
                    status = _gr_tower_lazy_special(img, r->kind, r->param, NULL, ctx);
                    break;
                case GR_TOWER_ALGEBRAIC:
                    if (r->absolute)
                        status = _gr_tower_lazy_set_qqbar(img, &r->q, ctx);
                    else if (r->has_modulus)
                    {
                        /* the roots of the mapped modulus; the one whose
                           enclosure meets z (unique: the roots are distinct
                           and their enclosures shrink) */
                        fmpz_mpoly_univar_t uni;
                        gr_poly_t P;
                        gr_vec_t roots;
                        fmpz_vec_t mult;
                        slong i, prec, found = -1;
                        fmpz_mpoly_univar_init(uni, mctx);
                        gr_poly_init(P, ctx);
                        gr_vec_init(roots, 0, ctx);
                        fmpz_vec_init(mult, 0);
                        fmpz_mpoly_to_univar(uni, &r->modulus, r->var, mctx);
                        for (i = 0; i < uni->length && status == GR_SUCCESS; i++)
                        {
                            status = _gr_tower_lazy_eval_poly(&num, uni->coeffs + i, mctx, nvars, imgs, ctx);
                            if (status == GR_SUCCESS)
                                status = gr_poly_set_coeff_scalar(P, fmpz_get_si(uni->exps + i), &num, ctx);
                        }
                        if (status == GR_SUCCESS)
                            status = _gr_tower_lazy_poly_roots(roots, mult, P, 0, ctx);
                        for (prec = GR_TOWER_DEFAULT_PREC; status == GR_SUCCESS && found < 0 && prec <= 256; prec *= 2)
                        {
                            acb_t w;
                            slong hits = 0, last = -1;
                            acb_init(w);
                            for (i = 0; i < roots->length; i++)
                            {
                                if (_gr_tower_lazy_get_acb_impl(w, gr_vec_entry_ptr(roots, i, ctx), prec, ctx) != GR_SUCCESS)
                                {
                                    hits = -1;
                                    break;
                                }
                                if (acb_overlaps(w, &r->z))
                                    hits++, last = i;
                            }
                            acb_clear(w);
                            if (hits == 1)
                                found = last;
                            else if (hits <= 0)
                                break;
                        }
                        if (status == GR_SUCCESS)
                            status = (found >= 0) ? gr_set(img, gr_vec_entry_ptr(roots, found, ctx), ctx) : GR_UNABLE;
                        fmpz_mpoly_univar_clear(uni, mctx);
                        gr_poly_clear(P, ctx);
                        gr_vec_clear(roots, ctx);
                        fmpz_vec_clear(mult);
                    }
                    else
                        status = GR_UNABLE;
                    break;
                default:
                    if (GR_TOWER_KIND_IS_SPECIAL(r->kind) && r->has_arg && r->num_xargs > 0)
                    {
                        /* the other arguments, evaluated at the images */
                        slong i, n = 1 + r->num_xargs;
                        gr_tower_lazy_elem_struct * xa = flint_malloc(sizeof(gr_tower_lazy_elem_struct) * n);
                        for (i = 0; i < n; i++)
                            _gr_tower_lazy_init(xa + i, ctx);
                        status = _gr_tower_lazy_set(xa, &u, ctx);
                        for (i = 1; i < n && status == GR_SUCCESS; i++)
                        {
                            status = _gr_tower_lazy_eval_poly(&num, fmpz_mpoly_q_numref(r->xargs + i - 1), mctx, nvars, imgs, ctx);
                            if (status == GR_SUCCESS)
                                status = _gr_tower_lazy_eval_poly(&den, fmpz_mpoly_q_denref(r->xargs + i - 1), mctx, nvars, imgs, ctx);
                            if (status == GR_SUCCESS)
                                status = _gr_tower_lazy_div(xa + i, &num, &den, ctx);
                        }
                        if (status == GR_SUCCESS)
                            status = _gr_tower_lazy_special_multi(img, r->kind, r->param, xa, n, ctx);
                        for (i = 0; i < n; i++)
                            _gr_tower_lazy_clear(xa + i, ctx);
                        flint_free(xa);
                    }
                    else if (GR_TOWER_KIND_IS_SPECIAL(r->kind) && r->has_arg)
                        status = _gr_tower_lazy_special(img, r->kind, r->param, &u, ctx);
                    else
                        status = GR_UNABLE;
            }
        }

        if (status == GR_SUCCESS)
        {
            status = _gr_tower_lazy_eval_poly(&num, fmpz_mpoly_q_numref(data), mctx, nvars, imgs, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_eval_poly(&den, fmpz_mpoly_q_denref(data), mctx, nvars, imgs, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_div(res, &num, &den, ctx);
        }

        for (v = 0; v < nvars; v++)
            if (imgs[v] != NULL)
            {
                _gr_tower_lazy_clear(imgs[v], ctx);
                flint_free(imgs[v]);
            }
        flint_free(imgs);
        _gr_tower_lazy_clear(&u, ctx);
        _gr_tower_lazy_clear(&num, ctx);
        _gr_tower_lazy_clear(&den, ctx);
    }

    for (d = 0; d < G; d++)
    {
        transfer_gen_struct * r = gens + d;
        if (r->has_arg)
            fmpz_mpoly_q_clear(&r->arg, mctx);
        if (r->num_xargs > 0)
        {
            slong i;
            for (i = 0; i < r->num_xargs; i++)
                fmpz_mpoly_q_clear(r->xargs + i, mctx);
            flint_free(r->xargs);
        }
        if (r->absolute)
            qqbar_clear(&r->q);
        if (r->has_modulus)
        {
            fmpz_mpoly_clear(&r->modulus, mctx);
            acb_clear(&r->z);
        }
    }
    flint_free(gens);
    fmpz_mpoly_q_clear(data, mctx);
    flint_free(used);
    flint_free(need);
    flint_free(clos);
    flint_free(var_of);
    return status;
}

/* |x| = sqrt(x conj(x)) */
int
_gr_tower_lazy_abs(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct t;
    int status;

    _gr_tower_lazy_init(&t, ctx);
    status = _gr_tower_lazy_conj(&t, x, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_mul(&t, &t, x, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_root_ui(res, &t, 2, ctx);
    _gr_tower_lazy_clear(&t, ctx);
    return status;
}

/*
    x^y: for a rational exponent p/q, the principal q-th root raised to
    the p-th power; otherwise exp(y log(x)).
*/
int
_gr_tower_lazy_pow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y_in, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * y = (gr_tower_lazy_elem_struct *) y_in;
    int status;

    y = _gr_tower_lazy_flat_view(y);

    if (y->F == LAZY(ctx)->trivial)
    {
        fmpq_t q;
        gr_tower_lazy_elem_struct t;

        fmpq_init(q);
        (void) fmpz_mpoly_q_get_fmpq(q, &y->elem.flat.data, y->elem.flat.mctx);

        /* a positive rational x to a negative fractional power: (1/x)^(-q),
           the root of a rational number (rather than the inverse of a
           root, a division in the tower of the radicals, which may be
           large) */
        if (fmpq_sgn(q) < 0 && !fmpz_is_one(fmpq_denref(q)))
        {
            fmpq_t c;
            fmpq_init(c);
            if (_gr_tower_lazy_rational_repr_locked(c, x, ctx) && fmpq_sgn(c) > 0)
            {
                fmpq_inv(c, c);
                fmpq_neg(q, q);
                _gr_tower_lazy_init(&t, ctx);
                status = _gr_tower_lazy_set_fmpq(&t, c, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_root_ui(&t, &t, fmpz_get_ui(fmpq_denref(q)), ctx);
                if (status == GR_SUCCESS)
                    status = gr_pow_fmpz(res, &t, fmpq_numref(q), ctx);
                _gr_tower_lazy_clear(&t, ctx);
                fmpq_clear(c);
                fmpq_clear(q);
                return status;
            }
            fmpq_clear(c);
        }

        if (fmpz_is_one(fmpq_denref(q)))
        {
            status = gr_pow_fmpz(res, x, fmpq_numref(q), ctx);
        }
        else if (fmpz_fits_si(fmpq_denref(q)))
        {
            _gr_tower_lazy_init(&t, ctx);
            status = _gr_tower_lazy_root_ui(&t, x, fmpz_get_ui(fmpq_denref(q)), ctx);
            if (status == GR_SUCCESS)
                status = gr_pow_fmpz(res, &t, fmpq_numref(q), ctx);
            _gr_tower_lazy_clear(&t, ctx);
        }
        else
            status = GR_UNABLE;

        fmpq_clear(q);
        return status;
    }

    {
        gr_tower_lazy_elem_struct t;
        _gr_tower_lazy_init(&t, ctx);
        status = _gr_tower_lazy_exp_log(&t, x, GR_TOWER_LOG, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_mul(&t, &t, y, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_exp_log(res, &t, GR_TOWER_EXP, ctx);
        _gr_tower_lazy_clear(&t, ctx);
        return status;
    }
}

int
_gr_tower_lazy_sqrt(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    return _gr_tower_lazy_root_ui(res, x, 2, ctx);
}

POP_OPTIONS
