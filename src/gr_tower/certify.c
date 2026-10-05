/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Richardson's algorithm: deciding whether an element of a tower with
    exp and log generators is zero, assuming Schanuel's conjecture (the
    conjecture is only needed for termination; every answer is proved).

    The log side of a tower consists of the arguments a of the generators
    exp(a), the generators log(u) themselves, and 2 pi i (when both pi and
    i are present). A Q-linear relation between these numbers (found by
    LLL from numerical values) is verified by a zero test in the tower
    below its highest generator (recursively), and then used to eliminate
    that generator: it becomes algebraic (of degree 1, or a radical of
    higher degree) over the generators preceding it. Since the generators
    are kept, flat elements remain valid; the reduction modulo the new
    modulus performs the substitution.

    If x is nonzero as a rational function of the transcendental
    generators but numerically indistinguishable from zero at the current
    precision, either a relation is found and eliminated, or the precision
    is doubled. Under Schanuel's conjecture, a zero x always leads to a
    relation.
*/

#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "ulong_extras.h"
#include "fmpz_mpoly_q.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpq_vec.h"
#include "arb.h"
#include <stdio.h>
#include "acb.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

int
_gr_tower_has_conjectural(const gr_tower_t T)
{
    slong d;
    for (d = 0; d < T->num_gens; d++)
        if (T->gens[d].kind != GR_TOWER_ALGEBRAIC && T->gens[d].status != GR_TOWER_STATUS_PROVEN)
            return 1;
    return 0;
}

/* marks (by definition order) the generators x involves, with those
   their definitions involve, recursively */
void
_gr_tower_flat_involved(int * mark, const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong d;

    for (d = 0; d < T->num_gens; d++)
    {
        slong v = GR_TOWER_FLAT_VAR_D(F, d);
        mark[d] = (fmpz_mpoly_degree_si(fmpz_mpoly_q_numref(x), v, F->mctx) > 0 ||
                   fmpz_mpoly_degree_si(fmpz_mpoly_q_denref(x), v, F->mctx) > 0);
    }
    _gr_tower_involved_gens_closure(mark, T);
}

/* whether x involves (through the definitions, recursively) a
   transcendental generator which is not proven transcendental: a
   nonzero element of the field is otherwise a nonzero number, whatever
   the other generators of the tower */
int
_gr_tower_flat_involves_conjectural(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    int * mark;
    slong d;
    int res = 0;

    if (!_gr_tower_has_conjectural(T))
        return 0;

    mark = flint_calloc(FLINT_MAX(T->num_gens, 1), sizeof(int));
    _gr_tower_flat_involved(mark, x, F);
    for (d = 0; d < T->num_gens && !res; d++)
        if (mark[d] && T->gens[d].kind != GR_TOWER_ALGEBRAIC && T->gens[d].status != GR_TOWER_STATUS_PROVEN)
            res = 1;
    flint_free(mark);
    return res;
}

/* -------------------------------------------------------------------- */
/* the log side                                                          */
/* -------------------------------------------------------------------- */

typedef struct
{
    fmpz_mpoly_q_struct val;     /* the log-side number, as a flat element */
    slong level;                 /* definition order of the generator it belongs to */
    slong gen;                   /* that generator, or -1 for 2 pi i */
}
logside_entry_struct;

/* Definition order of a generator equal to i (root of x^2 + 1 with
   positive imaginary part), or -1. */
static slong
_find_i(gr_tower_t T)
{
    slong k;

    for (k = 1; k <= T->length; k++)
    {
        const gr_poly_struct * m = gr_tower_step_minpoly(T, k);
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        const gr_tower_gen_struct * g = GR_TOWER_STEP(T, k - 1);

        /* a root of unity of order divisible by 4 has i as a power */
        if (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param % 4 == 0)
            return g->def_order;

        if (m->length == 3 && gr_is_one(gr_poly_coeff_srcptr(m, 0, below), below) == T_TRUE
                           && gr_is_zero(gr_poly_coeff_srcptr(m, 1, below), below) == T_TRUE
                           && arb_is_positive(acb_imagref(&g->enclosure)))
            return g->def_order;
    }

    return -1;
}

/* The flat element i, given the generator found by _find_i. */
static void
_flat_i(fmpz_mpoly_q_t res, gr_tower_flat_t F, slong di)
{
    const gr_tower_gen_struct * g = GR_TOWER_GEN(F->T, di);

    fmpz_mpoly_q_gen(res, GR_TOWER_FLAT_VAR_D(F, di), F->mctx);
    if (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param != 4)
    {
        fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), g->def_param / 4, F->mctx);
        GR_MUST_SUCCEED(gr_tower_flat_reduce(res, F));
    }
}

static slong
_find_pi(gr_tower_t T, slong limit)
{
    slong d;
    for (d = 0; d < limit; d++)
        if (T->gens[d].kind == GR_TOWER_PI)
            return d;
    return -1;
}

/* Builds the log side of the generators with definition order < limit. */
static slong
_logside(logside_entry_struct ** L_out, gr_tower_flat_t F, slong limit, const int * mask)
{
    gr_tower_struct * T = F->T;
    logside_entry_struct * L;
    slong d, n = 0, di, dpi;

    gr_tower_flat_ensure(F);

    L = flint_malloc(sizeof(logside_entry_struct) * (2 * limit + 1));

    for (d = 0; d < limit; d++)
    {
        const gr_tower_gen_struct * g = T->gens + d;

        if (mask != NULL && !mask[d])
            continue;

        if (g->kind == GR_TOWER_EXP)
        {
            fmpz_mpoly_q_init(&L[n].val, F->mctx);
            gr_tower_flat_convert(&L[n].val, &g->arg.data, g->arg.mctx, F);
            L[n].level = d;
            L[n].gen = d;
            n++;
        }
        else if (g->kind == GR_TOWER_LOG || g->kind == GR_TOWER_LAMBERTW)
        {
            /* (Lambert W: w = W(z) with e^w = z / w) */
            fmpz_mpoly_q_init(&L[n].val, F->mctx);
            fmpz_mpoly_q_gen(&L[n].val, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            L[n].level = d;
            L[n].gen = d;
            n++;
        }
    }

    di = _find_i(T);
    dpi = _find_pi(T, limit);

    /* the angles of the trigonometric generators, as i A with
       A = 2 u for tan(u) and A = 2 atan(u) (usable when i is present:
       e^{i A} = (1 + i q) / (1 - i q) with q = tan(u) or u) */
    if (di >= 0 && di < limit)
    {
        for (d = 0; d < limit; d++)
        {
            const gr_tower_gen_struct * g = T->gens + d;
            fmpz_mpoly_q_t ii;

            if (g->kind != GR_TOWER_TAN && g->kind != GR_TOWER_ATAN)
                continue;

            fmpz_mpoly_q_init(&L[n].val, F->mctx);
            fmpz_mpoly_q_init(ii, F->mctx);
            if (g->kind == GR_TOWER_TAN)
                gr_tower_flat_convert(&L[n].val, &g->arg.data, g->arg.mctx, F);
            else
                fmpz_mpoly_q_gen(&L[n].val, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            _flat_i(ii, F, di);
            fmpz_mpoly_q_mul(&L[n].val, &L[n].val, ii, F->mctx);
            fmpz_mpoly_q_mul_si(&L[n].val, &L[n].val, 2, F->mctx);
            fmpz_mpoly_q_clear(ii, F->mctx);
            L[n].level = FLINT_MAX(d, di);
            L[n].gen = d;
            n++;
        }
    }

    if (di >= 0 && di < limit && dpi >= 0)
    {
        fmpz_mpoly_q_t t;
        fmpz_mpoly_q_init(&L[n].val, F->mctx);
        fmpz_mpoly_q_init(t, F->mctx);
        _flat_i(&L[n].val, F, di);
        fmpz_mpoly_q_gen(t, GR_TOWER_FLAT_VAR_D(F, dpi), F->mctx);
        fmpz_mpoly_q_mul(&L[n].val, &L[n].val, t, F->mctx);
        fmpz_mpoly_q_mul_si(&L[n].val, &L[n].val, 2, F->mctx);
        fmpz_mpoly_q_clear(t, F->mctx);
        L[n].level = FLINT_MAX(di, dpi);
        L[n].gen = -1;
        n++;
    }

    *L_out = L;
    return n;
}

static void
_logside_clear(logside_entry_struct * L, slong n, gr_tower_flat_t F)
{
    slong i;
    for (i = 0; i < n; i++)
        fmpz_mpoly_q_clear(&L[i].val, F->mctx);
    flint_free(L);
}

/* e^l for a log-side entry: the generator exp(a), the argument of
   log(u), or 1 for 2 pi i. */
/* q = tan(u) (the generator) for tan(u), q = u for atan(u): the angle
   A of the generator has cos A = (1 - q^2) / (1 + q^2),
   sin A = 2 q / (1 + q^2), e^{i A} = (1 + i q) / (1 - i q) */
static void
_trig_q(fmpz_mpoly_q_t q, slong d, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    if (T->gens[d].kind == GR_TOWER_TAN)
        fmpz_mpoly_q_gen(q, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
    else
        gr_tower_flat_convert(q, &T->gens[d].arg.data, T->gens[d].arg.mctx, F);
}

static void
_expside(fmpz_mpoly_q_t res, const logside_entry_struct * e, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;

    if (e->gen < 0)
        fmpz_mpoly_q_one(res, F->mctx);
    else if (T->gens[e->gen].kind == GR_TOWER_EXP)
        fmpz_mpoly_q_gen(res, GR_TOWER_FLAT_VAR_D(F, e->gen), F->mctx);
    else if (T->gens[e->gen].kind == GR_TOWER_TAN || T->gens[e->gen].kind == GR_TOWER_ATAN)
    {
        /* (1 + i q) / (1 - i q) */
        fmpz_mpoly_q_t q, ii, a, b;
        fmpz_mpoly_q_init(q, F->mctx);
        fmpz_mpoly_q_init(ii, F->mctx);
        fmpz_mpoly_q_init(a, F->mctx);
        fmpz_mpoly_q_init(b, F->mctx);
        _trig_q(q, e->gen, F);
        _flat_i(ii, F, _find_i(T));
        fmpz_mpoly_q_mul(q, q, ii, F->mctx);
        fmpz_mpoly_q_add_si(a, q, 1, F->mctx);
        fmpz_mpoly_q_neg(b, q, F->mctx);
        fmpz_mpoly_q_add_si(b, b, 1, F->mctx);
        fmpz_mpoly_q_div(res, a, b, F->mctx);
        (void) gr_tower_flat_rationalize(res, F);
        fmpz_mpoly_q_clear(q, F->mctx);
        fmpz_mpoly_q_clear(ii, F->mctx);
        fmpz_mpoly_q_clear(a, F->mctx);
        fmpz_mpoly_q_clear(b, F->mctx);
    }
    else if (T->gens[e->gen].kind == GR_TOWER_LAMBERTW)
    {
        /* e^w = z / w */
        fmpz_mpoly_q_t w;
        fmpz_mpoly_q_init(w, F->mctx);
        gr_tower_flat_convert(res, &T->gens[e->gen].arg.data, T->gens[e->gen].arg.mctx, F);
        fmpz_mpoly_q_gen(w, GR_TOWER_FLAT_VAR_D(F, e->gen), F->mctx);
        fmpz_mpoly_q_div(res, res, w, F->mctx);
        fmpz_mpoly_q_clear(w, F->mctx);
    }
    else
        gr_tower_flat_convert(res, &T->gens[e->gen].arg.data, T->gens[e->gen].arg.mctx, F);
}

/* res = res * x^e (e may be negative) */
static void
_mul_pow_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong e, gr_tower_flat_t F)
{
    fmpz_mpoly_q_t t;
    slong i;

    if (e == 0)
        return;

    fmpz_mpoly_q_init(t, F->mctx);
    if (e < 0)
    {
        fmpz_mpoly_q_inv(t, x, F->mctx);
        /* an algebraic denominator (the inverse of a generator, say)
           is cleared through the moduli when possible: substituting the
           moduli into a denominator would blow it up */
        (void) gr_tower_flat_rationalize(t, F);
        e = -e;
    }
    else
        fmpz_mpoly_q_set(t, x, F->mctx);

    for (i = 0; i < e; i++)
        fmpz_mpoly_q_mul(res, res, t, F->mctx);

    fmpz_mpoly_q_clear(t, F->mctx);
}

/* -------------------------------------------------------------------- */
/* integer relations                                                     */
/* -------------------------------------------------------------------- */

/*
    All relations found by acb_lindep (the rows of the LLL-reduced lattice
    which give numerically zero combinations of vec), as the rows of the
    matrix rel, which is initialized by this function with that number of
    rows. Returns the number of rows.
*/
static slong
_lindep_all_core(fmpz_mat_t rel, acb_srcptr vec, slong n, slong prec)
{
    fmpz_mat_t A;
    slong i, j, num;

    fmpz_mat_init(A, n, n);
    num = acb_lindep(A, vec, n, prec);
    fmpz_mat_init(rel, num, n);
    for (i = 0; i < num; i++)
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(rel, i, j), fmpz_mat_entry(A, i, j));
    fmpz_mat_clear(A);
    return num;
}

/*
    Entries which could not be evaluated, or only to a poor accuracy,
    would abort or degrade the search: the lattice is built over the
    others, the relations found being extended by zero coefficients.
*/
static slong
_lindep_all(fmpz_mat_t rel, acb_srcptr vec, slong n, slong prec)
{
    slong i, j, nn = 0, num;
    slong * idx;
    acb_ptr sub;
    fmpz_mat_t subrel;

    idx = flint_malloc(sizeof(slong) * FLINT_MAX(n, 1));
    for (i = 0; i < n; i++)
        if (acb_is_finite(vec + i) && (acb_is_zero(vec + i) || acb_rel_accuracy_bits(vec + i) >= 16))
            idx[nn++] = i;

    if (nn == n)
    {
        flint_free(idx);
        return _lindep_all_core(rel, vec, n, prec);
    }

    if (nn == 0)
    {
        flint_free(idx);
        fmpz_mat_init(rel, 0, n);
        return 0;
    }

    sub = _acb_vec_init(nn);
    for (i = 0; i < nn; i++)
        acb_set(sub + i, vec + idx[i]);
    num = _lindep_all_core(subrel, sub, nn, prec);

    fmpz_mat_init(rel, num, n);
    for (i = 0; i < num; i++)
        for (j = 0; j < nn; j++)
            fmpz_set(fmpz_mat_entry(rel, i, idx[j]), fmpz_mat_entry(subrel, i, j));

    fmpz_mat_clear(subrel);
    _acb_vec_clear(sub, nn);
    flint_free(idx);
    return num;
}

/* The highest level of a relation, and the number of entries at that level. */
static slong
_relation_top(slong * count, const fmpz * c, const logside_entry_struct * L, slong n)
{
    slong j, top = -1, cnt = 0;

    for (j = 0; j < n; j++)
    {
        if (!fmpz_is_zero(c + j))
        {
            if (L[j].level > top)
            {
                top = L[j].level;
                cnt = 1;
            }
            else if (L[j].level == top)
                cnt++;
        }
    }

    *count = cnt;
    return top;
}

/*
    0: log substitution, 1: unit-coefficient monomial substitution,
    2: radical with only exp generators involved, 3: radical,
    4: unusable (2 pi i on top, or several entries at the top level).
*/
static int
_relation_kind(const fmpz * c, const logside_entry_struct * L, slong n, gr_tower_t T)
{
    slong j, jt = -1, cnt, top;

    top = _relation_top(&cnt, c, L, n);
    if (cnt != 1)
        return 4;
    for (j = 0; j < n; j++)
        if (!fmpz_is_zero(c + j) && L[j].level == top)
            jt = j;
    if (L[jt].gen < 0)
        return 4;
    if (T->gens[L[jt].gen].kind == GR_TOWER_LOG || T->gens[L[jt].gen].kind == GR_TOWER_ATAN)
        return 0;
    if (T->gens[L[jt].gen].kind == GR_TOWER_LAMBERTW)
        return fmpz_is_pm1(c + jt) ? 0 : 4;
    if (fmpz_is_pm1(c + jt))
        return 1;
    if (T->gens[L[jt].gen].kind == GR_TOWER_TAN)
        return 3;
    for (j = 0; j < n; j++)
        if (j != jt && !fmpz_is_zero(c + j) && (L[j].gen < 0 || T->gens[L[j].gen].kind != GR_TOWER_EXP))
            return 3;
    return 2;
}

/*
    Row-HNF of the candidate relations with columns ordered by decreasing
    level (2 pi i last). Each row then has, as its highest-level entry, the
    gcd of that coefficient over all relations living at or below that
    level, i.e. the best possible elimination for that generator. The rows
    are returned sorted by increasing top level, then by kind.
*/
static slong
_hnf_by_level(fmpz_mat_t rows, const fmpz_mat_t cands, const logside_entry_struct * L, slong n, gr_tower_t T)
{
    slong * order, * key;
    slong i, j, m = fmpz_mat_nrows(cands), num = 0;
    fmpz_mat_t M, H;

    if (m == 0)
    {
        fmpz_mat_init(rows, 0, n);
        return 0;
    }

    /* column order: decreasing level, 2 pi i last */
    order = flint_malloc(sizeof(slong) * n);
    key = flint_malloc(sizeof(slong) * n);
    for (j = 0; j < n; j++)
    {
        order[j] = j;
        key[j] = (L[j].gen < 0) ? -1 : L[j].level;
    }
    /* insertion sort by decreasing key */
    for (i = 1; i < n; i++)
    {
        slong o = order[i], k = key[o];
        j = i;
        while (j > 0 && key[order[j - 1]] < k)
        {
            order[j] = order[j - 1];
            j--;
        }
        order[j] = o;
    }

    fmpz_mat_init(M, m, n);
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(M, i, j), fmpz_mat_entry(cands, i, order[j]));

    fmpz_mat_init(H, m, n);
    fmpz_mat_hnf(H, M);

    fmpz_mat_init(rows, m, n);
    for (i = 0; i < m; i++)
    {
        int nonzero = 0;
        for (j = 0; j < n; j++)
            if (!fmpz_is_zero(fmpz_mat_entry(H, i, j)))
                nonzero = 1;
        if (!nonzero)
            continue;
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(rows, num, order[j]), fmpz_mat_entry(H, i, j));
        num++;
    }

    /* sort by (top level, kind) */
    {
        slong * lev = flint_malloc(sizeof(slong) * FLINT_MAX(num, 1));
        for (i = 0; i < num; i++)
        {
            slong cnt;
            lev[i] = 8 * _relation_top(&cnt, fmpz_mat_row(rows, i), L, n) + _relation_kind(fmpz_mat_row(rows, i), L, n, T);
        }
        for (i = 1; i < num; i++)
        {
            j = i;
            while (j > 0 && lev[j - 1] > lev[j])
            {
                slong t = lev[j - 1]; lev[j - 1] = lev[j]; lev[j] = t;
                fmpz_mat_swap_rows(rows, NULL, j - 1, j);
                j--;
            }
        }
        flint_free(lev);
    }

    fmpz_mat_clear(M);
    fmpz_mat_clear(H);
    flint_free(order);
    flint_free(key);

    return num;
}

/* -------------------------------------------------------------------- */
/* elimination                                                           */
/* -------------------------------------------------------------------- */

/*
    For a rational u, replaces the modulus X^c - u (stored in m, c + 1
    coefficients) by its irreducible factor over Q vanishing at the
    transcendental generator t_j (identified numerically), and returns
    the length of the new modulus (c + 1 if the factor could not be
    identified).
*/
/* p = r w^c with r an integer and w a polynomial (monomials, and squares
   up to the content); returns 1 on success */
static int
_mpoly_split_power(fmpz_t r, fmpz_mpoly_t w, const fmpz_mpoly_t p, slong c, const fmpz_mpoly_ctx_t mctx)
{
    slong nvars = mctx->minfo->nvars, i;

    if (fmpz_mpoly_is_zero(p, mctx))
        return 0;

    if (p->length == 1)
    {
        ulong * e = flint_malloc(sizeof(ulong) * FLINT_MAX(nvars, 1));
        int ok = 1;
        fmpz_mpoly_get_term_exp_ui(e, p, 0, mctx);
        for (i = 0; i < nvars; i++)
        {
            if (e[i] % c != 0)
                ok = 0;
            e[i] /= c;
        }
        if (ok)
        {
            fmpz_mpoly_get_term_coeff_fmpz(r, p, 0, mctx);
            fmpz_mpoly_zero(w, mctx);
            fmpz_mpoly_push_term_ui_ui(w, 1, e, mctx);
        }
        flint_free(e);
        return ok;
    }

    if (c == 2)
    {
        fmpz_mpoly_t q;
        int ok;
        fmpz_mpoly_init(q, mctx);
        _fmpz_vec_content(r, p->coeffs, p->length);
        if (fmpz_sgn(p->coeffs + 0) < 0)
            fmpz_neg(r, r);
        fmpz_mpoly_scalar_divexact_fmpz(q, p, r, mctx);
        ok = fmpz_mpoly_sqrt(w, q, mctx);
        fmpz_mpoly_clear(q, mctx);
        return ok;
    }

    return 0;
}

static slong _rational_radical_factor_scaled(fmpz_mpoly_q_struct * m, slong c, const fmpz_mpoly_q_t u, const fmpz_mpoly_q_t w, gr_tower_t T, slong j, gr_tower_flat_t F);

static slong
_rational_radical_factor(fmpz_mpoly_q_struct * m, slong c, const fmpz_mpoly_q_t u, gr_tower_t T, slong j, gr_tower_flat_t F)
{
    return _rational_radical_factor_scaled(m, c, u, NULL, T, j, F);
}

/*
    For u = r w^c with r rational and w a flat element with an exact
    c-th power (a monomial in the generators, say), replaces the modulus
    X^c - u by w^e f(X / w), f being the irreducible factor over Q of
    Y^c - r vanishing at t_j / w: the generator is w times a root of
    unity, of degree 1 when that root is +-1. Returns the length of the
    new modulus (c + 1 if nothing applies).
*/
static slong
_radical_factor(fmpz_mpoly_q_struct * m, slong c, const fmpz_mpoly_q_t u, gr_tower_t T, slong j, gr_tower_flat_t F)
{
    fmpz_t rn, rd;
    fmpz_mpoly_q_t w, r;
    slong len = c + 1;

    fmpz_init(rn);
    fmpz_init(rd);
    fmpz_mpoly_q_init(w, F->mctx);
    fmpz_mpoly_q_init(r, F->mctx);

    if (_mpoly_split_power(rn, fmpz_mpoly_q_numref(w), fmpz_mpoly_q_numref(u), c, F->mctx) &&
        _mpoly_split_power(rd, fmpz_mpoly_q_denref(w), fmpz_mpoly_q_denref(u), c, F->mctx) &&
        !fmpz_is_zero(rd))
    {
        fmpz_mpoly_q_canonicalise(w, F->mctx);
        fmpz_mpoly_q_set_fmpz(r, rn, F->mctx);
        fmpz_mpoly_q_div_fmpz(r, r, rd, F->mctx);
        if (!fmpz_mpoly_q_is_zero(r, F->mctx))
            len = _rational_radical_factor_scaled(m, c, r, w, T, j, F);
    }

    fmpz_mpoly_q_clear(w, F->mctx);
    fmpz_mpoly_q_clear(r, F->mctx);
    fmpz_clear(rn);
    fmpz_clear(rd);
    return len;
}

/* (the value of the generator j divided by w, for the factor selection) */
typedef struct { gr_tower_struct * T; gr_tower_flat_struct * F; slong j; const fmpz_mpoly_q_struct * w; } _get_z_scaled_arg;

static int
_get_z_scaled(acb_t z, slong prec, void * arg)
{
    _get_z_scaled_arg * a = (_get_z_scaled_arg *) arg;
    int st = gr_tower_trans_get_acb(z, a->T, a->j, prec);
    if (st == GR_SUCCESS && a->w != NULL)
    {
        acb_t wv;
        acb_init(wv);
        st = gr_tower_flat_get_acb(wv, a->w, prec, a->F);
        if (st == GR_SUCCESS)
            acb_div(z, z, wv, prec);
        acb_clear(wv);
    }
    return st;
}

static slong
_rational_radical_factor_scaled(fmpz_mpoly_q_struct * m, slong c, const fmpz_mpoly_q_t u, const fmpz_mpoly_q_t w, gr_tower_t T, slong j, gr_tower_flat_t F)
{
    fmpz_poly_t f;
    fmpz_poly_factor_t fac;
    fmpz_t num, den;
    slong i, len = c + 1;

    fmpz_init(num);
    fmpz_init(den);
    fmpz_mpoly_get_fmpz(num, fmpz_mpoly_q_numref(u), F->mctx);
    fmpz_mpoly_get_fmpz(den, fmpz_mpoly_q_denref(u), F->mctx);

    /* den X^c - num */
    fmpz_poly_init(f);
    fmpz_poly_set_coeff_fmpz(f, c, den);
    fmpz_neg(num, num);
    fmpz_poly_set_coeff_fmpz(f, 0, num);

    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, f);

    if (fac->num > 1)
    {
        _get_z_scaled_arg a = { T, F, j, w };
        slong found = _gr_tower_select_fmpz_factor(fac, _get_z_scaled, &a, T);

        if (found >= 0)
        {
            const fmpz_poly_struct * g = fac->p + found;
            fmpz_t lc;
            fmpz_init(lc);
            fmpz_set(lc, g->coeffs + g->length - 1);
            for (i = 0; i < g->length; i++)
            {
                fmpz_mpoly_q_set_fmpz(m + i, g->coeffs + i, F->mctx);
                fmpz_mpoly_q_div_fmpz(m + i, m + i, lc, F->mctx);
                /* (w^(e - i) for the modulus w^e f(X / w)) */
                if (w != NULL && i < g->length - 1)
                {
                    fmpz_mpoly_q_t wp;
                    fmpz_mpoly_q_init(wp, F->mctx);
                    fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(wp), fmpz_mpoly_q_numref(w), g->length - 1 - i, F->mctx);
                    fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(wp), fmpz_mpoly_q_denref(w), g->length - 1 - i, F->mctx);
                    fmpz_mpoly_q_mul(m + i, m + i, wp, F->mctx);
                    fmpz_mpoly_q_clear(wp, F->mctx);
                }
            }
            fmpz_clear(lc);
            len = g->length;
        }
    }

    fmpz_poly_factor_clear(fac);
    fmpz_poly_clear(f);
    fmpz_clear(num);
    fmpz_clear(den);

    return len;
}

/*
    A relation c_g u_g + c_j u_j = 0 between the arguments of exactly two
    exponential generators (no 2 pi i term), with |c_g| > 1: instead of
    the radical modulus X^{c_g} - exp(u_j)^{-c_j}, both exponentials are
    expressed as powers of the exponential of the primitive argument
    u* = u_j / m with m = |c_g| / gcd(c_g, c_j): u_j = m u* and
    u_g = e u* with e = -sign(c_g) c_j / gcd. The new generator is
    inserted before the generator j; the two exponentials get moduli of
    degree 1, so the degree of the tower does not grow (as it would with
    exp(x/2) adjoined over exp(x)). Returns 1 if this was done.
*/
static int
_primitive_exp(const fmpz * c, slong jt, const logside_entry_struct * L, slong n, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong j, jo = -1, others = 0;
    slong cg, cj, s, m, e, gid_j, gid_g, gid_new, p;
    ulong dd;
    fmpz_mpoly_q_t ustar;
    fmpz_mpoly_q_struct mm[2];
    int status;

    for (j = 0; j < n; j++)
    {
        if (j != jt && !fmpz_is_zero(c + j))
        {
            others++;
            if (L[j].gen >= 0 && T->gens[L[j].gen].kind == GR_TOWER_EXP)
                jo = j;
        }
    }

    if (others != 1 || jo < 0 || !fmpz_fits_si(c + jt) || !fmpz_fits_si(c + jo))
        return 0;

    cg = fmpz_get_si(c + jt);
    cj = fmpz_get_si(c + jo);
    if (FLINT_ABS(cg) <= 1)
        return 0;

    dd = n_gcd(FLINT_ABS(cg), FLINT_ABS(cj));
    m = FLINT_ABS(cg) / dd;
    s = (cg > 0) ? 1 : -1;
    e = -s * (cj / (slong) dd);

    gid_j = T->gens[L[jo].gen].gid;
    gid_g = T->gens[L[jt].gen].gid;

    fmpz_mpoly_q_init(ustar, F->mctx);
    fmpz_mpoly_q_div_si(ustar, &L[jo].val, m, F->mctx);
    status = gr_tower_adjoin_exp_flat(T, ustar, F->mctx, NULL);
    fmpz_mpoly_q_clear(ustar, F->mctx);
    if (status != GR_SUCCESS)
        return 0;

    gid_new = T->gens[T->num_gens - 1].gid;
    p = gr_tower_gid_order(T, gid_j);
    _gr_tower_move_gen(T, T->num_gens - 1, p);

    /* exp(u_j) = exp(u*)^m */
    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(mm + 0, F->mctx);
    fmpz_mpoly_q_init(mm + 1, F->mctx);
    fmpz_mpoly_q_gen(mm + 0, GR_TOWER_FLAT_VAR_D(F, gr_tower_gid_order(T, gid_new)), F->mctx);
    fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(mm + 0), fmpz_mpoly_q_numref(mm + 0), m, F->mctx);
    fmpz_mpoly_q_neg(mm + 0, mm + 0, F->mctx);
    fmpz_mpoly_q_one(mm + 1, F->mctx);
    _gr_tower_make_algebraic(T, gr_tower_gid_order(T, gid_j), mm, 2, F->mctx, GR_TOWER_STATUS_PROVEN);
    fmpz_mpoly_q_clear(mm + 0, F->mctx);
    fmpz_mpoly_q_clear(mm + 1, F->mctx);

    /* exp(u_g) = exp(u*)^e */
    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(mm + 0, F->mctx);
    fmpz_mpoly_q_init(mm + 1, F->mctx);
    fmpz_mpoly_q_gen(mm + 0, GR_TOWER_FLAT_VAR_D(F, gr_tower_gid_order(T, gid_new)), F->mctx);
    if (e >= 0)
        fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(mm + 0), fmpz_mpoly_q_numref(mm + 0), e, F->mctx);
    else
    {
        fmpz_mpoly_q_inv(mm + 0, mm + 0, F->mctx);
        fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(mm + 0), fmpz_mpoly_q_denref(mm + 0), -e, F->mctx);
    }
    fmpz_mpoly_q_neg(mm + 0, mm + 0, F->mctx);
    fmpz_mpoly_q_one(mm + 1, F->mctx);
    _gr_tower_make_algebraic(T, gr_tower_gid_order(T, gid_g), mm, 2, F->mctx, GR_TOWER_STATUS_PROVEN);
    fmpz_mpoly_q_clear(mm + 0, F->mctx);
    fmpz_mpoly_q_clear(mm + 1, F->mctx);

    return 1;
}

/*
    Whether f(w) = w e^w - z has a unique zero in the ball D (Krawczyk:
    K = c - f(c)/f'(c) + (1 - f'(D)/f'(c)) (D - c) inside D, c the center).
*/
static int
_lambertw_unique_root(const acb_t D, const acb_t z, slong prec)
{
    acb_t c, fc, dfc, dfD, K, t;
    int res;

    acb_init(c);
    acb_init(fc);
    acb_init(dfc);
    acb_init(dfD);
    acb_init(K);
    acb_init(t);

    acb_get_mid(c, D);

    /* f(c) = c e^c - z, f'(w) = e^w (1 + w) */
    acb_exp(t, c, prec);
    acb_mul(fc, c, t, prec);
    acb_sub(fc, fc, z, prec);
    acb_add_ui(dfc, c, 1, prec);
    acb_mul(dfc, dfc, t, prec);

    acb_exp(t, D, prec);
    acb_add_ui(dfD, D, 1, prec);
    acb_mul(dfD, dfD, t, prec);

    acb_div(K, fc, dfc, prec);
    acb_sub(K, c, K, prec);
    /* (1 - f'(D)/f'(c)) (D - c) */
    {
        acb_t q, dc;
        acb_init(q);
        acb_init(dc);
        acb_div(q, dfD, dfc, prec);
        acb_sub_ui(q, q, 1, prec);
        acb_neg(q, q);
        acb_sub(dc, D, c, prec);
        acb_mul(q, q, dc, prec);
        acb_add(K, K, q, prec);
        acb_clear(q);
        acb_clear(dc);
    }

    res = acb_is_finite(K) && acb_contains(D, K) && !acb_equal(D, K);

    acb_clear(c);
    acb_clear(fc);
    acb_clear(dfc);
    acb_clear(dfD);
    acb_clear(K);
    acb_clear(t);
    return res;
}

/*
    A relation c_g w + sum_{j != g} c_j l_j = 0 with w = W_k(z) on top
    (c_g = +- 1): the candidate w* = -(sum_{j != g} c_j l_j) / c_g lies in
    the tower below w. It is a solution of w e^w = z when
    z = w* prod_{j != g} (e^{l_j})^(-c_j / c_g), which is decided exactly;
    it is the value of the generator (the branch) when a ball containing
    both enclosures contains a unique solution.
*/
static int
_lambertw_verify(const fmpz * c, slong jt, const fmpz_mpoly_q_t relation, const logside_entry_struct * L, slong n,
    gr_tower_flat_t F, slong top, slong depth)
{
    gr_tower_struct * T = F->T;
    gr_tower_gen_struct * g = T->gens + L[jt].gen;
    fmpz_mpoly_q_t wstar, prod, t, z;
    slong j, cg = fmpz_get_si(c + jt), cost = 0;
    int ok;

    for (j = 0; j < n; j++)
    {
        if (j != jt && !fmpz_is_zero(c + j))
        {
            slong terms;
            fmpz_mpoly_q_init(t, F->mctx);
            _expside(t, L + j, F);
            terms = fmpz_mpoly_q_numref(t)->length + fmpz_mpoly_q_denref(t)->length;
            cost += FLINT_ABS(fmpz_get_si(c + j)) * terms;
            fmpz_mpoly_q_clear(t, F->mctx);
        }
    }
    if (cost > GR_TOWER_OPTION(F->T, GR_TOWER_OPT_RELATION_COST_LIMIT))
        return 0;

    fmpz_mpoly_q_init(wstar, F->mctx);
    fmpz_mpoly_q_init(prod, F->mctx);
    fmpz_mpoly_q_init(t, F->mctx);
    fmpz_mpoly_q_init(z, F->mctx);

    /* w* = -(relation - c_g w) / c_g */
    fmpz_mpoly_q_mul_fmpz(t, &L[jt].val, c + jt, F->mctx);
    fmpz_mpoly_q_sub(wstar, relation, t, F->mctx);
    fmpz_mpoly_q_div_si(wstar, wstar, -cg, F->mctx);

    /* w* e^{w*} - z */
    fmpz_mpoly_q_one(prod, F->mctx);
    for (j = 0; j < n; j++)
    {
        if (j != jt && !fmpz_is_zero(c + j))
        {
            _expside(t, L + j, F);
            _mul_pow_si(prod, t, -fmpz_get_si(c + j) * cg, F);
        }
    }
    fmpz_mpoly_q_mul(prod, prod, wstar, F->mctx);
    gr_tower_flat_convert(z, &g->arg.data, g->arg.mctx, F);
    fmpz_mpoly_q_sub(prod, prod, z, F->mctx);

    ok = (_gr_tower_decide_zero_flat(prod, F, top, depth + 1) == T_TRUE);

    if (ok)
    {
        /* the branch: a ball around both enclosures with a unique root */
        acb_t a, b, zz, D;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(zz);
        acb_init(D);
        ok = 0;
        gr_tower_flat_ensure(F);
        for (prec = 64; prec <= GR_TOWER_OPTION(F->T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT) && !ok; prec *= 2)
        {
            mag_t r;
            if (gr_tower_trans_get_acb(a, T, g->index, prec) != GR_SUCCESS)
                break;
            gr_tower_flat_convert(t, wstar, F->mctx, F);
            if (gr_tower_flat_get_acb(b, t, prec, F) != GR_SUCCESS)
                break;
            gr_tower_flat_convert(t, z, F->mctx, F);
            if (gr_tower_flat_get_acb(zz, t, prec, F) != GR_SUCCESS)
                break;
            if (!acb_overlaps(a, b))
                break;
            acb_union(D, a, b, prec);
            mag_init(r);
            {
                /* r = 2 max(rad(Re D), rad(Im D)) */
                mag_max(r, arb_radref(acb_realref(D)), arb_radref(acb_imagref(D)));
            }
            mag_mul_2exp_si(r, r, 1);
            mag_add_ui_2exp_si(r, r, 1, -prec / 2);
            arb_add_error_mag(acb_realref(D), r);
            arb_add_error_mag(acb_imagref(D), r);
            mag_clear(r);
            ok = _lambertw_unique_root(D, zz, prec);
        }
        acb_clear(a);
        acb_clear(b);
        acb_clear(zz);
        acb_clear(D);
    }

    fmpz_mpoly_q_clear(wstar, F->mctx);
    fmpz_mpoly_q_clear(prod, F->mctx);
    fmpz_mpoly_q_clear(t, F->mctx);
    fmpz_mpoly_q_clear(z, F->mctx);
    return ok;
}

/*
    Verifies the candidate relation sum c_j l_j = 0 by a zero test in the
    tower below its top generator; if it holds, eliminates that generator.
    Returns 1 if the tower was changed.
*/
static int
_try_eliminate(const fmpz * c, const logside_entry_struct * L, slong n, gr_tower_flat_t F, slong depth)
{
    gr_tower_struct * T = F->T;
    slong j, jt = -1, cnt, top;
    gr_tower_gen_struct * g;
    fmpz_mpoly_q_t relation, t;
    fmpz_mpoly_ctx_struct * mctx;
    ulong version, structure_version;
    int ok = 0;

    top = _relation_top(&cnt, c, L, n);
    if (cnt != 1)
        return 0;
    for (j = 0; j < n; j++)
        if (!fmpz_is_zero(c + j) && L[j].level == top)
            jt = j;
    if (L[jt].gen < 0)
        return 0;
    g = T->gens + L[jt].gen;
    if (g->kind != GR_TOWER_EXP && g->kind != GR_TOWER_LOG && g->kind != GR_TOWER_TAN && g->kind != GR_TOWER_ATAN &&
        g->kind != GR_TOWER_LAMBERTW)
        return 0;
    /* (W with a coefficient > 1: e^w would be a radical) */
    if (g->kind == GR_TOWER_LAMBERTW && !fmpz_is_pm1(c + jt))
        return 0;
    /* (tan(u) with a coefficient > 1: a modulus of higher degree over a
       field containing i; not done here) */
    if (g->kind == GR_TOWER_TAN && !fmpz_is_pm1(c + jt))
        return 0;

    mctx = F->mctx;   /* (the context of relation and t; F->mctx may change) */
    fmpz_mpoly_q_init(relation, mctx);
    fmpz_mpoly_q_init(t, mctx);

    for (j = 0; j < n; j++)
    {
        if (!fmpz_is_zero(c + j))
        {
            fmpz_mpoly_q_mul_fmpz(t, &L[j].val, c + j, mctx);
            fmpz_mpoly_q_add(relation, relation, t, mctx);
        }
    }

    version = T->version;
    structure_version = T->structure_version;

    if (g->kind == GR_TOWER_EXP || g->kind == GR_TOWER_TAN)
    {
        /* the relation lives in the tower below g */
        ok = (_gr_tower_decide_zero_flat(relation, F, top, depth + 1) == T_TRUE);
    }
    else if (g->kind == GR_TOWER_LAMBERTW)
    {
        ok = _lambertw_verify(c, jt, relation, L, n, F, top, depth);
    }
    else
    {
        /* u^{c_g} prod e^{l_j c_j} = 1 in the tower below g, where
           g = log(u); the branch is pinned down by |relation| < 2 */
        fmpz_mpoly_q_t prod;
        slong cost = 0;

        /* (the product is not formed when it would be enormous) */
        for (j = 0; j < n; j++)
        {
            if (!fmpz_is_zero(c + j))
            {
                slong terms;
                _expside(t, L + j, F);
                terms = fmpz_mpoly_q_numref(t)->length + fmpz_mpoly_q_denref(t)->length;
                cost += FLINT_ABS(fmpz_get_si(c + j)) * terms;
            }
        }
        if (cost > GR_TOWER_OPTION(F->T, GR_TOWER_OPT_RELATION_COST_LIMIT))
        {
            fmpz_mpoly_q_clear(relation, mctx);
            fmpz_mpoly_q_clear(t, mctx);
            return 0;
        }

        fmpz_mpoly_q_init(prod, F->mctx);
        fmpz_mpoly_q_one(prod, F->mctx);
        for (j = 0; j < n; j++)
        {
            if (!fmpz_is_zero(c + j))
            {
                _expside(t, L + j, F);
                _mul_pow_si(prod, t, fmpz_get_si(c + j), F);
            }
        }
        fmpz_mpoly_q_sub_si(prod, prod, 1, F->mctx);

        ok = (_gr_tower_decide_zero_flat(prod, F, top, depth + 1) == T_TRUE);

        /* (the zero test may have changed the tower, and with it the
           context of relation: restart below rather than evaluate it) */
        if (ok && (T->version != version || T->structure_version != structure_version || F->mctx != mctx))
            ok = 0;

        if (ok)
        {
            /* |relation| < 2: the value is an integer multiple of 2 pi i,
               so this pins it down to 0 (evaluated with increasing
               precision, since the terms may be huge and cancel) */
            acb_t z;
            mag_t m;
            slong wp;

            acb_init(z);
            mag_init(m);
            ok = 0;
            for (wp = 64; wp <= GR_TOWER_OPTION(F->T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT) * 4; wp *= 2)
            {
                if (gr_tower_flat_get_acb(z, relation, wp, F) != GR_SUCCESS)
                    break;
                acb_get_mag(m, z);
                if (mag_cmp_2exp_si(m, 1) < 0)
                {
                    ok = 1;
                    break;
                }
                acb_get_mag_lower(m, z);
                if (mag_cmp_2exp_si(m, 1) >= 0)
                    break;
            }
            acb_clear(z);
            mag_clear(m);
        }

        fmpz_mpoly_q_clear(prod, F->mctx);
    }

    /* the verification may have changed the tower (eliminations below,
       refinements, reorderings): then the log side is stale, but progress
       has been made; restart */
    if (T->version != version || T->structure_version != structure_version || F->mctx != mctx)
    {
        fmpz_mpoly_q_clear(relation, mctx);
        fmpz_mpoly_q_clear(t, mctx);
        return 1;
    }

    if (ok)
    {
        slong cg = fmpz_get_si(c + jt);
        slong d = g->def_order;

        if (g->kind == GR_TOWER_LOG || g->kind == GR_TOWER_ATAN || g->kind == GR_TOWER_LAMBERTW)
        {
            /* g = -(sum_{j != jt} c_j l_j) / c_g: modulus X - expr (for
               atan, l_g = 2 i g: g = -(sum) / (2 i c_g)) */
            fmpz_mpoly_q_struct m[2];

            fmpz_mpoly_q_init(m + 0, F->mctx);
            fmpz_mpoly_q_init(m + 1, F->mctx);
            fmpz_mpoly_q_mul_fmpz(t, &L[jt].val, c + jt, F->mctx);
            fmpz_mpoly_q_sub(m + 0, relation, t, F->mctx);    /* sum over j != jt */
            fmpz_mpoly_q_div_si(m + 0, m + 0, cg, F->mctx);   /* -expr */
            if (g->kind == GR_TOWER_ATAN)
            {
                /* -expr / (2 i) = -expr (-i) / 2 */
                _flat_i(t, F, _find_i(T));
                fmpz_mpoly_q_mul(m + 0, m + 0, t, F->mctx);
                fmpz_mpoly_q_div_si(m + 0, m + 0, -2, F->mctx);
                GR_MUST_SUCCEED(gr_tower_flat_reduce(m + 0, F));
            }
            fmpz_mpoly_q_one(m + 1, F->mctx);

            ok = _gr_tower_make_algebraic(T, d, m, 2, F->mctx, GR_TOWER_STATUS_PROVEN);

            fmpz_mpoly_q_clear(m + 0, F->mctx);
            fmpz_mpoly_q_clear(m + 1, F->mctx);
        }
        else if (g->kind == GR_TOWER_TAN)
        {
            /* e^{2 i u} = w = prod_{j != jt} e^{l_j}^{-c_j / c_g}, and
               tan(u) = -i (w - 1) / (w + 1) */
            fmpz_mpoly_q_t w, a, b;
            fmpz_mpoly_q_struct m[2];

            fmpz_mpoly_q_init(w, F->mctx);
            fmpz_mpoly_q_init(a, F->mctx);
            fmpz_mpoly_q_init(b, F->mctx);
            fmpz_mpoly_q_one(w, F->mctx);
            for (j = 0; j < n; j++)
            {
                if (j != jt && !fmpz_is_zero(c + j))
                {
                    _expside(t, L + j, F);
                    _mul_pow_si(w, t, -fmpz_get_si(c + j) * cg, F);   /* (c_g = +/- 1) */
                }
            }
            fmpz_mpoly_q_sub_si(a, w, 1, F->mctx);
            fmpz_mpoly_q_add_si(b, w, 1, F->mctx);
            fmpz_mpoly_q_div(a, a, b, F->mctx);
            _flat_i(t, F, _find_i(T));
            fmpz_mpoly_q_mul(a, a, t, F->mctx);          /* i (w - 1)/(w + 1) = -tan(u) */

            fmpz_mpoly_q_init(m + 0, F->mctx);
            fmpz_mpoly_q_init(m + 1, F->mctx);
            fmpz_mpoly_q_set(m + 0, a, F->mctx);         /* X + i (w-1)/(w+1) */
            fmpz_mpoly_q_one(m + 1, F->mctx);

            if (gr_tower_flat_rationalize(m + 0, F) == GR_SUCCESS)
            {
                GR_MUST_SUCCEED(gr_tower_flat_reduce(m + 0, F));
                ok = _gr_tower_make_algebraic(T, d, m, 2, F->mctx, GR_TOWER_STATUS_PROVEN);
            }
            else
                ok = 0;

            fmpz_mpoly_q_clear(m + 0, F->mctx);
            fmpz_mpoly_q_clear(m + 1, F->mctx);
            fmpz_mpoly_q_clear(w, F->mctx);
            fmpz_mpoly_q_clear(a, F->mctx);
            fmpz_mpoly_q_clear(b, F->mctx);
        }
        else if (_primitive_exp(c, jt, L, n, F))
        {
            /* c_g u_g + c_j u_j = 0 with two exponentials: both are
               powers of the new generator exp(u_j / (c_g / gcd)) */
        }
        else
        {
            /* x_g^{c_g} = u = prod_{j != jt} e^{l_j}^{-c_j} */
            fmpz_mpoly_q_t u;
            slong cc = FLINT_ABS(cg), i;
            fmpz_mpoly_q_struct * m;

            fmpz_mpoly_q_init(u, F->mctx);
            fmpz_mpoly_q_one(u, F->mctx);
            for (j = 0; j < n; j++)
            {
                if (j != jt && !fmpz_is_zero(c + j))
                {
                    _expside(t, L + j, F);
                    _mul_pow_si(u, t, -fmpz_get_si(c + j), F);
                }
            }
            if (cg < 0)
                fmpz_mpoly_q_inv(u, u, F->mctx);

            /* the modulus must have coefficients with denominators free
               of algebraic variables (a nested chain is built from it):
               an elimination which cannot be rationalized within the size
               limit is not performed */
            if (gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(u), F))
            {
                int st = gr_tower_flat_rationalize(u, F);

                /* a zero divisor met in the flat rationalization: the
                   nested division records it in the quotient ring, the
                   tower is refined (dynamic evaluation) and the
                   rationalization tried again */
                if (st == GR_DOMAIN)
                {
                    slong k = gr_tower_flat_alg_level(u, F);
                    gr_ctx_struct * Fk = gr_tower_field_at(T, k);
                    gr_ptr tmp;
                    GR_TMP_INIT(tmp, Fk);
                    (void) gr_tower_flat_get_nested_at(tmp, u, k, F);
                    GR_TMP_CLEAR(tmp, Fk);
                    if (gr_tower_refine(T))
                    {
                        gr_tower_flat_ensure(F);
                        st = gr_tower_flat_rationalize(u, F);
                    }
                }

                if (st != GR_SUCCESS)
                {
                    fmpz_mpoly_q_clear(u, F->mctx);
                    fmpz_mpoly_q_clear(relation, mctx);
                    fmpz_mpoly_q_clear(t, mctx);
                    return 0;
                }
            }

            /* modulus X^cc - u, or the factor over Q vanishing at the
               generator when u is a rational number (e.g. roots of unity) */
            m = flint_malloc(sizeof(fmpz_mpoly_q_struct) * (cc + 1));
            for (i = 0; i <= cc; i++)
                fmpz_mpoly_q_init(m + i, F->mctx);
            fmpz_mpoly_q_neg(m + 0, u, F->mctx);
            fmpz_mpoly_q_one(m + cc, F->mctx);

            {
                slong len = cc + 1;

                if (cc > 1 && fmpz_mpoly_q_is_fmpq(u, F->mctx))
                    len = _rational_radical_factor(m, cc, u, T, g->index, F);
                else if (cc > 1)
                    len = _radical_factor(m, cc, u, T, g->index, F);

                ok = _gr_tower_make_algebraic(T, d, m, len, F->mctx,
                    (len == 2) ? GR_TOWER_STATUS_PROVEN : GR_TOWER_STATUS_DYNAMIC);
            }

            for (i = 0; i <= cc; i++)
                fmpz_mpoly_q_clear(m + i, F->mctx);
            flint_free(m);
            fmpz_mpoly_q_clear(u, F->mctx);
        }
    }

    fmpz_mpoly_q_clear(relation, mctx);
    fmpz_mpoly_q_clear(t, mctx);

    return ok;
}

/* -------------------------------------------------------------------- */
/* the real angular search                                               */
/* -------------------------------------------------------------------- */

/*
    Without i in the tower, the relations involving the trigonometric
    generators are searched among the real angles: A = 2 u for tan(u),
    A = 2 atan(u), and pi (entry with gen = -1, which cannot be the top
    of a usable relation). A relation sum c_j A_j = 0 is verified through
    the points (cos A_j, sin A_j) on the unit circle, which are rational
    functions of the generators: prod (cos A_j + i sin A_j)^{c_j} = 1,
    computed with pairs of flat elements, plus |sum c_j A_j| < 1 for the
    branch; or, when the top generator is tan(u) (whose angle 2u lies in
    the tower below it), directly as sum c_j A_j = 0 below it.
*/

static slong
_angleside(logside_entry_struct ** L_out, gr_tower_flat_t F, slong limit, const int * mask)
{
    gr_tower_struct * T = F->T;
    logside_entry_struct * L;
    slong d, n = 0, ntrig = 0;

    gr_tower_flat_ensure(F);
    L = flint_malloc(sizeof(logside_entry_struct) * (limit + 1));

    for (d = 0; d < limit; d++)
    {
        const gr_tower_gen_struct * g = T->gens + d;

        if (mask != NULL && !mask[d] && g->kind != GR_TOWER_PI)
            continue;

        if (g->kind == GR_TOWER_TAN || g->kind == GR_TOWER_ATAN)
        {
            fmpz_mpoly_q_init(&L[n].val, F->mctx);
            if (g->kind == GR_TOWER_TAN)
                gr_tower_flat_convert(&L[n].val, &g->arg.data, g->arg.mctx, F);
            else
                fmpz_mpoly_q_gen(&L[n].val, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            fmpz_mpoly_q_mul_si(&L[n].val, &L[n].val, 2, F->mctx);
            L[n].level = d;
            L[n].gen = d;
            n++;
            ntrig++;
        }
        else if (g->kind == GR_TOWER_PI)
        {
            fmpz_mpoly_q_init(&L[n].val, F->mctx);
            fmpz_mpoly_q_gen(&L[n].val, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            L[n].level = d;
            L[n].gen = -1;
            n++;
        }
    }

    if (ntrig == 0)
    {
        slong i;
        for (i = 0; i < n; i++)
            fmpz_mpoly_q_clear(&L[i].val, F->mctx);
        n = 0;
    }

    *L_out = L;
    return n;
}

/* (c, s) for the angle entry e: (-1, 0) for pi */
static void
_angle_point(fmpz_mpoly_q_t c, fmpz_mpoly_q_t s, const logside_entry_struct * e, gr_tower_flat_t F)
{
    if (e->gen < 0)
    {
        fmpz_mpoly_q_set_si(c, -1, F->mctx);
        fmpz_mpoly_q_zero(s, F->mctx);
    }
    else
    {
        fmpz_mpoly_q_t q, q2, den;
        fmpz_mpoly_q_init(q, F->mctx);
        fmpz_mpoly_q_init(q2, F->mctx);
        fmpz_mpoly_q_init(den, F->mctx);
        _trig_q(q, e->gen, F);
        fmpz_mpoly_q_mul(q2, q, q, F->mctx);
        fmpz_mpoly_q_add_si(den, q2, 1, F->mctx);
        fmpz_mpoly_q_neg(c, q2, F->mctx);
        fmpz_mpoly_q_add_si(c, c, 1, F->mctx);
        fmpz_mpoly_q_div(c, c, den, F->mctx);
        fmpz_mpoly_q_mul_si(s, q, 2, F->mctx);
        fmpz_mpoly_q_div(s, s, den, F->mctx);
        fmpz_mpoly_q_clear(q, F->mctx);
        fmpz_mpoly_q_clear(q2, F->mctx);
        fmpz_mpoly_q_clear(den, F->mctx);
    }
}

/* (c, s) = (c, s) * (a, b)^e on the unit circle (e may be negative:
   the inverse is the conjugate) */
static void
_point_mul_pow(fmpz_mpoly_q_t c, fmpz_mpoly_q_t s, const fmpz_mpoly_q_t a, const fmpz_mpoly_q_t b, slong e, gr_tower_flat_t F)
{
    fmpz_mpoly_q_t bb, t1, t2, t3;
    slong k;

    fmpz_mpoly_q_init(bb, F->mctx);
    fmpz_mpoly_q_init(t1, F->mctx);
    fmpz_mpoly_q_init(t2, F->mctx);
    fmpz_mpoly_q_init(t3, F->mctx);
    if (e < 0)
        fmpz_mpoly_q_neg(bb, b, F->mctx);
    else
        fmpz_mpoly_q_set(bb, b, F->mctx);

    for (k = 0; k < FLINT_ABS(e); k++)
    {
        /* (c + i s)(a + i bb) = (c a - s bb) + i (c bb + s a) */
        fmpz_mpoly_q_mul(t1, c, a, F->mctx);
        fmpz_mpoly_q_mul(t2, s, bb, F->mctx);
        fmpz_mpoly_q_sub(t1, t1, t2, F->mctx);
        fmpz_mpoly_q_mul(t2, c, bb, F->mctx);
        fmpz_mpoly_q_mul(t3, s, a, F->mctx);
        fmpz_mpoly_q_add(s, t2, t3, F->mctx);
        fmpz_mpoly_q_swap(c, t1, F->mctx);
        (void) gr_tower_flat_reduce(c, F);
        (void) gr_tower_flat_reduce(s, F);
    }

    fmpz_mpoly_q_clear(bb, F->mctx);
    fmpz_mpoly_q_clear(t1, F->mctx);
    fmpz_mpoly_q_clear(t2, F->mctx);
    fmpz_mpoly_q_clear(t3, F->mctx);
}

/*
    The coefficients (low to high) of the real polynomial
    C Im((1 + i X)^m) - S Re((1 + i X)^m), m = 2 |c|, whose roots include
    tan(u) when cos(2 u c) = C, sin(2 u c) = S. Returns its length.
*/
static slong
_tan_multiple_poly(fmpz_mpoly_q_struct * p, const fmpz_mpoly_q_t C, const fmpz_mpoly_q_t S, slong m, gr_tower_flat_t F)
{
    slong k;
    fmpz_t b;
    fmpz_init(b);

    for (k = 0; k <= m; k++)
    {
        /* (i X)^k binomial(m, k): real for even k (sign (-1)^(k/2)),
           imaginary for odd k (sign (-1)^((k-1)/2)) */
        fmpz_bin_uiui(b, m, k);
        if (k % 2 == 0)
        {
            if ((k / 2) % 2 == 1)
                fmpz_neg(b, b);
            fmpz_mpoly_q_mul_fmpz(p + k, S, b, F->mctx);
            fmpz_mpoly_q_neg(p + k, p + k, F->mctx);
        }
        else
        {
            if (((k - 1) / 2) % 2 == 1)
                fmpz_neg(b, b);
            fmpz_mpoly_q_mul_fmpz(p + k, C, b, F->mctx);
        }
    }

    fmpz_clear(b);
    return m + 1;
}

/*
    A relation sum_j c_j A_j = 0 among the angles of tangents tan(u_j)
    and arctangents (no pi) whose top entry is a tangent tan(u_g) with
    |c_g| > 1 and no other tangent to swap it with: the lattice spanned
    by the angles (of rank k - 1 for k entries) gets a basis B_1, ...,
    B_{k-1} (Hermite normal form), each a rational combination of the
    angles; new generators tan(B_r / 2) are inserted before the tangents
    involved, and each of these becomes a rational function of the new
    generators through its integer coordinates (addition formulas on the
    unit circle), so that the tower degree does not grow (rather than
    tan(u_g) becoming a root of a polynomial of degree 2 |c_g| with
    extraneous roots). Returns 1 if this was done.
*/
static int
_basis_tan(const fmpz * c, slong jt, const logside_entry_struct * L, slong n, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong * J, k = 0, i, r, p, * gid_old, * gid_new;
    slong cg;
    fmpz_mat_t A, H, U, Ui;
    fmpz_t den;
    int ok = 1, has_atan;

    for (i = 0; i < n; i++)
    {
        if (!fmpz_is_zero(c + i))
        {
            if (L[i].gen < 0 || !fmpz_fits_si(c + i) || FLINT_ABS(fmpz_get_si(c + i)) > 64)
                return 0;
            k++;
        }
    }
    if (k < 2)
        return 0;

    J = flint_malloc(sizeof(slong) * k);
    gid_old = flint_malloc(sizeof(slong) * k);
    gid_new = flint_malloc(sizeof(slong) * k);
    for (i = 0, r = 0; i < n; i++)
        if (!fmpz_is_zero(c + i))
            J[r++] = i;

    cg = fmpz_get_si(c + jt);

    /* progress: the k - 1 new tangents replace the tangents involved
       (the arctangents stay). Without an arctangent the count of
       transcendental generators decreases, whatever the new arguments;
       with one arctangent it is unchanged, and the new arguments must
       differ from the old ones (see below); with more it would grow. */
    {
        slong ntan = 0;
        for (r = 0; r < k; r++)
            if (T->gens[L[J[r]].gen].kind == GR_TOWER_TAN)
                ntan++;
        if (ntan == 0 || k - ntan > 1)
        {
            flint_free(J); flint_free(gid_old); flint_free(gid_new);
            return 0;
        }
        has_atan = (k - ntan == 1);
    }

    /* the position: before the first tangent involved; the arguments of
       the new generators (combinations of the arguments of all the
       entries) must lie below it */
    p = WORD_MAX;
    for (r = 0; r < k; r++)
        if (T->gens[L[J[r]].gen].kind == GR_TOWER_TAN)
            p = FLINT_MIN(p, L[J[r]].gen);
    for (r = 0; r < k && ok; r++)
    {
        slong lev;
        fmpz_mpoly_q_t w;
        fmpz_mpoly_q_init(w, F->mctx);
        if (T->gens[L[J[r]].gen].kind == GR_TOWER_TAN)
        {
            gr_tower_flat_convert(w, &T->gens[L[J[r]].gen].arg.data, T->gens[L[J[r]].gen].arg.mctx, F);
            lev = gr_tower_flat_level(w, F);
        }
        else
            lev = L[J[r]].gen + 1;   /* (an arctangent enters through its value) */
        if (lev > p)
            ok = 0;
        fmpz_mpoly_q_clear(w, F->mctx);
    }

    if (!ok)
    {
        flint_free(J); flint_free(gid_old); flint_free(gid_new);
        return 0;
    }

    /* coordinates over the other entries E = J \ {jt}, scaled by c_g:
       the entry j != jt is c_g e_j, the entry jt is -(c_j)_j */
    fmpz_mat_init(A, k, k - 1);
    fmpz_mat_init(H, k, k - 1);
    fmpz_mat_init(U, k, k);
    fmpz_mat_init(Ui, k, k);
    fmpz_init(den);
    {
        slong col = 0, rt = -1;
        for (r = 0; r < k; r++)
            if (J[r] == jt)
                rt = r;
        fmpz_mat_zero(A);
        for (r = 0; r < k; r++)
        {
            if (r == rt)
                continue;
            fmpz_set_si(fmpz_mat_entry(A, r, col), cg);
            fmpz_neg(fmpz_mat_entry(A, rt, col), c + J[r]);
            col++;
        }
    }

    fmpz_mat_hnf_transform(H, U, A);
    if (!fmpz_mat_inv(Ui, den, U) || !fmpz_is_pm1(den))
        ok = 0;
    if (ok && fmpz_is_one(den) == 0)
        fmpz_mat_neg(Ui, Ui);

    for (r = 0; r < k; r++)
        gid_old[r] = T->gens[L[J[r]].gen].gid;

    /* the new generators tan(B_r / 2), B_r = (1/c_g) sum_j H[r][j] A_j */
    for (r = 0; r < k - 1 && ok; r++)
    {
        fmpz_mpoly_q_t w, t;
        fmpz_mpoly_ctx_struct * wm;
        slong col = 0, i2;

        gr_tower_flat_ensure(F);
        wm = F->mctx;
        fmpz_mpoly_q_init(w, wm);
        fmpz_mpoly_q_init(t, wm);
        for (i2 = 0; i2 < k; i2++)
        {
            if (J[i2] == jt)
                continue;
            if (!fmpz_is_zero(fmpz_mat_entry(H, r, col)))
            {
                /* the angle A_j = L[j].val, converted to the current
                   context through its generator (the log side may be
                   stale after an adjunction) */
                slong d = gr_tower_gid_order(T, gid_old[i2]);
                const gr_tower_gen_struct * gg = T->gens + d;
                if (gg->kind == GR_TOWER_TAN || gg->def_kind == GR_TOWER_TAN)
                    gr_tower_flat_convert(t, &gg->arg.data, gg->arg.mctx, F);
                else
                    fmpz_mpoly_q_gen(t, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
                /* t = A_j / 2 */
                fmpz_mpoly_q_mul_fmpz(t, t, fmpz_mat_entry(H, r, col), F->mctx);
                fmpz_mpoly_q_add(w, w, t, F->mctx);
            }
            col++;
        }
        fmpz_mpoly_q_div_si(w, w, cg, F->mctx);

        /* (with an arctangent involved, no progress when the new argument
           is that of a tangent involved: tan(atan(v)/3) with the relation
           3 A_g = A_a has the basis A_a/3, the generator itself; the
           degree-3 modulus is the outcome then. Without one, a new
           tangent with an old argument is harmless: the old generator
           becomes algebraic over it, and the count of transcendental
           generators still decreases.) */
        if (has_atan)
        {
            slong i2;
            gr_tower_flat_reduce(w, F);
            for (i2 = 0; i2 < k && ok; i2++)
            {
                slong d = gr_tower_gid_order(T, gid_old[i2]);
                const gr_tower_gen_struct * gg = T->gens + d;
                if (gg->kind == GR_TOWER_TAN && gg->arg.mctx != NULL)
                {
                    fmpz_mpoly_q_t t2;
                    fmpz_mpoly_q_init(t2, F->mctx);
                    gr_tower_flat_convert(t2, &gg->arg.data, gg->arg.mctx, F);
                    gr_tower_flat_reduce(t2, F);
                    if (fmpz_mpoly_q_equal(t2, w, F->mctx))
                        ok = 0;
                    fmpz_mpoly_q_clear(t2, F->mctx);
                }
            }
        }
        if (!ok)
        {
            /* (generators already adjoined in this round stay: they are
               harmless, and the relation search goes on without them) */
            fmpz_mpoly_q_clear(w, wm);
            fmpz_mpoly_q_clear(t, wm);
            break;
        }

        if (gr_tower_adjoin_tan_flat(T, w, wm, NULL) != GR_SUCCESS)
            ok = 0;
        else
            gid_new[r] = T->gens[T->num_gens - 1].gid;

        fmpz_mpoly_q_clear(w, wm);
        fmpz_mpoly_q_clear(t, wm);
    }

    if (ok)
    {
        /* the new generators before the first tangent involved */
        for (r = 0; r < k - 1; r++)
        {
            slong pos = WORD_MAX, i2;
            for (i2 = 0; i2 < k; i2++)
            {
                slong d = gr_tower_gid_order(T, gid_old[i2]);
                if (T->gens[d].kind == GR_TOWER_TAN)
                    pos = FLINT_MIN(pos, d);
            }
            _gr_tower_move_gen(T, gr_tower_gid_order(T, gid_new[r]), pos);
        }

        /* each tangent involved: A_j = sum_r n_jr B_r, n = rows of U^-1 */
        for (i = 0; i < k; i++)
        {
            slong d = gr_tower_gid_order(T, gid_old[i]);
            fmpz_mpoly_q_t pc, ps, a, b, q, q2, dd;
            fmpz_mpoly_q_struct m[2];

            if (T->gens[d].kind != GR_TOWER_TAN)
                continue;

            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_init(pc, F->mctx); fmpz_mpoly_q_init(ps, F->mctx);
            fmpz_mpoly_q_init(a, F->mctx); fmpz_mpoly_q_init(b, F->mctx);
            fmpz_mpoly_q_init(q, F->mctx); fmpz_mpoly_q_init(q2, F->mctx);
            fmpz_mpoly_q_init(dd, F->mctx);
            fmpz_mpoly_q_one(pc, F->mctx);
            fmpz_mpoly_q_zero(ps, F->mctx);

            for (r = 0; r < k - 1; r++)
            {
                slong e = fmpz_get_si(fmpz_mat_entry(Ui, i, r));
                if (e == 0)
                    continue;
                /* the point of B_r: t = tan(B_r / 2) */
                fmpz_mpoly_q_gen(q, GR_TOWER_FLAT_VAR_D(F, gr_tower_gid_order(T, gid_new[r])), F->mctx);
                fmpz_mpoly_q_mul(q2, q, q, F->mctx);
                fmpz_mpoly_q_add_si(dd, q2, 1, F->mctx);
                fmpz_mpoly_q_neg(a, q2, F->mctx);
                fmpz_mpoly_q_add_si(a, a, 1, F->mctx);
                fmpz_mpoly_q_div(a, a, dd, F->mctx);
                fmpz_mpoly_q_mul_si(b, q, 2, F->mctx);
                fmpz_mpoly_q_div(b, b, dd, F->mctx);
                _point_mul_pow(pc, ps, a, b, e, F);
            }

            /* tan(A_j / 2) = sin / (1 + cos) */
            fmpz_mpoly_q_init(m + 0, F->mctx);
            fmpz_mpoly_q_init(m + 1, F->mctx);
            fmpz_mpoly_q_add_si(a, pc, 1, F->mctx);
            fmpz_mpoly_q_div(m + 0, ps, a, F->mctx);
            fmpz_mpoly_q_neg(m + 0, m + 0, F->mctx);
            fmpz_mpoly_q_one(m + 1, F->mctx);
            ok = _gr_tower_make_algebraic(T, d, m, 2, F->mctx, GR_TOWER_STATUS_PROVEN);
            fmpz_mpoly_q_clear(m + 0, F->mctx);
            fmpz_mpoly_q_clear(m + 1, F->mctx);

            fmpz_mpoly_q_clear(pc, F->mctx); fmpz_mpoly_q_clear(ps, F->mctx);
            fmpz_mpoly_q_clear(a, F->mctx); fmpz_mpoly_q_clear(b, F->mctx);
            fmpz_mpoly_q_clear(q, F->mctx); fmpz_mpoly_q_clear(q2, F->mctx);
            fmpz_mpoly_q_clear(dd, F->mctx);
        }
    }

    fmpz_mat_clear(A);
    fmpz_mat_clear(H);
    fmpz_mat_clear(U);
    fmpz_mat_clear(Ui);
    fmpz_clear(den);
    flint_free(J);
    flint_free(gid_old);
    flint_free(gid_new);
    return ok;
}

static int
_try_eliminate_angle(const fmpz * c, const logside_entry_struct * L, slong n, gr_tower_flat_t F, slong depth)
{
    gr_tower_struct * T = F->T;
    slong j, jt = -1, cnt, top, cg, d;
    gr_tower_gen_struct * g;
    fmpz_mpoly_q_t relation, t, pc, ps, a, b;
    ulong version, structure_version;
    int ok = 0;

    top = _relation_top(&cnt, c, L, n);
    if (cnt != 1)
        return 0;
    for (j = 0; j < n; j++)
        if (!fmpz_is_zero(c + j) && L[j].level == top)
            jt = j;
    if (L[jt].gen < 0 || !fmpz_fits_si(c + jt))
        return 0;
    g = T->gens + L[jt].gen;
    if (g->kind != GR_TOWER_TAN && g->kind != GR_TOWER_ATAN)
        return 0;
    for (j = 0; j < n; j++)
        if (!fmpz_fits_si(c + j) || FLINT_ABS(fmpz_get_si(c + j)) > 1000)
            return 0;

    cg = fmpz_get_si(c + jt);
    d = g->def_order;

    fmpz_mpoly_q_init(relation, F->mctx);
    fmpz_mpoly_q_init(t, F->mctx);
    fmpz_mpoly_q_init(pc, F->mctx);
    fmpz_mpoly_q_init(ps, F->mctx);
    fmpz_mpoly_q_init(a, F->mctx);
    fmpz_mpoly_q_init(b, F->mctx);

    for (j = 0; j < n; j++)
    {
        if (!fmpz_is_zero(c + j))
        {
            fmpz_mpoly_q_mul_fmpz(t, &L[j].val, c + j, F->mctx);
            fmpz_mpoly_q_add(relation, relation, t, F->mctx);
        }
    }

    version = T->version;
    structure_version = T->structure_version;

    if (g->kind == GR_TOWER_TAN)
    {
        /* the angles live in the tower below tan(u) */
        ok = (_gr_tower_decide_zero_flat(relation, F, top, depth + 1) == T_TRUE);
    }
    else
    {
        /* prod (cos A_j + i sin A_j)^{c_j} = 1, and |relation| < 1 */
        fmpz_mpoly_q_one(pc, F->mctx);
        fmpz_mpoly_q_zero(ps, F->mctx);
        for (j = 0; j < n; j++)
        {
            if (!fmpz_is_zero(c + j))
            {
                _angle_point(a, b, L + j, F);
                _point_mul_pow(pc, ps, a, b, fmpz_get_si(c + j), F);
            }
        }
        fmpz_mpoly_q_sub_si(pc, pc, 1, F->mctx);
        ok = (_gr_tower_decide_zero_flat(ps, F, top, depth + 1) == T_TRUE);
        if (ok && T->version == version && T->structure_version == structure_version)
            ok = (_gr_tower_decide_zero_flat(pc, F, top, depth + 1) == T_TRUE);

        if (ok && T->version == version && T->structure_version == structure_version)
        {
            acb_t z;
            mag_t m;
            slong wp;

            acb_init(z);
            mag_init(m);
            ok = 0;
            for (wp = 64; wp <= GR_TOWER_OPTION(F->T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT) * 4; wp *= 2)
            {
                if (gr_tower_flat_get_acb(z, relation, wp, F) != GR_SUCCESS)
                    break;
                acb_get_mag(m, z);
                if (mag_cmp_2exp_si(m, 0) < 0)
                {
                    ok = 1;
                    break;
                }
                acb_get_mag_lower(m, z);
                if (mag_cmp_2exp_si(m, 0) >= 0)
                    break;
            }
            acb_clear(z);
            mag_clear(m);
        }
    }

    if (T->version != version || T->structure_version != structure_version)
    {
        ok = 1;   /* (progress: restart) */
        goto cleanup;
    }

    if (ok)
    {
        if (g->kind == GR_TOWER_ATAN)
        {
            /* 2 g c_g = -sum_{j != jt} c_j A_j */
            fmpz_mpoly_q_struct m[2];
            fmpz_mpoly_q_init(m + 0, F->mctx);
            fmpz_mpoly_q_init(m + 1, F->mctx);
            fmpz_mpoly_q_mul_fmpz(t, &L[jt].val, c + jt, F->mctx);
            fmpz_mpoly_q_sub(m + 0, relation, t, F->mctx);
            fmpz_mpoly_q_div_si(m + 0, m + 0, 2 * cg, F->mctx);
            fmpz_mpoly_q_one(m + 1, F->mctx);
            ok = _gr_tower_make_algebraic(T, d, m, 2, F->mctx, GR_TOWER_STATUS_PROVEN);
            fmpz_mpoly_q_clear(m + 0, F->mctx);
            fmpz_mpoly_q_clear(m + 1, F->mctx);
        }
        else
        {
            /* the point Q of c_g A_g = -sum_{j != jt} c_j A_j */
            slong cc = FLINT_ABS(cg), len, i;
            fmpz_mpoly_q_struct * m;

            /* another tan(v) with a unit coefficient, which tan(u) can
               precede (its argument lying below v): the order is
               swapped, so that tan(v) is eliminated linearly (tan(2x) in
               terms of tan(x), rather than tan(x) as a root of a
               polynomial over tan(2x)) */
            if (cc > 1)
            {
                slong jbest = -1, lev;
                gr_tower_flat_ensure(F);
                {
                    fmpz_mpoly_q_t w;
                    fmpz_mpoly_q_init(w, F->mctx);
                    gr_tower_flat_convert(w, &g->arg.data, g->arg.mctx, F);
                    lev = gr_tower_flat_level(w, F);
                    fmpz_mpoly_q_clear(w, F->mctx);
                }
                for (j = 0; j < n; j++)
                    if (j != jt && L[j].gen >= 0 && fmpz_is_pm1(c + j) &&
                        T->gens[L[j].gen].kind == GR_TOWER_TAN && lev <= L[j].gen &&
                        (jbest < 0 || L[j].gen < L[jbest].gen))
                        jbest = j;
                if (jbest >= 0)
                {
                    _gr_tower_move_gen(T, d, L[jbest].gen);
                    ok = 1;
                    goto cleanup;
                }
            }

            /* tangents with non-unit coefficients: a basis of their
               angles */
            if (cc > 1 && _basis_tan(c, jt, L, n, F))
            {
                ok = 1;
                goto cleanup;
            }

            fmpz_mpoly_q_one(pc, F->mctx);
            fmpz_mpoly_q_zero(ps, F->mctx);
            for (j = 0; j < n; j++)
            {
                if (j != jt && !fmpz_is_zero(c + j))
                {
                    _angle_point(a, b, L + j, F);
                    _point_mul_pow(pc, ps, a, b, -fmpz_get_si(c + j), F);
                }
            }
            /* for c_g < 0: the angle of |c_g| A_g is the opposite */
            if (cg < 0)
                fmpz_mpoly_q_neg(ps, ps, F->mctx);

            if (cc == 1)
            {
                /* tan(u) = sin(2u) / (1 + cos(2u)) */
                m = flint_malloc(sizeof(fmpz_mpoly_q_struct) * 2);
                fmpz_mpoly_q_init(m + 0, F->mctx);
                fmpz_mpoly_q_init(m + 1, F->mctx);
                fmpz_mpoly_q_add_si(t, pc, 1, F->mctx);
                fmpz_mpoly_q_div(m + 0, ps, t, F->mctx);
                fmpz_mpoly_q_neg(m + 0, m + 0, F->mctx);
                fmpz_mpoly_q_one(m + 1, F->mctx);
                len = 2;
            }
            else
            {
                m = flint_malloc(sizeof(fmpz_mpoly_q_struct) * (2 * cc + 1));
                for (i = 0; i <= 2 * cc; i++)
                    fmpz_mpoly_q_init(m + i, F->mctx);
                len = _tan_multiple_poly(m, pc, ps, 2 * cc, F);
                while (len > 0 && fmpz_mpoly_q_is_zero(m + len - 1, F->mctx))
                    len--;
                for (i = 0; i < len - 1; i++)
                    fmpz_mpoly_q_div(m + i, m + i, m + len - 1, F->mctx);
                if (len >= 2)
                    fmpz_mpoly_q_one(m + len - 1, F->mctx);
            }

            /* denominators free of algebraic variables */
            for (i = 0; i < len && ok; i++)
            {
                if (gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(m + i), F) &&
                        gr_tower_flat_rationalize(m + i, F) != GR_SUCCESS)
                    ok = 0;
                else
                    (void) gr_tower_flat_reduce(m + i, F);
            }

            if (ok && len >= 2)
                ok = _gr_tower_make_algebraic(T, d, m, len, F->mctx,
                    (len == 2) ? GR_TOWER_STATUS_PROVEN : GR_TOWER_STATUS_DYNAMIC);
            else
                ok = 0;

            for (i = 0; i < ((cc == 1) ? 2 : 2 * cc + 1); i++)
                fmpz_mpoly_q_clear(m + i, F->mctx);
            flint_free(m);
        }
    }

cleanup:
    fmpz_mpoly_q_clear(relation, F->mctx);
    fmpz_mpoly_q_clear(t, F->mctx);
    fmpz_mpoly_q_clear(pc, F->mctx);
    fmpz_mpoly_q_clear(ps, F->mctx);
    fmpz_mpoly_q_clear(a, F->mctx);
    fmpz_mpoly_q_clear(b, F->mctx);
    return ok;
}

/* One round of the real angular search (towers without i). */
static int
_angle_round(gr_tower_flat_t F, slong limit, slong prec, slong depth, const int * mask)
{
    gr_tower_struct * T = F->T;
    logside_entry_struct * L;
    slong n, i, j, num_cands, num_rows;
    acb_ptr vals;
    fmpz_mat_t cands, rows;
    int eliminated = 0;

    n = _angleside(&L, F, limit, mask);
    if (n < 1)
    {
        _logside_clear(L, n, F);
        return 0;
    }

    /* without pi in the tower, pi is included numerically: a relation
       involving it makes pi adjoined (in front), and the search restart */
    {
        int has_pi = 0;
        for (i = 0; i < n; i++)
            if (L[i].gen < 0)
                has_pi = 1;

        vals = _acb_vec_init(n + 1);
        for (i = 0; i < n; i++)
            if (gr_tower_flat_get_acb(vals + i, &L[i].val, prec, F) != GR_SUCCESS)
                acb_indeterminate(vals + i);

        if (!has_pi)
        {
            acb_const_pi(vals + n, prec);
            num_cands = _lindep_all(cands, vals, n + 1, prec);
            for (i = 0; i < num_cands; i++)
            {
                slong bits = 0;
                for (j = 0; j <= n; j++)
                    bits = FLINT_MAX(bits, fmpz_bits(fmpz_mat_entry(cands, i, j)));
                if (bits <= FLINT_MAX(3, prec / (2 * (n + 1))) && !fmpz_is_zero(fmpz_mat_entry(cands, i, n)))
                {
                    fmpz_mat_clear(cands);
                    _acb_vec_clear(vals, n + 1);
                    _logside_clear(L, n, F);
                    GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
                    _gr_tower_move_to_front(T, T->num_gens - 1);
                    return 1;
                }
            }
            fmpz_mat_clear(cands);
        }

        if (n < 2)
        {
            _acb_vec_clear(vals, n + 1);
            _logside_clear(L, n, F);
            return 0;
        }
    }

    num_cands = _lindep_all(cands, vals, n, prec);

    {
        slong maxbits = FLINT_MAX(3, prec / (2 * n)), kept = 0;
        for (i = 0; i < num_cands; i++)
        {
            slong bits = 0;
            for (j = 0; j < n; j++)
                bits = FLINT_MAX(bits, fmpz_bits(fmpz_mat_entry(cands, i, j)));
            if (bits <= maxbits)
            {
                if (kept != i)
                    fmpz_mat_swap_rows(cands, NULL, kept, i);
                kept++;
            }
        }
        num_cands = kept;
    }

    {
        fmpz_mat_t cands2;
        fmpz_mat_init(cands2, num_cands, n);
        for (i = 0; i < num_cands; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(cands2, i, j), fmpz_mat_entry(cands, i, j));
        num_rows = _hnf_by_level(rows, cands2, L, n, T);
        fmpz_mat_clear(cands2);
    }

    for (i = 0; i < num_rows && !eliminated; i++)
    {
        slong cnt, top = _relation_top(&cnt, fmpz_mat_row(rows, i), L, n), jt = -1;
        for (j = 0; j < n; j++)
            if (!fmpz_is_zero(fmpz_mat_entry(rows, i, j)) && L[j].level == top)
                jt = j;
        if (cnt == 1 && jt >= 0 && L[jt].gen < 0)
        {
            /* pi on top: pi (which depends on nothing) moves to the front */
            slong dpi = _find_pi(T, limit);
            if (dpi > 0)
            {
                _gr_tower_move_to_front(T, dpi);
                eliminated = 1;
                break;
            }
            continue;
        }
        eliminated = _try_eliminate_angle(fmpz_mat_row(rows, i), L, n, F, depth);
    }

    fmpz_mat_clear(rows);
    fmpz_mat_clear(cands);
    _acb_vec_clear(vals, n + 1);
    _logside_clear(L, n, F);
    return eliminated;
}

/* -------------------------------------------------------------------- */
/* helpers shared with gamma_relations.c                                 */
/* -------------------------------------------------------------------- */

void
_gr_tower_certify_mul_pow_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong e, gr_tower_flat_t F)
{
    _mul_pow_si(res, x, e, F);
}

slong
_gr_tower_certify_find_i(gr_tower_t T)
{
    return _find_i(T);
}

slong
_gr_tower_certify_find_pi(gr_tower_t T, slong limit)
{
    return _find_pi(T, limit);
}

slong
_gr_tower_certify_rational_radical_factor(fmpz_mpoly_q_struct * m, slong c, const fmpz_mpoly_q_t u, gr_tower_t T, slong j, gr_tower_flat_t F)
{
    return _rational_radical_factor(m, c, u, T, j, F);
}

slong
_gr_tower_certify_radical_factor(fmpz_mpoly_q_struct * m, slong c, const fmpz_mpoly_q_t u, gr_tower_t T, slong j, gr_tower_flat_t F)
{
    return _radical_factor(m, c, u, T, j, F);
}

slong
_gr_tower_certify_lindep_all(fmpz_mat_t rel, acb_srcptr vec, slong n, slong prec)
{
    return _lindep_all(rel, vec, n, prec);
}

/*
    One round of the relation search among the generators with definition
    order < limit, at the given precision: returns 1 if the tower was
    changed (a relation eliminated, or a reordering), 0 otherwise.
*/
static int
_search_round(gr_tower_flat_t F, slong limit, slong prec, slong depth, const int * mask)
{
    gr_tower_struct * T = F->T;
    logside_entry_struct * L;
    slong n, i, j;
    acb_ptr vals;
    fmpz_mat_t cands, rows;
    slong num_cands, num_rows;
    int eliminated = 0;

    gr_tower_flat_ensure(F);

    /* Hurwitz zeta (polygamma) values at rational arguments, gamma
       values at rational arguments and on lines: exact relations */
    if (_gr_tower_hurwitz_round(F, limit, depth))
        return 1;
    gr_tower_flat_ensure(F);
    if (_gr_tower_elliptic_round(F, limit, depth))
        return 1;
    gr_tower_flat_ensure(F);
    if (_gr_tower_dilog_round(F, limit, depth))
        return 1;
    gr_tower_flat_ensure(F);
    if (_gr_tower_gamma_round(F, limit, depth))
        return 1;
    gr_tower_flat_ensure(F);

    n = _logside(&L, F, limit, mask);

    if (n == 0)
    {
        _logside_clear(L, n, F);
        return _angle_round(F, limit, prec, depth, mask);
    }

    vals = _acb_vec_init(n);
    for (i = 0; i < n; i++)
        if (gr_tower_flat_get_acb(vals + i, &L[i].val, prec, F) != GR_SUCCESS)
            acb_indeterminate(vals + i);

    num_cands = _lindep_all(cands, vals, n, prec);

    /* discard relations with unreasonably large coefficients: a
       genuine relation among n numbers known to prec bits has
       coefficients of about prec / n bits at most, and verifying a
       spurious one with large coefficients (powers of large elements)
       is expensive; more precision admits larger coefficients */
    {
        slong maxbits = FLINT_MAX(3, prec / (2 * n)), kept = 0;
        for (i = 0; i < num_cands; i++)
        {
            slong bits = 0;
            for (j = 0; j < n; j++)
                bits = FLINT_MAX(bits, fmpz_bits(fmpz_mat_entry(cands, i, j)));
            if (bits <= maxbits)
            {
                if (kept != i)
                    fmpz_mat_swap_rows(cands, NULL, kept, i);
                kept++;
            }
        }
        num_cands = kept;
    }

    {
        fmpz_mat_t cands2;
        fmpz_mat_init(cands2, num_cands, n);
        for (i = 0; i < num_cands; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(cands2, i, j), fmpz_mat_entry(cands, i, j));
        num_rows = _hnf_by_level(rows, cands2, L, n, T);
        fmpz_mat_clear(cands2);
    }

    for (i = 0; i < num_rows && !eliminated; i++)
    {
        if (_relation_kind(fmpz_mat_row(rows, i), L, n, T) == 4)
        {
            /* a relation with 2 pi i as its unique top entry becomes
               usable once pi and i (which depend on nothing) are
               moved to the front of the definition order */
            slong cnt, top = _relation_top(&cnt, fmpz_mat_row(rows, i), L, n);
            slong jt = -1;
            for (j = 0; j < n; j++)
                if (!fmpz_is_zero(fmpz_mat_entry(rows, i, j)) && L[j].level == top)
                    jt = j;
            if (cnt == 1 && L[jt].gen < 0)
            {
                slong dpi = _find_pi(T, limit), di = _find_i(T);
                if (dpi > 1 || di > 1)
                {
                    slong g_first = T->gens[FLINT_MAX(dpi, di)].gid;
                    slong g_second = T->gens[FLINT_MIN(dpi, di)].gid;
                    _gr_tower_move_to_front(T, gr_tower_gid_order(T, g_first));
                    _gr_tower_move_to_front(T, gr_tower_gid_order(T, g_second));
                    eliminated = 1;   /* restart with the new order */
                    break;
                }
            }
            continue;
        }
        eliminated = _try_eliminate(fmpz_mat_row(rows, i), L, n, F, depth);
    }

    fmpz_mat_clear(rows);
    fmpz_mat_clear(cands);
    _acb_vec_clear(vals, n);
    _logside_clear(L, n, F);

    /* the angles of the trigonometric generators among themselves (also
       with i present: the complex search does not handle tangents with
       coefficients other than +/- 1) */
    if (!eliminated)
        eliminated = _angle_round(F, limit, prec, depth, mask);

    return eliminated;
}

int
_gr_tower_search_relations(gr_tower_flat_t F, slong prec)
{
    int changed = 0;
    while (_search_round(F, F->T->num_gens, prec, 0, NULL))
        changed = 1;
    return changed;
}

/* -------------------------------------------------------------------- */
/* the main loop                                                         */
/* -------------------------------------------------------------------- */

/* Whether every generator has a proven status (irreducible modulus, or
   proven transcendental), so that the reduced flat representation is
   canonical. */
int
_gr_tower_all_proven(const gr_tower_t T)
{
    slong d;
    for (d = 0; d < T->num_gens; d++)
        if (T->gens[d].status != GR_TOWER_STATUS_PROVEN)
            return 0;
    return 1;
}

/* Whether the generators with definition order < limit include
   conjecturally transcendental ones. */
int
_gr_tower_has_conjectural_below(const gr_tower_t T, slong limit)
{
    slong d;
    for (d = 0; d < limit; d++)
        if (T->gens[d].kind != GR_TOWER_ALGEBRAIC && T->gens[d].status != GR_TOWER_STATUS_PROVEN)
            return 1;
    return 0;
}

truth_t
_gr_tower_decide_zero_flat(const fmpz_mpoly_q_t x_in, gr_tower_flat_t F, slong limit, slong depth)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_q_t x;
    fmpz_mpoly_ctx_struct * xctx;
    slong prec;
    truth_t res = T_UNKNOWN;

    if (depth > 64)
        return T_UNKNOWN;

    slong limit_gid = -1;

    gr_tower_flat_ensure(F);
    xctx = F->mctx;
    fmpz_mpoly_q_init(x, xctx);
    fmpz_mpoly_q_set(x, x_in, xctx);

    /* the limit is the generator at that position (eliminations may
       insert generators below it, shifting the definition order) */
    if (limit < T->num_gens)
        limit_gid = T->gens[limit].gid;

    for (prec = GR_TOWER_DEFAULT_PREC; prec <= GR_TOWER_OPTION(F->T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT); )
    {
        int code;
        int eliminated = 0;

        limit = (limit_gid >= 0) ? gr_tower_gid_order(T, limit_gid) : T->num_gens;

        /* the layout may have changed (reorderings) */
        gr_tower_flat_ensure(F);
        if (xctx != F->mctx)
        {
            fmpz_mpoly_q_t t;
            fmpz_mpoly_q_init(t, F->mctx);
            gr_tower_flat_convert(t, x, xctx, F);
            fmpz_mpoly_q_clear(x, xctx);
            xctx = F->mctx;
            fmpz_mpoly_q_init(x, xctx);
            fmpz_mpoly_q_swap(x, t, xctx);
            fmpz_mpoly_q_clear(t, xctx);
        }

        if (gr_tower_flat_reduce(x, F) != GR_SUCCESS)
            break;

        if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), F->mctx))
        {
            res = T_TRUE;
            break;
        }

        /* numerical separation */
        {
            acb_t z;
            acb_init(z);
            if (gr_tower_flat_get_acb(z, x, prec, F) == GR_SUCCESS && !acb_contains_zero(z))
            {
                acb_clear(z);
                res = T_FALSE;
                break;
            }
            acb_clear(z);
        }

        /* exact test over the algebraic part (may refine the tower) */
        {
            slong k = gr_tower_flat_alg_level(x, F);
            gr_ctx_struct * Fk = gr_tower_field_at(T, k);
            gr_ptr t;

            GR_TMP_INIT(t, Fk);
            if (gr_tower_flat_poly_get_nested_at(t, fmpz_mpoly_q_numref(x), k, F) == GR_SUCCESS)
                code = _gr_tower_field_is_zero_at(t, k, T);
            else
                code = GR_TOWER_UNKNOWN;
            GR_TMP_CLEAR(t, Fk);
        }

        if (code == GR_TOWER_ZERO)
        {
            res = T_TRUE;
            break;
        }
        if (code == GR_TOWER_NONZERO)
        {
            res = T_FALSE;
            break;
        }
        if (code == GR_TOWER_FIELD_NONZERO &&
            (!_gr_tower_has_conjectural_below(T, limit) || !_gr_tower_flat_involves_conjectural(x, F)))
        {
            res = T_FALSE;
            break;
        }

        /* search for a relation on the log side, among the generators
           x involves (with those their definitions involve): a relation
           between other generators is irrelevant, and may be a costly
           one (exp(260515 i) = exp(i)^260515 when exp(i) is around) */
        {
            int * mask = flint_calloc(FLINT_MAX(T->num_gens, 1), sizeof(int));
            _gr_tower_flat_involved(mask, x, F);
            eliminated = _search_round(F, limit, prec, depth, mask);
            flint_free(mask);
        }

        if (!eliminated)
            prec *= 2;
        /* else: restart at the same precision with the reduced tower */
    }

    fmpz_mpoly_q_clear(x, xctx);
    return res;
}

/*
    Certification for the nested representation: the tower is copied
    (the caller's nested elements must remain valid) and the search is run
    on the copy.
*/
truth_t
_gr_tower_certify_nonzero_at(gr_srcptr x, slong k, gr_tower_t T)
{
    gr_tower_t U;
    fmpz_mpoly_q_t xf, xu;
    truth_t res;

    gr_tower_flat_ensure(&T->flat);
    fmpz_mpoly_q_init(xf, T->flat.mctx);
    if (gr_tower_flat_set_nested_at(xf, x, k, &T->flat) != GR_SUCCESS)
    {
        fmpz_mpoly_q_clear(xf, T->flat.mctx);
        return T_UNKNOWN;
    }

    gr_tower_init(U, T->consts);
    gr_tower_set(U, T);
    U->keep_retired = 0;   /* (internal: no nested elements across rebuilds) */
    gr_tower_flat_ensure(&U->flat);
    fmpz_mpoly_q_init(xu, U->flat.mctx);
    _gr_tower_flat_transport(xu, xf, &T->flat, &U->flat);

    res = _gr_tower_decide_zero_flat(xu, &U->flat, U->num_gens, 0);

    fmpz_mpoly_q_clear(xu, U->flat.mctx);
    gr_tower_clear(U);
    fmpz_mpoly_q_clear(xf, T->flat.mctx);

    return res;
}

/* -------------------------------------------------------------------- */
/* shared steps of the relation rounds                                   */
/* -------------------------------------------------------------------- */

/* the largest definition order of a generator occurring in x (-1: none) */
slong
_gr_tower_flat_max_dep(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    slong d, res = -1;
    for (d = 0; d < F->T->num_gens; d++)
    {
        slong v = GR_TOWER_FLAT_VAR_D(F, d);
        if (fmpz_mpoly_degree_si(fmpz_mpoly_q_numref(x), v, F->mctx) > 0 ||
            fmpz_mpoly_degree_si(fmpz_mpoly_q_denref(x), v, F->mctx) > 0)
            res = d;
    }
    return res;
}

/* the safeguard of an elimination: whether the value of expr is that of
   the transcendental generator d (overlapping enclosures with at least
   bits (and 32) accurate bits, at precisions up to 1024, or as many as
   bits requires). When the two values are known to agree up to a P-th
   root of unity, bits > log2(P) + 2 decides: distinct P-th roots of
   unity differ by at least 4/P in absolute value. */
int
_gr_tower_flat_matches_gen_bits(const fmpz_mpoly_q_t expr, slong d, slong bits, gr_tower_flat_t F)
{
    acb_t u, v;
    slong prec, limit;
    int same = 0;

    bits = FLINT_MAX(bits, 32);
    limit = FLINT_MAX(1024, 4 * bits);

    acb_init(u);
    acb_init(v);
    for (prec = 64; prec <= limit && !same; prec *= 2)
    {
        if (gr_tower_flat_get_acb(u, expr, prec, F) != GR_SUCCESS ||
            gr_tower_trans_get_acb(v, F->T, F->T->gens[d].index, prec) != GR_SUCCESS)
            break;
        if (!acb_overlaps(u, v))
            break;
        if (acb_rel_accuracy_bits(u) > bits && acb_rel_accuracy_bits(v) > bits)
            same = 1;
    }
    acb_clear(u);
    acb_clear(v);
    return same;
}

int
_gr_tower_flat_matches_gen(const fmpz_mpoly_q_t expr, slong d, gr_tower_flat_t F)
{
    return _gr_tower_flat_matches_gen_bits(expr, d, 32, F);
}

/* makes the generator d algebraic of degree one, g_d = expr (expr free of
   g_d and of the generators after it; an algebraic denominator is
   rationalized first, in place); returns 1 on success */
int
_gr_tower_flat_eliminate_gen(slong d, fmpz_mpoly_q_t expr, gr_tower_flat_t F)
{
    fmpz_mpoly_q_struct mod[2];
    int ok;

    if (gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(expr), F) &&
        gr_tower_flat_rationalize(expr, F) != GR_SUCCESS)
        return 0;

    fmpz_mpoly_q_init(mod + 0, F->mctx);
    fmpz_mpoly_q_init(mod + 1, F->mctx);
    fmpz_mpoly_q_neg(mod + 0, expr, F->mctx);
    fmpz_mpoly_q_one(mod + 1, F->mctx);
    ok = _gr_tower_make_algebraic(F->T, d, mod, 2, F->mctx, GR_TOWER_STATUS_PROVEN);
    fmpz_mpoly_q_clear(mod + 0, F->mctx);
    fmpz_mpoly_q_clear(mod + 1, F->mctx);
    return ok;
}
