/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Legendre's relation between complete elliptic integrals, found by
    the zero test (called from the relation search in certify.c):

        E(m) K(1-m) + E(1-m) K(m) - K(m) K(1-m) = pi/2

    for m off the cuts (nonreal, or 0 < m < 1). The lazy field writes
    K(m), E(m) for |m - 1| > 1 through a = m/(m - 1) (the imaginary-modulus
    transformation: K(m) = s(a) K(a), E(m) = E(a)/s(a) with
    s(x) = sqrt(1 - x)), so a generator at a stands for m = a or
    m = a/(a - 1), and the partner for 1 - m sits at b = 1 - m or
    b = (m - 1)/m. With alpha = s(a) or 1 and beta = s(b) or 1
    accordingly, Legendre's relation reads

        beta^2 E(a) K(b) + alpha^2 E(b) K(a) - alpha^2 beta^2 K(a) K(b)
            = (pi/2) alpha beta,

    linear in each of the four generators. When they are all present,
    the latest one is eliminated (degree one); pi and the square roots
    (root generators of 1 - a, 1 - b) are adjoined as needed, and the
    result is checked numerically.

    Landen's transformation relates generators at a and b = k1^2 with
    k1 = (1 - s)/(1 + s), s = sqrt(1 - a):

        K(a) = (1 + k1) K(b),    E(a) = (1 + s) E(b) - s K(a).

    Pairs are detected by the necessary condition
    a^2 (b + 1)^2 = 4 (2 - a)^2 b (eliminating the square roots; tested
    numerically, then exactly), then s is adjoined (a root generator) and
    k1^2 = b verified exactly; the latest generator of the relation is
    eliminated.
*/

#include "fmpz_mpoly_q.h"
#include "acb.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

/* S = sqrt(c) (principal) as a flat element: a rational square directly,
   otherwise through a root generator with definition order < before
   (returns 1 if the tower changed, -1 if impossible, 0 when S is set) */
static int
_ell_sqrt(fmpz_mpoly_q_t S, gr_tower_flat_t F, const fmpz_mpoly_q_t c, slong before)
{
    slong ds, mult;
    int st;

    if (fmpz_mpoly_q_is_fmpq(c, F->mctx))
    {
        fmpz_t p, q;
        int sq;
        fmpz_init(p);
        fmpz_init(q);
        fmpz_mpoly_get_fmpz(p, fmpz_mpoly_q_numref(c), F->mctx);
        fmpz_mpoly_get_fmpz(q, fmpz_mpoly_q_denref(c), F->mctx);
        if (fmpz_sgn(q) < 0)
        {
            fmpz_neg(p, p);
            fmpz_neg(q, q);
        }
        sq = (fmpz_sgn(p) >= 0 && fmpz_is_square(p) && fmpz_is_square(q));
        if (sq)
        {
            fmpz_sqrt(p, p);
            fmpz_sqrt(q, q);
            fmpz_mpoly_q_set_fmpz(S, p, F->mctx);
            fmpz_mpoly_q_div_fmpz(S, S, q, F->mctx);
        }
        fmpz_clear(p);
        fmpz_clear(q);
        if (sq)
            return 0;
    }

    st = _gr_tower_special_root_gen(&ds, &mult, F, c, 2, before);
    if (st != 0)
        return st;
    fmpz_mpoly_q_gen(S, GR_TOWER_FLAT_VAR_D(F, ds), F->mctx);
    if (mult != 1)
        fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(S), fmpz_mpoly_q_numref(S), mult, F->mctx);
    return 0;
}

/* a generator with the definition K or E (transcendental, or eliminated
   already: algebraic with that definition) */
static int
_ell_is(const gr_tower_gen_struct * g, int kind)
{
    return (g->kind == kind || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == kind)) && g->arg.mctx != NULL;
}

static int
_ell_gen(const gr_tower_gen_struct * g, slong limit)
{
    return (g->kind == GR_TOWER_ELLIPTIC_K || g->kind == GR_TOWER_ELLIPTIC_E) && g->arg.mctx != NULL &&
        g->def_order < limit;
}

/* the generator of the given kind at argument u (flat), or -1 */
static slong
_ell_find(gr_tower_flat_t F, int kind, const fmpz_mpoly_q_t u, slong limit)
{
    gr_tower_struct * T = F->T;
    slong d;
    fmpz_mpoly_q_t x;

    fmpz_mpoly_q_init(x, F->mctx);
    for (d = 0; d < limit; d++)
    {
        const gr_tower_gen_struct * g = T->gens + d;
        if (!_ell_is(g, kind))
            continue;
        gr_tower_flat_convert(x, &g->arg.data, g->arg.mctx, F);
        fmpz_mpoly_q_sub(x, x, u, F->mctx);
        if (gr_tower_flat_reduce(x, F) == GR_SUCCESS && fmpz_mpoly_q_is_zero(x, F->mctx))
        {
            fmpz_mpoly_q_clear(x, F->mctx);
            return d;
        }
    }
    fmpz_mpoly_q_clear(x, F->mctx);
    return -1;
}

/* whether a is off the cuts: nonreal, or 0 < a < 1 (from enclosures;
   undecided cases are skipped) */
static int
_ell_off_cuts(const fmpz_mpoly_q_t a, gr_tower_flat_t F)
{
    acb_t z;
    int ok = 0;
    acb_init(z);
    if (gr_tower_flat_get_acb(z, a, 128, F) == GR_SUCCESS)
    {
        if (!arb_contains_zero(acb_imagref(z)))
            ok = 1;
        else if (arb_is_zero(acb_imagref(z)) && arb_is_positive(acb_realref(z)))
        {
            arb_sub_ui(acb_realref(z), acb_realref(z), 1, 128);
            ok = arb_is_negative(acb_realref(z));
        }
    }
    acb_clear(z);
    return ok;
}

/* the Legendre relation for the generators at a and b (ta, tb: whether
   they stand for the transformed values); returns 1 if the tower changed */
static int
_ell_legendre(gr_tower_flat_t F, slong dKa, slong dEa, slong dKb, slong dEb,
    const fmpz_mpoly_q_t a, const fmpz_mpoly_q_t b, int ta, int tb)
{
    gr_tower_struct * T = F->T;
    slong dtop = -1, dpi, dd[4], i;
    fmpz_mpoly_q_t Ka, Ea, Kb, Eb, R, A2, B2, num, den, expr, t, c;
    int st, ok = 1, changed = 0;

    /* the latest generator which is still transcendental */
    dd[0] = dKa; dd[1] = dEa; dd[2] = dKb; dd[3] = dEb;
    for (i = 0; i < 4; i++)
        if (T->gens[dd[i]].kind != GR_TOWER_ALGEBRAIC && dd[i] > dtop)
            dtop = dd[i];
    if (dtop < 0)
        return 0;
    /* (the others must come before it) */
    for (i = 0; i < 4; i++)
        if (dd[i] > dtop)
            return 0;

    st = _gr_tower_special_prepare_constants(T, 1, 1, dtop);
    if (st != 0)
        return (st > 0);

    /* the square roots s(a) = sqrt(1 - a), s(b) (checked below to come
       before dtop) */
    if (ta || tb)
    {
        fmpz_mpoly_q_t sa, sb;
        fmpz_mpoly_q_init(c, F->mctx);
        fmpz_mpoly_q_init(sa, F->mctx);
        fmpz_mpoly_q_init(sb, F->mctx);
        if (ta)
        {
            fmpz_mpoly_q_sub_si(c, a, 1, F->mctx);
            fmpz_mpoly_q_neg(c, c, F->mctx);
            st = _ell_sqrt(sa, F, c, dtop);
        }
        if (st == 0 && tb)
        {
            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_sub_si(c, b, 1, F->mctx);
            fmpz_mpoly_q_neg(c, c, F->mctx);
            st = _ell_sqrt(sb, F, c, dtop);
        }
        fmpz_mpoly_q_clear(c, F->mctx);
        fmpz_mpoly_q_clear(sa, F->mctx);
        fmpz_mpoly_q_clear(sb, F->mctx);
        if (st != 0)
            return (st > 0);
    }

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(Ka, F->mctx);
    fmpz_mpoly_q_init(Ea, F->mctx);
    fmpz_mpoly_q_init(Kb, F->mctx);
    fmpz_mpoly_q_init(Eb, F->mctx);
    fmpz_mpoly_q_init(R, F->mctx);
    fmpz_mpoly_q_init(A2, F->mctx);
    fmpz_mpoly_q_init(B2, F->mctx);
    fmpz_mpoly_q_init(num, F->mctx);
    fmpz_mpoly_q_init(den, F->mctx);
    fmpz_mpoly_q_init(expr, F->mctx);
    fmpz_mpoly_q_init(t, F->mctx);

    fmpz_mpoly_q_gen(Ka, GR_TOWER_FLAT_VAR_D(F, dKa), F->mctx);
    fmpz_mpoly_q_gen(Ea, GR_TOWER_FLAT_VAR_D(F, dEa), F->mctx);
    fmpz_mpoly_q_gen(Kb, GR_TOWER_FLAT_VAR_D(F, dKb), F->mctx);
    fmpz_mpoly_q_gen(Eb, GR_TOWER_FLAT_VAR_D(F, dEb), F->mctx);

    /* R = (pi/2) alpha beta; A2 = alpha^2, B2 = beta^2 */
    dpi = _gr_tower_certify_find_pi(T, T->num_gens);
    fmpz_mpoly_q_gen(R, GR_TOWER_FLAT_VAR_D(F, dpi), F->mctx);
    fmpz_mpoly_q_div_si(R, R, 2, F->mctx);
    fmpz_mpoly_q_one(A2, F->mctx);
    fmpz_mpoly_q_one(B2, F->mctx);
    if (ta)
    {
        gr_tower_flat_convert(A2, a, F->mctx, F);
        fmpz_mpoly_q_sub_si(A2, A2, 1, F->mctx);
        fmpz_mpoly_q_neg(A2, A2, F->mctx);
        /* (present now: no change to the tower) */
        if (_ell_sqrt(t, F, A2, dtop) != 0 || _gr_tower_flat_max_dep(t, F) >= dtop)
            ok = 0;
        fmpz_mpoly_q_mul(R, R, t, F->mctx);
    }
    if (tb && ok)
    {
        gr_tower_flat_convert(B2, b, F->mctx, F);
        fmpz_mpoly_q_sub_si(B2, B2, 1, F->mctx);
        fmpz_mpoly_q_neg(B2, B2, F->mctx);
        if (_ell_sqrt(t, F, B2, dtop) != 0 || _gr_tower_flat_max_dep(t, F) >= dtop)
            ok = 0;
        fmpz_mpoly_q_mul(R, R, t, F->mctx);
    }

    /* B2 E_a K_b + A2 E_b K_a - A2 B2 K_a K_b = R */
    if (dtop == dEa || dtop == dEb)
    {
        /* E_a = (R - A2 E_b K_a + A2 B2 K_a K_b) / (B2 K_b), and symmetrically */
        int sw = (dtop == dEb);
        fmpz_mpoly_q_struct * X2 = sw ? B2 : A2, * Y2 = sw ? A2 : B2;
        fmpz_mpoly_q_struct * Ey = sw ? Ea : Eb, * Kx = sw ? Kb : Ka, * Ky = sw ? Ka : Kb;
        fmpz_mpoly_q_mul(num, X2, Ey, F->mctx);
        fmpz_mpoly_q_mul(num, num, Kx, F->mctx);
        fmpz_mpoly_q_sub(num, R, num, F->mctx);
        fmpz_mpoly_q_mul(t, A2, B2, F->mctx);
        fmpz_mpoly_q_mul(t, t, Ka, F->mctx);
        fmpz_mpoly_q_mul(t, t, Kb, F->mctx);
        fmpz_mpoly_q_add(num, num, t, F->mctx);
        fmpz_mpoly_q_mul(den, Y2, Ky, F->mctx);
    }
    else
    {
        /* K_a = (R - B2 E_a K_b) / (A2 (E_b - B2 K_b)), and symmetrically */
        int sw = (dtop == dKb);
        fmpz_mpoly_q_struct * X2 = sw ? B2 : A2, * Y2 = sw ? A2 : B2;
        fmpz_mpoly_q_struct * Ex = sw ? Eb : Ea, * Ey = sw ? Ea : Eb, * Ky = sw ? Ka : Kb;
        fmpz_mpoly_q_mul(num, Y2, Ex, F->mctx);
        fmpz_mpoly_q_mul(num, num, Ky, F->mctx);
        fmpz_mpoly_q_sub(num, R, num, F->mctx);
        fmpz_mpoly_q_mul(den, Y2, Ky, F->mctx);
        fmpz_mpoly_q_sub(den, Ey, den, F->mctx);
        fmpz_mpoly_q_mul(den, den, X2, F->mctx);
    }

    if (ok && fmpz_mpoly_q_is_zero(den, F->mctx))
        ok = 0;
    if (ok)
    {
        fmpz_mpoly_q_div(expr, num, den, F->mctx);
        ok = (gr_tower_flat_reduce(expr, F) == GR_SUCCESS);
    }

    /* safeguard: the value of expr is that of the generator */
    if (ok)
        ok = _gr_tower_flat_matches_gen(expr, dtop, F);

    if (ok && _gr_tower_flat_eliminate_gen(dtop, expr, F))
        changed = 1;

    fmpz_mpoly_q_clear(Ka, F->mctx);
    fmpz_mpoly_q_clear(Ea, F->mctx);
    fmpz_mpoly_q_clear(Kb, F->mctx);
    fmpz_mpoly_q_clear(Eb, F->mctx);
    fmpz_mpoly_q_clear(R, F->mctx);
    fmpz_mpoly_q_clear(A2, F->mctx);
    fmpz_mpoly_q_clear(B2, F->mctx);
    fmpz_mpoly_q_clear(num, F->mctx);
    fmpz_mpoly_q_clear(den, F->mctx);
    fmpz_mpoly_q_clear(expr, F->mctx);
    fmpz_mpoly_q_clear(t, F->mctx);
    return changed;
}

/* the value m represented by a generator at a: a, or a/(a - 1) */
static int
_ell_value(fmpz_mpoly_q_t m, const fmpz_mpoly_q_t a, int t, gr_tower_flat_t F)
{
    if (!t)
    {
        fmpz_mpoly_q_set(m, a, F->mctx);
        return 1;
    }
    fmpz_mpoly_q_sub_si(m, a, 1, F->mctx);
    if (fmpz_mpoly_q_is_zero(m, F->mctx))
        return 0;
    fmpz_mpoly_q_div(m, a, m, F->mctx);
    if (gr_tower_flat_reduce(m, F) != GR_SUCCESS)
        return 0;
    if (gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(m), F) && gr_tower_flat_rationalize(m, F) != GR_SUCCESS)
        return 0;
    return 1;
}

/*
    Landen's transformation between the generators at a and b standing
    for m and m1 (ta, tb: transformed, m = a/(a - 1), with K(m) = A K(a),
    E(m) = E(a)/A, A = sqrt(1 - a); likewise B for b):

        A K(a) = (1 + k1) B K(b),
        E(a)/A = (1 + s) E(b)/B - s A K(a),

    s = sqrt(1 - m) (= 1/A when transformed), k1 = (1 - s)^2/m, verified
    k1^2 = m1. Eliminates the latest generator of the K relation if it is
    still transcendental, else of the E relation; returns 1 if the tower
    changed.
*/
static int
_ell_landen(gr_tower_flat_t F, slong dKa, slong dKb, slong dEa, slong dEb,
    const fmpz_mpoly_q_t a, const fmpz_mpoly_q_t b, int ta, int tb)
{
    gr_tower_struct * T = F->T;
    slong dtop;
    fmpz_mpoly_q_t c, A, B, S, m, m1, k1, x, y, expr;
    int st = 0, ok = 1, changed = 0, useE = 0;

    dtop = -1;
    if (dKa >= 0 && dKb >= 0 && T->gens[FLINT_MAX(dKa, dKb)].kind != GR_TOWER_ALGEBRAIC)
        dtop = FLINT_MAX(dKa, dKb);
    else if (dEa >= 0 && dEb >= 0 && dKa >= 0)
    {
        dtop = FLINT_MAX(FLINT_MAX(dEa, dEb), dKa);
        useE = 1;
        if (T->gens[dtop].kind == GR_TOWER_ALGEBRAIC)
            return 0;
    }
    if (dtop < 0)
        return 0;

    /* the square roots (adjoined before dtop when needed: restart) */
    fmpz_mpoly_q_init(c, F->mctx);
    fmpz_mpoly_q_init(S, F->mctx);
    if (ta)
    {
        fmpz_mpoly_q_sub_si(c, a, 1, F->mctx);
        fmpz_mpoly_q_neg(c, c, F->mctx);
        st = _ell_sqrt(S, F, c, dtop);
    }
    if (st == 0 && tb)
    {
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_clear(c, F->mctx);
        fmpz_mpoly_q_init(c, F->mctx);
        fmpz_mpoly_q_clear(S, F->mctx);
        fmpz_mpoly_q_init(S, F->mctx);
        fmpz_mpoly_q_sub_si(c, b, 1, F->mctx);
        fmpz_mpoly_q_neg(c, c, F->mctx);
        st = _ell_sqrt(S, F, c, dtop);
    }
    if (st == 0 && !ta)
    {
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_clear(c, F->mctx);
        fmpz_mpoly_q_init(c, F->mctx);
        fmpz_mpoly_q_clear(S, F->mctx);
        fmpz_mpoly_q_init(S, F->mctx);
        fmpz_mpoly_q_sub_si(c, a, 1, F->mctx);
        fmpz_mpoly_q_neg(c, c, F->mctx);
        st = _ell_sqrt(S, F, c, dtop);
    }
    fmpz_mpoly_q_clear(c, F->mctx);
    fmpz_mpoly_q_clear(S, F->mctx);
    if (st != 0)
        return (st > 0);

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(c, F->mctx);
    fmpz_mpoly_q_init(A, F->mctx);
    fmpz_mpoly_q_init(B, F->mctx);
    fmpz_mpoly_q_init(S, F->mctx);
    fmpz_mpoly_q_init(m, F->mctx);
    fmpz_mpoly_q_init(m1, F->mctx);
    fmpz_mpoly_q_init(k1, F->mctx);
    fmpz_mpoly_q_init(x, F->mctx);
    fmpz_mpoly_q_init(y, F->mctx);
    fmpz_mpoly_q_init(expr, F->mctx);

    /* (present now: no change to the tower) */
    fmpz_mpoly_q_one(A, F->mctx);
    fmpz_mpoly_q_one(B, F->mctx);
    if (ta)
    {
        fmpz_mpoly_q_sub_si(c, a, 1, F->mctx);
        fmpz_mpoly_q_neg(c, c, F->mctx);
        ok = ok && (_ell_sqrt(A, F, c, dtop) == 0);
    }
    if (tb)
    {
        fmpz_mpoly_q_sub_si(c, b, 1, F->mctx);
        fmpz_mpoly_q_neg(c, c, F->mctx);
        ok = ok && (_ell_sqrt(B, F, c, dtop) == 0);
    }
    ok = ok && _ell_value(m, a, ta, F) && _ell_value(m1, b, tb, F);
    if (ok)
    {
        if (ta)
            fmpz_mpoly_q_inv(S, A, F->mctx);
        else
        {
            fmpz_mpoly_q_sub_si(c, a, 1, F->mctx);
            fmpz_mpoly_q_neg(c, c, F->mctx);
            ok = (_ell_sqrt(S, F, c, dtop) == 0);
        }
    }
    if (ok)
        ok = (_gr_tower_flat_max_dep(A, F) < dtop && _gr_tower_flat_max_dep(B, F) < dtop && _gr_tower_flat_max_dep(S, F) < dtop);

    /* k1 = (1 - s)^2 / m; k1^2 = m1 exactly */
    if (ok)
    {
        fmpz_mpoly_q_sub_si(k1, S, 1, F->mctx);
        fmpz_mpoly_q_mul(k1, k1, k1, F->mctx);
        fmpz_mpoly_q_div(k1, k1, m, F->mctx);
        ok = (gr_tower_flat_reduce(k1, F) == GR_SUCCESS);
        if (ok && gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(k1), F))
            ok = (gr_tower_flat_rationalize(k1, F) == GR_SUCCESS);
    }
    if (ok)
    {
        fmpz_mpoly_q_mul(x, k1, k1, F->mctx);
        fmpz_mpoly_q_sub(x, x, m1, F->mctx);
        ok = (gr_tower_flat_num_is_zero(x, F) == T_TRUE);
    }

    if (ok)
    {
        gr_tower_flat_ensure(F);
        if (!useE)
        {
            /* A K(a) = (1 + k1) B K(b) */
            fmpz_mpoly_q_t Kx;
            fmpz_mpoly_q_init(Kx, F->mctx);
            fmpz_mpoly_q_add_si(x, k1, 1, F->mctx);
            fmpz_mpoly_q_mul(x, x, B, F->mctx);         /* (1 + k1) B */
            if (dtop == dKa)
            {
                fmpz_mpoly_q_gen(Kx, GR_TOWER_FLAT_VAR_D(F, dKb), F->mctx);
                fmpz_mpoly_q_mul(expr, x, Kx, F->mctx);
                fmpz_mpoly_q_div(expr, expr, A, F->mctx);
            }
            else
            {
                fmpz_mpoly_q_gen(Kx, GR_TOWER_FLAT_VAR_D(F, dKa), F->mctx);
                fmpz_mpoly_q_mul(expr, A, Kx, F->mctx);
                fmpz_mpoly_q_div(expr, expr, x, F->mctx);
            }
            fmpz_mpoly_q_clear(Kx, F->mctx);
        }
        else
        {
            /* E(a)/A = (1 + s) E(b)/B - s A K(a) */
            fmpz_mpoly_q_t Ea, Eb, Ka;
            fmpz_mpoly_q_init(Ea, F->mctx);
            fmpz_mpoly_q_init(Eb, F->mctx);
            fmpz_mpoly_q_init(Ka, F->mctx);
            fmpz_mpoly_q_gen(Ea, GR_TOWER_FLAT_VAR_D(F, dEa), F->mctx);
            fmpz_mpoly_q_gen(Eb, GR_TOWER_FLAT_VAR_D(F, dEb), F->mctx);
            fmpz_mpoly_q_gen(Ka, GR_TOWER_FLAT_VAR_D(F, dKa), F->mctx);
            fmpz_mpoly_q_add_si(x, S, 1, F->mctx);       /* 1 + s */
            fmpz_mpoly_q_mul(y, S, A, F->mctx);          /* s A */
            if (dtop == dEa)
            {
                /* E(a) = A ((1 + s) E(b)/B - s A K(a)) */
                fmpz_mpoly_q_mul(expr, x, Eb, F->mctx);
                fmpz_mpoly_q_div(expr, expr, B, F->mctx);
                fmpz_mpoly_q_mul(Ka, Ka, y, F->mctx);
                fmpz_mpoly_q_sub(expr, expr, Ka, F->mctx);
                fmpz_mpoly_q_mul(expr, expr, A, F->mctx);
            }
            else if (dtop == dEb)
            {
                /* E(b) = B (E(a)/A + s A K(a)) / (1 + s) */
                fmpz_mpoly_q_div(expr, Ea, A, F->mctx);
                fmpz_mpoly_q_mul(Ka, Ka, y, F->mctx);
                fmpz_mpoly_q_add(expr, expr, Ka, F->mctx);
                fmpz_mpoly_q_mul(expr, expr, B, F->mctx);
                fmpz_mpoly_q_div(expr, expr, x, F->mctx);
            }
            else
            {
                /* K(a) = ((1 + s) E(b)/B - E(a)/A) / (s A) */
                fmpz_mpoly_q_mul(expr, x, Eb, F->mctx);
                fmpz_mpoly_q_div(expr, expr, B, F->mctx);
                fmpz_mpoly_q_div(Ea, Ea, A, F->mctx);
                fmpz_mpoly_q_sub(expr, expr, Ea, F->mctx);
                if (fmpz_mpoly_q_is_zero(y, F->mctx))
                    ok = 0;
                else
                    fmpz_mpoly_q_div(expr, expr, y, F->mctx);
            }
            fmpz_mpoly_q_clear(Ea, F->mctx);
            fmpz_mpoly_q_clear(Eb, F->mctx);
            fmpz_mpoly_q_clear(Ka, F->mctx);
        }
        ok = ok && (gr_tower_flat_reduce(expr, F) == GR_SUCCESS);
    }

    /* safeguard: the value of expr is that of the generator */
    if (ok)
        ok = _gr_tower_flat_matches_gen(expr, dtop, F);

    if (ok && _gr_tower_flat_eliminate_gen(dtop, expr, F))
        changed = 1;

    fmpz_mpoly_q_clear(c, F->mctx);
    fmpz_mpoly_q_clear(A, F->mctx);
    fmpz_mpoly_q_clear(B, F->mctx);
    fmpz_mpoly_q_clear(S, F->mctx);
    fmpz_mpoly_q_clear(m, F->mctx);
    fmpz_mpoly_q_clear(m1, F->mctx);
    fmpz_mpoly_q_clear(k1, F->mctx);
    fmpz_mpoly_q_clear(x, F->mctx);
    fmpz_mpoly_q_clear(y, F->mctx);
    fmpz_mpoly_q_clear(expr, F->mctx);
    return changed;
}

/* the Landen pairs among the arguments of the K and E generators */
static int
_ell_landen_round(gr_tower_flat_t F, slong limit)
{
    gr_tower_struct * T = F->T;
    slong i, j;
    int changed = 0;

    for (i = limit - 1; i >= 0 && !changed; i--)
    {
        for (j = limit - 1; j >= 0 && !changed; j--)
        {
            const gr_tower_gen_struct * gi = T->gens + i, * gj = T->gens + j;
            fmpz_mpoly_q_t a, b, x, y;
            acb_t za, zb, w, v;
            int cand = 0;

            if (i == j)
                continue;
            if (!_ell_is(gi, GR_TOWER_ELLIPTIC_K) && !_ell_is(gi, GR_TOWER_ELLIPTIC_E))
                continue;
            if (!_ell_is(gj, GR_TOWER_ELLIPTIC_K) && !_ell_is(gj, GR_TOWER_ELLIPTIC_E))
                continue;

            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_init(a, F->mctx);
            fmpz_mpoly_q_init(b, F->mctx);
            fmpz_mpoly_q_init(x, F->mctx);
            fmpz_mpoly_q_init(y, F->mctx);
            acb_init(za);
            acb_init(zb);
            acb_init(w);
            acb_init(v);
            gr_tower_flat_convert(a, &gi->arg.data, gi->arg.mctx, F);
            gr_tower_flat_convert(b, &gj->arg.data, gj->arg.mctx, F);

            /* for the values m, m1 represented (a or a/(a - 1), ...):
               m^2 (m1 + 1)^2 - 4 (2 - m)^2 m1, numerically, then exactly */
            {
                int ta, tb;
                fmpz_mpoly_q_t mm, mm1;
                fmpz_mpoly_q_init(mm, F->mctx);
                fmpz_mpoly_q_init(mm1, F->mctx);
                for (ta = 0; ta < 2 && !changed; ta++)
                for (tb = 0; tb < 2 && !changed; tb++)
                {
                    if (!_ell_value(mm, a, ta, F) || !_ell_value(mm1, b, tb, F))
                        continue;
                    cand = 0;
                    if (gr_tower_flat_get_acb(za, mm, 128, F) == GR_SUCCESS &&
                        gr_tower_flat_get_acb(zb, mm1, 128, F) == GR_SUCCESS)
                    {
                        acb_add_ui(w, zb, 1, 128);
                        acb_mul(w, w, za, 128);
                        acb_sqr(w, w, 128);
                        acb_sub_ui(v, za, 2, 128);
                        acb_sqr(v, v, 128);
                        acb_mul(v, v, zb, 128);
                        acb_mul_2exp_si(v, v, 2);
                        acb_sub(w, w, v, 128);
                        cand = acb_contains_zero(w);
                    }
                    if (!cand || !_ell_off_cuts(mm, F))
                        continue;
                    fmpz_mpoly_q_add_si(x, mm1, 1, F->mctx);
                    fmpz_mpoly_q_mul(x, x, mm, F->mctx);
                    fmpz_mpoly_q_mul(x, x, x, F->mctx);
                    fmpz_mpoly_q_sub_si(y, mm, 2, F->mctx);
                    fmpz_mpoly_q_mul(y, y, y, F->mctx);
                    fmpz_mpoly_q_mul(y, y, mm1, F->mctx);
                    fmpz_mpoly_q_mul_si(y, y, 4, F->mctx);
                    fmpz_mpoly_q_sub(x, x, y, F->mctx);
                    if (gr_tower_flat_num_is_zero(x, F) == T_TRUE)
                    {
                        slong dKa = _ell_find(F, GR_TOWER_ELLIPTIC_K, a, limit);
                        slong dKb = _ell_find(F, GR_TOWER_ELLIPTIC_K, b, limit);
                        slong dEa = _ell_find(F, GR_TOWER_ELLIPTIC_E, a, limit);
                        slong dEb = _ell_find(F, GR_TOWER_ELLIPTIC_E, b, limit);
                        changed = _ell_landen(F, dKa, dKb, dEa, dEb, a, b, ta, tb);
                        /* (the flat context may have changed) */
                        if (!changed)
                        {
                            gr_tower_flat_ensure(F);
                        }
                    }
                }
                fmpz_mpoly_q_clear(mm, F->mctx);
                fmpz_mpoly_q_clear(mm1, F->mctx);
            }

            fmpz_mpoly_q_clear(a, F->mctx);
            fmpz_mpoly_q_clear(b, F->mctx);
            fmpz_mpoly_q_clear(x, F->mctx);
            fmpz_mpoly_q_clear(y, F->mctx);
            acb_clear(za);
            acb_clear(zb);
            acb_clear(w);
            acb_clear(v);
        }
    }
    return changed;
}

int
_gr_tower_elliptic_round(gr_tower_flat_t F, slong limit, slong depth)
{
    gr_tower_struct * T = F->T;
    slong j, n = 0;
    int changed = 0;

    (void) depth;

    {
        slong ntrans = 0;
        for (j = 0; j < limit; j++)
        {
            if (_ell_is(T->gens + j, GR_TOWER_ELLIPTIC_K) || _ell_is(T->gens + j, GR_TOWER_ELLIPTIC_E))
                n++;
            if (_ell_gen(T->gens + j, limit))
                ntrans++;
        }
        if (n < 2 || ntrans == 0)
            return 0;
    }

    if (_ell_landen_round(F, limit))
        return 1;
    if (n < 4)
        return 0;

    /* for each K(a): the value m it stands for (a or a/(a - 1)), and the
       partners for 1 - m (1 - m or (m - 1)/m) */
    for (j = limit - 1; j >= 0 && !changed; j--)
    {
        const gr_tower_gen_struct * g = T->gens + j;
        fmpz_mpoly_q_t a, m, b, x;
        slong dEa, dKb, dEb;
        int ta, tb;

        if (!_ell_is(g, GR_TOWER_ELLIPTIC_K))
            continue;

        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(a, F->mctx);
        fmpz_mpoly_q_init(m, F->mctx);
        fmpz_mpoly_q_init(b, F->mctx);
        fmpz_mpoly_q_init(x, F->mctx);
        gr_tower_flat_convert(a, &g->arg.data, g->arg.mctx, F);
        dEa = _ell_find(F, GR_TOWER_ELLIPTIC_E, a, limit);

        for (ta = 0; ta < 2 && dEa >= 0 && !changed; ta++)
        {
            if (ta == 0)
                fmpz_mpoly_q_set(m, a, F->mctx);
            else
            {
                fmpz_mpoly_q_sub_si(x, a, 1, F->mctx);
                fmpz_mpoly_q_div(m, a, x, F->mctx);
                if (gr_tower_flat_reduce(m, F) != GR_SUCCESS)
                    continue;
            }
            if (!_ell_off_cuts(m, F))
                continue;

            for (tb = 0; tb < 2 && !changed; tb++)
            {
                /* b = 1 - m, or (m - 1)/m */
                fmpz_mpoly_q_sub_si(b, m, 1, F->mctx);
                if (tb == 0)
                    fmpz_mpoly_q_neg(b, b, F->mctx);
                else
                    fmpz_mpoly_q_div(b, b, m, F->mctx);
                if (gr_tower_flat_reduce(b, F) != GR_SUCCESS)
                    continue;
                if (gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(b), F) &&
                    gr_tower_flat_rationalize(b, F) != GR_SUCCESS)
                    continue;
                dKb = _ell_find(F, GR_TOWER_ELLIPTIC_K, b, limit);
                dEb = _ell_find(F, GR_TOWER_ELLIPTIC_E, b, limit);
                if (dKb < 0 || dEb < 0 || dKb == j)
                    continue;
                changed = _ell_legendre(F, j, dEa, dKb, dEb, a, b, ta, tb);
            }
        }

        fmpz_mpoly_q_clear(a, F->mctx);
        fmpz_mpoly_q_clear(m, F->mctx);
        fmpz_mpoly_q_clear(b, F->mctx);
        fmpz_mpoly_q_clear(x, F->mctx);
    }

    return changed;
}
