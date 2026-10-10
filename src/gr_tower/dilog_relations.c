/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    The distribution relations of the dilogarithm, found by the zero test
    (called from the relation search in certify.c):

        Li_2(y^n) = n sum_{k<n} Li_2(zeta_n^k y)       (|y| <= 1, y^n != 1).

    For |y| <= 1 all the arguments lie in the closed unit disk, where the
    power series converge, so the relation holds exactly. The lazy field
    keeps Li_2 generators only at canonical arguments under the
    anharmonic group (lazy_special.c), so each value Li_2(v) above is
    sigma Li_2(x) + E with x the argument of a generator, sigma = +-1 and E
    elementary (logarithms and pi^2) along a word in the involutions

        s: Li_2(v) = -Li_2(1 - v) + pi^2/6 - log(v) log(1 - v),
        t: Li_2(v) = -Li_2(1/v) - pi^2/6 - log(-v)^2/2,
        r: Li_2(v) = -Li_2(v/(v - 1)) - log(1 - v)^2/2,
        c: Li_2(v) = -Li_2(1/v) + pi^2/3 - log(v)^2/2 - i pi log(v)   (v > 1),

    with the words of lazy_special.c (for nonreal v all of s, t, r are
    valid; real v uses r for v < 0, c for v > 1, then s). The candidate
    y runs over the orbit points with |y| <= 1 of each generator's
    argument and n over 2, 3, 4; when every value is found among the
    generators (special values such as roots of unity are not
    generators, and are skipped), the latest generator with a nonzero
    coefficient is eliminated (degree one), with the logarithms, pi, i
    and zeta_n adjoined as needed. The result is checked numerically.
*/

#include <string.h>
#include "ulong_extras.h"
#include "fmpz_mpoly_q.h"
#include "acb.h"
#include "gr_tower.h"
#include "fmpq.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

#define DILOG_MAX_N 4

/* the generators Li_2(x) */
static int
_dl_gen(const gr_tower_gen_struct * g)
{
    return g->kind == GR_TOWER_POLYLOG && g->def_param == 2 && g->arg.mctx != NULL;
}

/* one step of a word: w = the orbit point (flat); returns 0 on failure */
static int
_dl_point(fmpz_mpoly_q_t w, int step, const fmpz_mpoly_q_t v, gr_tower_flat_t F)
{
    fmpz_mpoly_q_t t;
    int ok = 1;
    fmpz_mpoly_q_init(t, F->mctx);
    if (step == 's')
    {
        fmpz_mpoly_q_sub_si(t, v, 1, F->mctx);
        fmpz_mpoly_q_neg(w, t, F->mctx);
    }
    else if (step == 't' || step == 'c')
    {
        if (fmpz_mpoly_q_is_zero(v, F->mctx))
            ok = 0;
        else
            fmpz_mpoly_q_inv(w, v, F->mctx);
    }
    else
    {
        fmpz_mpoly_q_sub_si(t, v, 1, F->mctx);
        if (fmpz_mpoly_q_is_zero(t, F->mctx))
            ok = 0;
        else
            fmpz_mpoly_q_div(w, v, t, F->mctx);
    }
    fmpz_mpoly_q_clear(t, F->mctx);
    if (ok)
        ok = (gr_tower_flat_reduce(w, F) == GR_SUCCESS);
    if (ok && gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(w), F))
        ok = (gr_tower_flat_rationalize(w, F) == GR_SUCCESS);
    return ok;
}

/* the words applicable to v (from its enclosure: nonreal, or real with
   v < 0, 0 < v < 1, v > 1); returns the number of words, 0 if undecided */
static int
_dl_words(const char ** words, const fmpz_mpoly_q_t v, gr_tower_flat_t F)
{
    acb_t z;
    int n = 0;
    acb_init(z);
    if (gr_tower_flat_get_acb(z, v, 128, F) == GR_SUCCESS)
    {
        if (!arb_contains_zero(acb_imagref(z)))
        {
            words[0] = ""; words[1] = "s"; words[2] = "t"; words[3] = "st"; words[4] = "ts"; words[5] = "r";
            n = 6;
        }
        else if (arb_is_zero(acb_imagref(z)))
        {
            arb_t u;
            arb_init(u);
            arb_sub_ui(u, acb_realref(z), 1, 128);
            if (arb_is_negative(acb_realref(z)))
            {
                words[0] = "r"; words[1] = "rs";
                n = 2;
            }
            else if (arb_is_positive(acb_realref(z)) && arb_is_negative(u))
            {
                words[0] = ""; words[1] = "s";
                n = 2;
            }
            else if (arb_is_positive(u))
            {
                words[0] = "c"; words[1] = "cs";
                n = 2;
            }
            arb_clear(u);
        }
    }
    acb_clear(z);
    return n;
}

/* the generator (index into gens) whose argument is the orbit point of v
   under one of its words: sets *word; -1 if none */
static slong
_dl_find(const char ** word, const fmpz_mpoly_q_t v, gr_tower_flat_t F, const slong * gens, slong ngens)
{
    const char * words[6];
    int nw = _dl_words(words, v, F), i;
    slong j, res = -1;
    fmpz_mpoly_q_t w, x;
    acb_t a, b;

    fmpz_mpoly_q_init(w, F->mctx);
    fmpz_mpoly_q_init(x, F->mctx);
    acb_init(a);
    acb_init(b);

    for (i = 0; i < nw && res < 0; i++)
    {
        const char * p;
        int ok = 1;
        fmpz_mpoly_q_set(w, v, F->mctx);
        for (p = words[i]; *p && ok; p++)
            ok = _dl_point(w, *p, w, F);
        if (!ok || gr_tower_flat_get_acb(a, w, 128, F) != GR_SUCCESS)
            continue;
        for (j = 0; j < ngens && res < 0; j++)
        {
            const gr_tower_gen_struct * g = F->T->gens + gens[j];
            gr_tower_flat_convert(x, &g->arg.data, g->arg.mctx, F);
            if (gr_tower_flat_get_acb(b, x, 128, F) != GR_SUCCESS || !acb_overlaps(a, b))
                continue;
            fmpz_mpoly_q_sub(x, x, w, F->mctx);
            if (gr_tower_flat_reduce(x, F) == GR_SUCCESS && fmpz_mpoly_q_is_zero(x, F->mctx))
            {
                res = j;
                *word = words[i];
            }
        }
    }

    fmpz_mpoly_q_clear(w, F->mctx);
    fmpz_mpoly_q_clear(x, F->mctx);
    acb_clear(a);
    acb_clear(b);
    return res;
}

/* E for Li_2(v) = sigma Li_2(word(v)) + E, as a flat element; returns 1 if
   the tower changed (a logarithm, pi or i adjoined), -1 if impossible, 0
   when E is set */
static int
_dl_elementary(fmpz_mpoly_q_t E, int * sigma, const fmpz_mpoly_q_t v_in, const char * word,
    gr_tower_flat_t F, slong before)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_q_t v, w, e, la, lb, P, x;
    const char * p;
    slong dpi, di, d;
    int st = 0, sg = 1;

    /* pi (and i for the cut formula) first */
    if (*word)
    {
        st = _gr_tower_special_prepare_constants(T, strchr(word, 'c') ? 4 : 1, 1, before);
        if (st != 0)
            return st;
    }

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(v, F->mctx);
    fmpz_mpoly_q_init(w, F->mctx);
    fmpz_mpoly_q_init(e, F->mctx);
    fmpz_mpoly_q_init(la, F->mctx);
    fmpz_mpoly_q_init(lb, F->mctx);
    fmpz_mpoly_q_init(P, F->mctx);
    fmpz_mpoly_q_init(x, F->mctx);

    fmpz_mpoly_q_set(v, v_in, F->mctx);
    fmpz_mpoly_q_zero(E, F->mctx);
    if (*word)
    {
        dpi = _gr_tower_certify_find_pi(T, T->num_gens);
        fmpz_mpoly_q_gen(P, GR_TOWER_FLAT_VAR_D(F, dpi), F->mctx);
        fmpz_mpoly_q_mul(P, P, P, F->mctx);         /* pi^2 */
    }

    for (p = word; *p && st == 0; p++)
    {
        if (!_dl_point(w, *p, v, F))
        {
            st = -1;
            break;
        }
        if (*p == 's')
        {
            /* pi^2/6 - log(v) log(1 - v) */
            st = _gr_tower_special_trans_gen(&d, F, GR_TOWER_LOG, v, before);
            if (st != 0) break;
            fmpz_mpoly_q_gen(la, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            st = _gr_tower_special_trans_gen(&d, F, GR_TOWER_LOG, w, before);
            if (st != 0) break;
            fmpz_mpoly_q_gen(lb, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            fmpz_mpoly_q_mul(e, la, lb, F->mctx);
            fmpz_mpoly_q_div_si(x, P, 6, F->mctx);
            fmpz_mpoly_q_sub(e, x, e, F->mctx);
        }
        else if (*p == 't')
        {
            /* -pi^2/6 - log(-v)^2/2 */
            fmpz_mpoly_q_neg(x, v, F->mctx);
            st = _gr_tower_special_trans_gen(&d, F, GR_TOWER_LOG, x, before);
            if (st != 0) break;
            fmpz_mpoly_q_gen(la, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            fmpz_mpoly_q_mul(e, la, la, F->mctx);
            fmpz_mpoly_q_div_si(e, e, 2, F->mctx);
            fmpz_mpoly_q_div_si(x, P, 6, F->mctx);
            fmpz_mpoly_q_add(e, e, x, F->mctx);
            fmpz_mpoly_q_neg(e, e, F->mctx);
        }
        else if (*p == 'r')
        {
            /* -log(1 - v)^2/2 */
            fmpz_mpoly_q_sub_si(x, v, 1, F->mctx);
            fmpz_mpoly_q_neg(x, x, F->mctx);
            st = _gr_tower_special_trans_gen(&d, F, GR_TOWER_LOG, x, before);
            if (st != 0) break;
            fmpz_mpoly_q_gen(la, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            fmpz_mpoly_q_mul(e, la, la, F->mctx);
            fmpz_mpoly_q_div_si(e, e, 2, F->mctx);
            fmpz_mpoly_q_neg(e, e, F->mctx);
        }
        else
        {
            /* pi^2/3 - log(v)^2/2 - i pi log(v) */
            st = _gr_tower_special_trans_gen(&d, F, GR_TOWER_LOG, v, before);
            if (st != 0) break;
            fmpz_mpoly_q_gen(la, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
            fmpz_mpoly_q_mul(e, la, la, F->mctx);
            fmpz_mpoly_q_div_si(e, e, 2, F->mctx);
            fmpz_mpoly_q_div_si(x, P, 3, F->mctx);
            fmpz_mpoly_q_sub(e, x, e, F->mctx);
            di = _gr_tower_certify_find_i(T);
            if (di < 0) { st = -1; break; }
            fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, di), F->mctx);
            {
                const gr_tower_gen_struct * gi = T->gens + di;
                if (gi->def_kind == GR_TOWER_ROOT_OF_UNITY && gi->def_param != 4)
                    fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(x), fmpz_mpoly_q_numref(x), gi->def_param / 4, F->mctx);
            }
            fmpz_mpoly_q_mul(x, x, la, F->mctx);
            fmpz_mpoly_q_gen(lb, GR_TOWER_FLAT_VAR_D(F, _gr_tower_certify_find_pi(T, T->num_gens)), F->mctx);
            fmpz_mpoly_q_mul(x, x, lb, F->mctx);
            fmpz_mpoly_q_sub(e, e, x, F->mctx);
        }
        /* Li_2(v0) = sg Li_2(v) + E = sg (-Li_2(w) + e) + E */
        if (sg < 0)
            fmpz_mpoly_q_neg(e, e, F->mctx);
        fmpz_mpoly_q_add(E, E, e, F->mctx);
        if (gr_tower_flat_reduce(E, F) != GR_SUCCESS)
            st = -1;
        sg = -sg;
        fmpz_mpoly_q_swap(v, w, F->mctx);
    }

    *sigma = sg;

    fmpz_mpoly_q_clear(v, F->mctx);
    fmpz_mpoly_q_clear(w, F->mctx);
    fmpz_mpoly_q_clear(e, F->mctx);
    fmpz_mpoly_q_clear(la, F->mctx);
    fmpz_mpoly_q_clear(lb, F->mctx);
    fmpz_mpoly_q_clear(P, F->mctx);
    fmpz_mpoly_q_clear(x, F->mctx);
    return st;
}

/* special values: v in {0, 1, -1, 1/2, 2} (rational, exactly); sets the
   code (0..4) and returns 1 */
static int
_dl_special(int * code, const fmpz_mpoly_q_t v, gr_tower_flat_t F)
{
    fmpz_t p, q;
    int ok = 0;
    if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(v), F->mctx) || !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(v), F->mctx))
        return 0;
    fmpz_init(p);
    fmpz_init(q);
    fmpz_mpoly_get_fmpz(p, fmpz_mpoly_q_numref(v), F->mctx);
    fmpz_mpoly_get_fmpz(q, fmpz_mpoly_q_denref(v), F->mctx);
    if (fmpz_sgn(q) < 0)
    {
        fmpz_neg(p, p);
        fmpz_neg(q, q);
    }
    if (fmpz_is_zero(p)) { *code = 0; ok = 1; }
    else if (fmpz_equal(p, q)) { *code = 1; ok = 1; }
    else if (fmpz_is_one(q) && fmpz_equal_si(p, -1)) { *code = 2; ok = 1; }
    else if (fmpz_is_one(p) && fmpz_equal_ui(q, 2)) { *code = 3; ok = 1; }
    else if (fmpz_is_one(q) && fmpz_equal_ui(p, 2)) { *code = 4; ok = 1; }
    fmpz_clear(p);
    fmpz_clear(q);
    return ok;
}

/* res = Li_2 at the special value: 0, pi^2/6, -pi^2/12,
   pi^2/12 - log(2)^2/2, pi^2/4 - i pi log(2) (returns 1 if the tower
   changed, -1 if impossible, 0 when set) */
static int
_dl_special_value(fmpz_mpoly_q_t res, int code, gr_tower_flat_t F, slong before)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_q_t P, L, x;
    slong d;
    int st;

    if (code == 0)
    {
        fmpz_mpoly_q_zero(res, F->mctx);
        return 0;
    }
    st = _gr_tower_special_prepare_constants(T, (code == 4) ? 4 : 1, 1, before);
    if (st != 0)
        return st;
    if (code >= 3)
    {
        fmpz_mpoly_q_t two;
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(two, F->mctx);
        fmpz_mpoly_q_set_si(two, 2, F->mctx);
        st = _gr_tower_special_trans_gen(&d, F, GR_TOWER_LOG, two, before);
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_clear(two, F->mctx);
        if (st != 0)
            return st;
    }
    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(P, F->mctx);
    fmpz_mpoly_q_init(L, F->mctx);
    fmpz_mpoly_q_init(x, F->mctx);
    fmpz_mpoly_q_gen(P, GR_TOWER_FLAT_VAR_D(F, _gr_tower_certify_find_pi(T, T->num_gens)), F->mctx);
    if (code >= 3)
        fmpz_mpoly_q_gen(L, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
    fmpz_mpoly_q_mul(res, P, P, F->mctx);
    if (code == 1)
        fmpz_mpoly_q_div_si(res, res, 6, F->mctx);
    else if (code == 2)
        fmpz_mpoly_q_div_si(res, res, -12, F->mctx);
    else if (code == 3)
    {
        fmpz_mpoly_q_div_si(res, res, 12, F->mctx);
        fmpz_mpoly_q_mul(x, L, L, F->mctx);
        fmpz_mpoly_q_div_si(x, x, 2, F->mctx);
        fmpz_mpoly_q_sub(res, res, x, F->mctx);
    }
    else
    {
        slong di = _gr_tower_certify_find_i(T);
        fmpz_mpoly_q_div_si(res, res, 4, F->mctx);
        if (di < 0)
            st = -1;
        else
        {
            const gr_tower_gen_struct * gi = T->gens + di;
            fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, di), F->mctx);
            if (gi->def_kind == GR_TOWER_ROOT_OF_UNITY && gi->def_param != 4)
                fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(x), fmpz_mpoly_q_numref(x), gi->def_param / 4, F->mctx);
            fmpz_mpoly_q_mul(x, x, P, F->mctx);
            fmpz_mpoly_q_mul(x, x, L, F->mctx);
            fmpz_mpoly_q_sub(res, res, x, F->mctx);
        }
    }
    if (st == 0 && gr_tower_flat_reduce(res, F) != GR_SUCCESS)
        st = -1;
    fmpz_mpoly_q_clear(P, F->mctx);
    fmpz_mpoly_q_clear(L, F->mctx);
    fmpz_mpoly_q_clear(x, F->mctx);
    return st;
}

/* whether v (numerically) is an orbit point of some generator's argument */
static int
_dl_match_numeric(const acb_t v, gr_tower_flat_t F, const slong * gens, slong ngens)
{
    acb_t x, w, t;
    slong j;
    int i, found = 0;
    acb_init(x);
    acb_init(w);
    acb_init(t);
    for (j = 0; j < ngens && !found; j++)
    {
        const gr_tower_gen_struct * g = F->T->gens + gens[j];
        fmpz_mpoly_q_t a;
        fmpz_mpoly_q_init(a, F->mctx);
        gr_tower_flat_convert(a, &g->arg.data, g->arg.mctx, F);
        if (gr_tower_flat_get_acb(x, a, 128, F) == GR_SUCCESS)
        {
            /* the six orbit points of x */
            for (i = 0; i < 6 && !found; i++)
            {
                switch (i)
                {
                    case 0: acb_set(w, x); break;
                    case 1: acb_sub_ui(w, x, 1, 128); acb_neg(w, w); break;
                    case 2: acb_inv(w, x, 128); break;
                    case 3: acb_sub_ui(w, x, 1, 128); acb_neg(w, w); acb_inv(w, w, 128); break;
                    case 4: acb_inv(w, x, 128); acb_sub_ui(w, w, 1, 128); acb_neg(w, w); break;
                    default: acb_sub_ui(t, x, 1, 128); acb_div(w, x, t, 128); break;
                }
                found = acb_overlaps(w, v);
            }
        }
        fmpz_mpoly_q_clear(a, F->mctx);
    }
    acb_clear(x);
    acb_clear(w);
    acb_clear(t);
    return found;
}

/* whether v is numerically one of 0, 1, -1, 1/2, 2 */
static int
_dl_special_numeric(const acb_t v)
{
    static const slong num[5] = { 0, 1, -1, 1, 2 }, den[5] = { 1, 1, 1, 2, 1 };
    acb_t w;
    int i, found = 0;
    acb_init(w);
    for (i = 0; i < 5 && !found; i++)
    {
        acb_set_si(w, num[i]);
        acb_div_si(w, w, den[i], 128);
        found = acb_overlaps(w, v);
    }
    acb_clear(w);
    return found;
}

/* x = p/q exactly (q > 0)? */
static int
_dl_rat(fmpz_t p, fmpz_t q, const fmpz_mpoly_q_t x, gr_tower_flat_t F)
{
    if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(x), F->mctx) || !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(x), F->mctx))
        return 0;
    fmpz_mpoly_get_fmpz(p, fmpz_mpoly_q_numref(x), F->mctx);
    fmpz_mpoly_get_fmpz(q, fmpz_mpoly_q_denref(x), F->mctx);
    if (fmpz_sgn(q) < 0)
    {
        fmpz_neg(p, p);
        fmpz_neg(q, q);
    }
    return 1;
}

/* the canonical orbit point of a rational p/q (not 0, 1): in (0, 1/2]
   (as in the lazy field), in place; h = max(|p|, q) of the result */
static void
_dl_rat_canonical(fmpz_t h, fmpz_t p, fmpz_t q)
{
    fmpz_t t;
    fmpz_init(t);
    if (fmpz_sgn(p) < 0)
    {
        /* z / (z - 1) = p / (p - q) */
        fmpz_sub(t, q, p);
        fmpz_neg(p, p);
        fmpz_swap(q, t);
    }
    else if (fmpz_cmp(p, q) > 0)
        fmpz_swap(p, q);
    fmpz_mul_2exp(t, p, 1);
    if (fmpz_cmp(t, q) > 0)
        fmpz_sub(p, q, p);
    fmpz_gcd(t, p, q);
    fmpz_divexact(p, p, t);
    fmpz_divexact(q, q, t);
    fmpz_set(h, q);
    fmpz_clear(t);
}

/*
    expr = R (the relation, flat): solved for the latest transcendental
    generator among gens occurring linearly in it, which is eliminated
    after a numerical check; returns 1 if eliminated
*/
static int
_dl_solve_eliminate(fmpz_mpoly_q_t expr, gr_tower_flat_t F, const slong * gens, slong ngens)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_t C, R0, t;
    slong j, dtop = -1, v;
    int ok;

    gr_tower_flat_ensure(F);
    fmpz_mpoly_init(C, F->mctx);
    fmpz_mpoly_init(R0, F->mctx);
    fmpz_mpoly_init(t, F->mctx);
    ok = (gr_tower_flat_reduce(expr, F) == GR_SUCCESS);
    for (j = 0; j < ngens && ok; j++)
    {
        slong d = gens[j];
        if (T->gens[d].kind == GR_TOWER_ALGEBRAIC || d <= dtop)
            continue;
        v = GR_TOWER_FLAT_VAR_D(F, d);
        if (fmpz_mpoly_degree_si(fmpz_mpoly_q_numref(expr), v, F->mctx) == 1 &&
            fmpz_mpoly_degree_si(fmpz_mpoly_q_denref(expr), v, F->mctx) <= 0)
            dtop = d;
    }
    if (dtop < 0)
        ok = 0;
    if (ok)
    {
        v = GR_TOWER_FLAT_VAR_D(F, dtop);
        fmpz_mpoly_derivative(C, fmpz_mpoly_q_numref(expr), v, F->mctx);
        fmpz_mpoly_gen(t, v, F->mctx);
        fmpz_mpoly_mul(t, t, C, F->mctx);
        fmpz_mpoly_sub(R0, fmpz_mpoly_q_numref(expr), t, F->mctx);
        fmpz_mpoly_neg(R0, R0, F->mctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_numref(expr), R0, F->mctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_denref(expr), C, F->mctx);
        ok = !fmpz_mpoly_is_zero(fmpz_mpoly_q_denref(expr), F->mctx);
        if (ok)
        {
            fmpz_mpoly_q_canonicalise(expr, F->mctx);
            ok = (gr_tower_flat_reduce(expr, F) == GR_SUCCESS) && _gr_tower_flat_max_dep(expr, F) < dtop;
        }
    }
    fmpz_mpoly_clear(C, F->mctx);
    fmpz_mpoly_clear(R0, F->mctx);
    fmpz_mpoly_clear(t, F->mctx);

    if (ok)
        ok = _gr_tower_flat_matches_gen(expr, dtop, F);
    return ok && _gr_tower_flat_eliminate_gen(dtop, expr, F);
}

/* the relation for y (flat, |y| <= 1) and n: eliminates a generator
   (returns 1, also when the tower changed otherwise), or returns 0 */
static int
_dl_relation(gr_tower_flat_t F, const fmpz_mpoly_q_t y_in, slong n, const slong * gens, slong ngens)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_q_struct vals[DILOG_MAX_N + 1];
    const char * words[DILOG_MAX_N + 1];
    slong idx[DILOG_MAX_N + 1], * coef, k, j, top = -1, dtop = -1;
    fmpz_mpoly_q_t y, z, zk, expr, E;
    int st = 0, ok = 1, changed = 0, rat = 0;
    slong miss = -1;

    /* numerically first: y^n and zeta^k y at generators (then zeta_n is
       adjoined if missing, and the round restarts) */
    {
        acb_t a, u, zz, v;
        int all = 1;
        acb_init(a);
        acb_init(u);
        acb_init(zz);
        acb_init(v);
        {
            fmpz_t p, q;
            fmpz_init(p);
            fmpz_init(q);
            rat = _dl_rat(p, q, y_in, F);
            fmpz_clear(p);
            fmpz_clear(q);
        }
        if (gr_tower_flat_get_acb(a, y_in, 128, F) != GR_SUCCESS)
            all = 0;
        if (all)
        {
            acb_pow_ui(u, a, n, 128);
            all = _dl_match_numeric(u, F, gens, ngens) || _dl_special_numeric(u);
        }
        for (k = 0; k < n && all; k++)
        {
            acb_set_si(zz, 2 * k);
            acb_div_si(zz, zz, n, 128);
            acb_exp_pi_i(zz, zz, 128);
            acb_mul(v, zz, a, 128);
            all = _dl_match_numeric(v, F, gens, ngens) || _dl_special_numeric(v) ||
                (rat && arb_is_zero(acb_imagref(v)));
        }
        acb_clear(a);
        acb_clear(u);
        acb_clear(zz);
        acb_clear(v);
        if (!all)
            return 0;
        st = _gr_tower_special_prepare_constants(T, (n % 4 == 0) ? n : n * (4 / n_gcd(n, 4)), 1, T->num_gens);
        if (st != 0)
            return (st > 0);
    }

    coef = flint_malloc(sizeof(slong) * FLINT_MAX(ngens, 1));
    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(y, F->mctx);
    fmpz_mpoly_q_set(y, y_in, F->mctx);

    /* zeta_n (and i) present? (not adjoined unless a relation is found) */
    for (k = 0; k <= n; k++)
        fmpz_mpoly_q_init(vals + k, F->mctx);
    fmpz_mpoly_q_init(z, F->mctx);
    fmpz_mpoly_q_init(zk, F->mctx);
    fmpz_mpoly_q_init(expr, F->mctx);
    fmpz_mpoly_q_init(E, F->mctx);

    /* vals[0] = y^n, vals[1 + k] = zeta^k y */
    if (!_gr_tower_special_flat_zeta(z, n, F))
        ok = 0;
    if (ok)
    {
        fmpz_mpoly_q_one(vals + 0, F->mctx);
        for (k = 0; k < n; k++)
            fmpz_mpoly_q_mul(vals + 0, vals + 0, y, F->mctx);
        ok = (gr_tower_flat_reduce(vals + 0, F) == GR_SUCCESS);
        fmpz_mpoly_q_set(zk, y, F->mctx);
        for (k = 0; k < n && ok; k++)
        {
            fmpz_mpoly_q_set(vals + 1 + k, zk, F->mctx);
            fmpz_mpoly_q_mul(zk, zk, z, F->mctx);
            ok = (gr_tower_flat_reduce(zk, F) == GR_SUCCESS);
        }
    }

    /* every value at a generator, or a special value (idx = -1 - code) */
    for (k = 0; k <= n && ok; k++)
    {
        int code;
        idx[k] = _dl_find(words + k, vals + k, F, gens, ngens);
        if (idx[k] < 0)
        {
            if (_dl_special(&code, vals + k, F))
                idx[k] = -1 - code;
            else if (k >= 1 && rat && miss < 0 && idx[0] >= 0 && T->gens[gens[idx[0]]].kind != GR_TOWER_ALGEBRAIC)
                miss = k;
            else
                ok = 0;
        }
    }

    /* a rational value without a generator: Li_2 at its canonical orbit
       point is adjoined when of smaller height than the generator to
       eliminate (so that this terminates), then the round restarts */
    if (ok && miss >= 0)
    {
        const gr_tower_gen_struct * g0 = T->gens + gens[idx[0]];
        fmpz_t p, q, h, p0, q0, h0;
        fmpz_mpoly_q_t a0, w;
        fmpz_init(p);
        fmpz_init(q);
        fmpz_init(h);
        fmpz_init(p0);
        fmpz_init(q0);
        fmpz_init(h0);
        fmpz_mpoly_q_init(a0, F->mctx);
        fmpz_mpoly_q_init(w, F->mctx);
        gr_tower_flat_convert(a0, &g0->arg.data, g0->arg.mctx, F);
        ok = 0;
        if (_dl_rat(p, q, vals + miss, F) && _dl_rat(p0, q0, a0, F) &&
            !fmpz_is_zero(p0) && !fmpz_equal(p0, q0))
        {
            _dl_rat_canonical(h, p, q);
            _dl_rat_canonical(h0, p0, q0);
            if (fmpz_cmp(h, h0) < 0)
            {
                fmpq_t r;
                fmpq_init(r);
                fmpq_set_fmpz_frac(r, p, q);
                fmpz_mpoly_q_set_fmpq(w, r, F->mctx);
                fmpq_clear(r);
                /* (right before the generator to eliminate: after the
                   constants already there) */
                slong d0 = gens[idx[0]];
                if (gr_tower_adjoin_special_flat(T, GR_TOWER_POLYLOG, 2, w, F->mctx, NULL) == GR_SUCCESS)
                {
                    _gr_tower_move_gen(T, T->num_gens - 1, d0);
                    st = 1;
                }
            }
        }
        /* (cleared in the context they were made in) */
        fmpz_mpoly_q_clear(a0, F->mctx);
        fmpz_mpoly_q_clear(w, F->mctx);
        fmpz_clear(p);
        fmpz_clear(q);
        fmpz_clear(h);
        fmpz_clear(p0);
        fmpz_clear(q0);
        fmpz_clear(h0);
    }

    if (ok)
    {
        /* the coefficients of the generators (sigma determined by the
           word length) */
        for (j = 0; j < ngens; j++)
            coef[j] = 0;
        for (k = 0; k <= n; k++)
        {
            slong c = (k == 0) ? 1 : -n;
            if (idx[k] < 0)
                continue;
            if (strlen(words[k]) % 2 == 1)
                c = -c;
            coef[idx[k]] += c;
        }
        for (j = 0; j < ngens; j++)
            if (coef[j] != 0 && T->gens[gens[j]].kind != GR_TOWER_ALGEBRAIC && gens[j] > dtop)
            {
                dtop = gens[j];
                top = j;
            }
        if (top < 0)
            ok = 0;
    }

    /* the elementary parts, with zeta_n before dtop */
    if (ok)
    {
        st = _gr_tower_special_prepare_constants(T, (n % 4 == 0) ? n : n * (4 / n_gcd(n, 4)), 1, dtop);
        if (st != 0)
            ok = 0;
    }

    if (ok)
    {
        /* sum_k c_k (sigma_k Li_2(x_k) + E_k) = 0 */
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_zero(expr, F->mctx);
        for (k = 0; k <= n && st == 0; k++)
        {
            int sigma;
            slong c = (k == 0) ? 1 : -n;
            gr_tower_flat_ensure(F);
            if (idx[k] < 0)
                st = _dl_special_value(E, -1 - idx[k], F, dtop);
            else
                st = _dl_elementary(E, &sigma, vals + k, words[k], F, dtop);
            if (st != 0)
                break;
            fmpz_mpoly_q_mul_si(E, E, c, F->mctx);
            fmpz_mpoly_q_add(expr, expr, E, F->mctx);
        }
        if (st != 0)
            ok = 0;
    }

    if (ok)
    {
        /* R = sum_k c_k (sigma_k g_k + E_k), reduced (eliminated
           generators substituted); then solved for the latest
           transcendental generator occurring linearly in it */
        fmpz_mpoly_q_t x;
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(x, F->mctx);
        for (j = 0; j < ngens; j++)
        {
            if (coef[j] == 0)
                continue;
            fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, gens[j]), F->mctx);
            fmpz_mpoly_q_mul_si(x, x, coef[j], F->mctx);
            fmpz_mpoly_q_add(expr, expr, x, F->mctx);
        }
        fmpz_mpoly_q_clear(x, F->mctx);
        if (_dl_solve_eliminate(expr, F, gens, ngens))
            changed = 1;
    }

    if (st > 0)
        changed = 1;

    /* (the flat context may have changed: clear in the current one) */
    gr_tower_flat_ensure(F);
    for (k = 0; k <= n; k++)
        fmpz_mpoly_q_clear(vals + k, F->mctx);
    fmpz_mpoly_q_clear(y, F->mctx);
    fmpz_mpoly_q_clear(z, F->mctx);
    fmpz_mpoly_q_clear(zk, F->mctx);
    fmpz_mpoly_q_clear(expr, F->mctx);
    fmpz_mpoly_q_clear(E, F->mctx);
    flint_free(coef);
    return changed;
}


/*
    Abel's five-term relation, for real x, y in (0, 1):

        Li_2(x) + Li_2(y) - Li_2(x y) - Li_2(x (1 - y)/(1 - x y)) - Li_2(y (1 - x)/(1 - x y))
            = log((1 - x)/(1 - x y)) log((1 - y)/(1 - x y)),

    with x and y among the points in (0, 1) of the orbits of the
    generators' arguments, and every value found in an orbit (or
    special): the latest generator occurring linearly is eliminated, as
    for the distribution relations.
*/
static int
_dl_five(gr_tower_flat_t F, const slong * gens, slong ngens)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_q_struct * pts;
    acb_ptr apts;
    slong npts = 0, i, j, k;
    int changed = 0;

    gr_tower_flat_ensure(F);
    pts = flint_malloc(sizeof(fmpz_mpoly_q_struct) * 2 * ngens);
    apts = _acb_vec_init(2 * ngens);
    for (j = 0; j < ngens; j++)
    {
        const gr_tower_gen_struct * g = T->gens + gens[j];
        fmpz_mpoly_q_struct * x = pts + npts;
        arb_t u;
        fmpz_mpoly_q_init(x, F->mctx);
        gr_tower_flat_convert(x, &g->arg.data, g->arg.mctx, F);
        arb_init(u);
        if (gr_tower_flat_get_acb(apts + npts, x, 128, F) == GR_SUCCESS && arb_is_zero(acb_imagref(apts + npts)) &&
            arb_is_positive(acb_realref(apts + npts)))
        {
            arb_sub_ui(u, acb_realref(apts + npts), 1, 128);
            if (arb_is_negative(u))
            {
                /* x and 1 - x */
                fmpz_mpoly_q_init(pts + npts + 1, F->mctx);
                fmpz_mpoly_q_sub_si(pts + npts + 1, x, 1, F->mctx);
                fmpz_mpoly_q_neg(pts + npts + 1, pts + npts + 1, F->mctx);
                acb_sub_ui(apts + npts + 1, apts + npts, 1, 128);
                acb_neg(apts + npts + 1, apts + npts + 1);
                npts += 2;
                arb_clear(u);
                continue;
            }
        }
        arb_clear(u);
        fmpz_mpoly_q_clear(x, F->mctx);
    }

    for (i = 0; i < npts && !changed; i++)
    {
        for (j = i; j < npts && !changed; j++)
        {
            acb_t a, b, c;
            fmpz_mpoly_q_struct vals[5];
            const char * words[5];
            slong idx[5];
            int ok = 1, st = 0;
            static const slong cf[5] = { 1, 1, -1, -1, -1 };

            /* numerically: x y, x (1 - y)/(1 - x y), y (1 - x)/(1 - x y) */
            acb_init(a);
            acb_init(b);
            acb_init(c);
            acb_mul(a, apts + i, apts + j, 128);
            acb_sub_ui(c, a, 1, 128);
            acb_neg(c, c);
            ok = _dl_match_numeric(a, F, gens, ngens) || _dl_special_numeric(a);
            if (ok)
            {
                acb_sub_ui(b, apts + j, 1, 128);
                acb_neg(b, b);
                acb_mul(b, b, apts + i, 128);
                acb_div(b, b, c, 128);
                ok = _dl_match_numeric(b, F, gens, ngens) || _dl_special_numeric(b);
            }
            if (ok)
            {
                acb_sub_ui(b, apts + i, 1, 128);
                acb_neg(b, b);
                acb_mul(b, b, apts + j, 128);
                acb_div(b, b, c, 128);
                ok = _dl_match_numeric(b, F, gens, ngens) || _dl_special_numeric(b);
            }
            acb_clear(a);
            acb_clear(b);
            acb_clear(c);
            if (!ok)
                continue;

            /* exactly */
            gr_tower_flat_ensure(F);
            for (k = 0; k < 5; k++)
                fmpz_mpoly_q_init(vals + k, F->mctx);
            {
                fmpz_mpoly_q_t d, t;
                fmpz_mpoly_q_init(d, F->mctx);
                fmpz_mpoly_q_init(t, F->mctx);
                fmpz_mpoly_q_set(vals + 0, pts + i, F->mctx);
                fmpz_mpoly_q_set(vals + 1, pts + j, F->mctx);
                fmpz_mpoly_q_mul(vals + 2, pts + i, pts + j, F->mctx);
                fmpz_mpoly_q_sub_si(d, vals + 2, 1, F->mctx);
                fmpz_mpoly_q_neg(d, d, F->mctx);        /* 1 - x y */
                fmpz_mpoly_q_sub_si(t, pts + j, 1, F->mctx);
                fmpz_mpoly_q_neg(t, t, F->mctx);
                fmpz_mpoly_q_mul(t, t, pts + i, F->mctx);
                fmpz_mpoly_q_div(vals + 3, t, d, F->mctx);
                fmpz_mpoly_q_sub_si(t, pts + i, 1, F->mctx);
                fmpz_mpoly_q_neg(t, t, F->mctx);
                fmpz_mpoly_q_mul(t, t, pts + j, F->mctx);
                fmpz_mpoly_q_div(vals + 4, t, d, F->mctx);
                fmpz_mpoly_q_clear(d, F->mctx);
                fmpz_mpoly_q_clear(t, F->mctx);
            }
            for (k = 0; k < 5 && ok; k++)
            {
                int code;
                ok = (gr_tower_flat_reduce(vals + k, F) == GR_SUCCESS);
                if (!ok)
                    break;
                idx[k] = _dl_find(words + k, vals + k, F, gens, ngens);
                if (idx[k] < 0)
                {
                    if (_dl_special(&code, vals + k, F))
                        idx[k] = -1 - code;
                    else
                        ok = 0;
                }
            }

            if (ok)
            {
                /* sum_k c_k (sigma_k g_k + E_k) - L1 L2 */
                fmpz_mpoly_q_t expr, E, x, L1, L2;
                slong dtop = -1, dl;
                gr_tower_flat_ensure(F);
                fmpz_mpoly_q_init(expr, F->mctx);
                fmpz_mpoly_q_init(E, F->mctx);
                fmpz_mpoly_q_init(x, F->mctx);
                fmpz_mpoly_q_init(L1, F->mctx);
                fmpz_mpoly_q_init(L2, F->mctx);

                for (k = 0; k < 5; k++)
                    if (idx[k] >= 0 && T->gens[gens[idx[k]]].kind != GR_TOWER_ALGEBRAIC)
                        dtop = FLINT_MAX(dtop, gens[idx[k]]);
                if (dtop < 0)
                    ok = 0;

                for (k = 0; k < 5 && ok && st == 0; k++)
                {
                    int sigma;
                    gr_tower_flat_ensure(F);
                    if (idx[k] < 0)
                        st = _dl_special_value(E, -1 - idx[k], F, dtop);
                    else
                    {
                        st = _dl_elementary(E, &sigma, vals + k, words[k], F, dtop);
                        if (st == 0)
                        {
                            fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, gens[idx[k]]), F->mctx);
                            if (sigma < 0)
                                fmpz_mpoly_q_neg(x, x, F->mctx);
                            fmpz_mpoly_q_add(E, E, x, F->mctx);
                        }
                    }
                    if (st == 0)
                    {
                        fmpz_mpoly_q_mul_si(E, E, cf[k], F->mctx);
                        fmpz_mpoly_q_add(expr, expr, E, F->mctx);
                    }
                }

                /* L1 = log((1 - x)/(1 - x y)), L2 = log((1 - y)/(1 - x y)) */
                for (k = 0; k < 2 && ok && st == 0; k++)
                {
                    fmpz_mpoly_q_t u, d;
                    gr_tower_flat_ensure(F);
                    fmpz_mpoly_q_init(u, F->mctx);
                    fmpz_mpoly_q_init(d, F->mctx);
                    fmpz_mpoly_q_mul(d, pts + i, pts + j, F->mctx);
                    fmpz_mpoly_q_sub_si(d, d, 1, F->mctx);
                    fmpz_mpoly_q_sub_si(u, (k == 0) ? pts + i : pts + j, 1, F->mctx);
                    fmpz_mpoly_q_div(u, u, d, F->mctx);
                    if (gr_tower_flat_reduce(u, F) != GR_SUCCESS)
                        st = -1;
                    else
                        st = _gr_tower_special_trans_gen(&dl, F, GR_TOWER_LOG, u, dtop);
                    if (st == 0)
                        fmpz_mpoly_q_gen((k == 0) ? L1 : L2, GR_TOWER_FLAT_VAR_D(F, dl), F->mctx);
                    fmpz_mpoly_q_clear(u, F->mctx);
                    fmpz_mpoly_q_clear(d, F->mctx);
                }

                if (ok && st == 0)
                {
                    fmpz_mpoly_q_mul(x, L1, L2, F->mctx);
                    fmpz_mpoly_q_sub(expr, expr, x, F->mctx);
                    changed = _dl_solve_eliminate(expr, F, gens, ngens);
                }
                else if (st > 0)
                    changed = 1;

                gr_tower_flat_ensure(F);
                fmpz_mpoly_q_clear(expr, F->mctx);
                fmpz_mpoly_q_clear(E, F->mctx);
                fmpz_mpoly_q_clear(x, F->mctx);
                fmpz_mpoly_q_clear(L1, F->mctx);
                fmpz_mpoly_q_clear(L2, F->mctx);
            }

            gr_tower_flat_ensure(F);
            for (k = 0; k < 5; k++)
                fmpz_mpoly_q_clear(vals + k, F->mctx);
        }
    }

    gr_tower_flat_ensure(F);
    for (i = 0; i < npts; i++)
        fmpz_mpoly_q_clear(pts + i, F->mctx);
    flint_free(pts);
    _acb_vec_clear(apts, 2 * ngens);
    return changed;
}

int
_gr_tower_dilog_round(gr_tower_flat_t F, slong limit, slong depth)
{
    gr_tower_struct * T = F->T;
    slong * gens, ngens = 0, j, n;
    int changed = 0;

    (void) depth;

    gens = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));
    for (j = 0; j < limit; j++)
        if (_dl_gen(T->gens + j) || (T->gens[j].kind == GR_TOWER_ALGEBRAIC &&
                T->gens[j].def_kind == GR_TOWER_POLYLOG && T->gens[j].def_param == 2 && T->gens[j].arg.mctx != NULL))
            gens[ngens++] = j;

    /* (one generator suffices: Li_2(y^2) = 2 Li_2(y) + 2 Li_2(-y) with
       y = (sqrt 5 - 1)/2 relates three points of the orbit of (3 - sqrt 5)/2
       and determines that value) */
    if (ngens < 1)
    {
        flint_free(gens);
        return 0;
    }

    /* y: the orbit points with |y| <= 1 of each generator's argument */
    for (j = ngens - 1; j >= 0 && !changed; j--)
    {
        const gr_tower_gen_struct * g = T->gens + gens[j];
        const char * words[6];
        fmpz_mpoly_q_t x, y;
        int nw, i;

        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(x, F->mctx);
        fmpz_mpoly_q_init(y, F->mctx);
        gr_tower_flat_convert(x, &g->arg.data, g->arg.mctx, F);
        {
            /* all six orbit points (the words for nonreal arguments; for
               real ones the subset that stays real) */
            static const char * all[6] = { "", "s", "t", "st", "ts", "r" };
            nw = 6;
            for (i = 0; i < 6; i++)
                words[i] = all[i];
        }

        for (i = 0; i < nw && !changed; i++)
        {
            const char * p;
            int ok = 1;
            acb_t a;
            fmpz_mpoly_q_set(y, x, F->mctx);
            for (p = words[i]; *p && ok; p++)
                ok = _dl_point(y, *p, y, F);
            if (!ok)
                continue;
            /* |y| <= 1 (numerically clear; the boundary is included only
               when exact) */
            acb_init(a);
            ok = (gr_tower_flat_get_acb(a, y, 128, F) == GR_SUCCESS);
            if (ok)
            {
                arb_t m;
                arb_init(m);
                acb_abs(m, a, 128);
                arb_sub_ui(m, m, 1, 128);
                ok = arb_is_nonpositive(m) || arb_is_negative(m);
                arb_clear(m);
            }
            acb_clear(a);
            if (!ok)
                continue;

            for (n = 2; n <= DILOG_MAX_N && !changed; n++)
                changed = _dl_relation(F, y, n, gens, ngens);
        }

        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_clear(x, F->mctx);
        fmpz_mpoly_q_clear(y, F->mctx);
    }

    if (!changed && ngens >= 2)
        changed = _dl_five(F, gens, ngens);

    flint_free(gens);
    return changed;
}

/* -------------------------------------------------------------------- */
/* higher polylogarithms                                                 */
/* -------------------------------------------------------------------- */

/*
    The distribution relations of Li_s for s >= 3,

        Li_s(y^n) = n^(s-1) sum_{k<n} Li_s(zeta_n^k y)       (|y| <= 1),

    between generators: the lazy field keeps them at |x| <= 1 (inversion
    formula, lazy_special.c), so that the relation is linear in the
    generators with no elementary part. (The other functional equations,
    Landen's for Li_3 say, are not used.)
    The candidate y runs over the generators' arguments and n over 2, 3,
    4 (zeta_n present); when every value is a generator (or 0), the latest
    generator with a nonzero coefficient is eliminated (degree one). The
    result is checked numerically.
*/

#define POLYLOG_MAX_N 4

static int
_pl_gen(const gr_tower_gen_struct * g, slong s)
{
    return (g->kind == GR_TOWER_POLYLOG || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_POLYLOG)) &&
        g->def_param == s && g->arg.mctx != NULL;
}

/* the generator (index into gens) with argument v exactly; -2 for v = 0;
   -1 if none */
static slong
_pl_find(const fmpz_mpoly_q_t v, gr_tower_flat_t F, const slong * gens, slong ngens)
{
    fmpz_mpoly_q_t x;
    acb_t a, b;
    slong j, res = -1;

    if (fmpz_mpoly_q_is_zero(v, F->mctx))
        return -2;

    fmpz_mpoly_q_init(x, F->mctx);
    acb_init(a);
    acb_init(b);
    if (gr_tower_flat_get_acb(a, v, 128, F) == GR_SUCCESS)
    {
        for (j = 0; j < ngens && res < 0; j++)
        {
            const gr_tower_gen_struct * g = F->T->gens + gens[j];
            gr_tower_flat_convert(x, &g->arg.data, g->arg.mctx, F);
            if (gr_tower_flat_get_acb(b, x, 128, F) != GR_SUCCESS || !acb_overlaps(a, b))
                continue;
            fmpz_mpoly_q_sub(x, x, v, F->mctx);
            if (gr_tower_flat_reduce(x, F) == GR_SUCCESS && fmpz_mpoly_q_is_zero(x, F->mctx))
                res = j;
        }
    }
    fmpz_mpoly_q_clear(x, F->mctx);
    acb_clear(a);
    acb_clear(b);
    return res;
}

/* the relation of order s for y = the argument of gens[iy] and n */
static int
_pl_relation(gr_tower_flat_t F, slong s, slong iy, slong n, const slong * gens, slong ngens)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_q_t y, z, v, expr, x;
    fmpz * coef;
    fmpz_t npow;
    slong k, j, idx, top = -1, dtop = -1;
    int ok = 1, changed = 0;

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(y, F->mctx);
    fmpz_mpoly_q_init(z, F->mctx);
    fmpz_mpoly_q_init(v, F->mctx);
    fmpz_mpoly_q_init(expr, F->mctx);
    fmpz_mpoly_q_init(x, F->mctx);
    coef = _fmpz_vec_init(ngens);
    fmpz_init(npow);

    gr_tower_flat_convert(y, &T->gens[gens[iy]].arg.data, T->gens[gens[iy]].arg.mctx, F);

    /* zeta_n present (not adjoined: the values must be generators) */
    ok = _gr_tower_special_flat_zeta(z, n, F);

    /* Li_s(y^n) */
    if (ok)
    {
        fmpz_mpoly_q_one(v, F->mctx);
        for (k = 0; k < n; k++)
            fmpz_mpoly_q_mul(v, v, y, F->mctx);
        ok = (gr_tower_flat_reduce(v, F) == GR_SUCCESS);
    }
    if (ok)
    {
        idx = _pl_find(v, F, gens, ngens);
        if (idx == -1)
            ok = 0;
        else if (idx >= 0)
            fmpz_add_ui(coef + idx, coef + idx, 1);
    }

    /* - n^(s-1) Li_s(zeta^k y) */
    fmpz_set_ui(npow, n);
    fmpz_pow_ui(npow, npow, s - 1);
    fmpz_mpoly_q_set(v, y, F->mctx);
    for (k = 0; k < n && ok; k++)
    {
        idx = _pl_find(v, F, gens, ngens);
        if (idx == -1)
            ok = 0;
        else if (idx >= 0)
            fmpz_sub(coef + idx, coef + idx, npow);
        fmpz_mpoly_q_mul(v, v, z, F->mctx);
        ok = ok && (gr_tower_flat_reduce(v, F) == GR_SUCCESS);
    }

    /* the latest transcendental generator, solved for */
    for (j = 0; j < ngens && ok; j++)
        if (!fmpz_is_zero(coef + j) && T->gens[gens[j]].kind != GR_TOWER_ALGEBRAIC && gens[j] > dtop)
        {
            dtop = gens[j];
            top = j;
        }
    if (top < 0)
        ok = 0;

    if (ok)
    {
        fmpq_t c;
        fmpq_init(c);
        fmpz_mpoly_q_zero(expr, F->mctx);
        for (j = 0; j < ngens && ok; j++)
        {
            if (j == top || fmpz_is_zero(coef + j))
                continue;
            if (gens[j] >= dtop)
            {
                ok = 0;
                break;
            }
            fmpq_set_fmpz_frac(c, coef + j, coef + top);
            fmpq_neg(c, c);
            fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, gens[j]), F->mctx);
            fmpz_mpoly_q_mul_fmpq(x, x, c, F->mctx);
            fmpz_mpoly_q_add(expr, expr, x, F->mctx);
        }
        fmpq_clear(c);
        if (ok)
            ok = (gr_tower_flat_reduce(expr, F) == GR_SUCCESS) && _gr_tower_flat_max_dep(expr, F) < dtop;
    }

    /* safeguard: the value of expr is that of the generator */
    if (ok)
        ok = _gr_tower_flat_matches_gen(expr, dtop, F);

    if (ok && _gr_tower_flat_eliminate_gen(dtop, expr, F))
        changed = 1;

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_clear(y, F->mctx);
    fmpz_mpoly_q_clear(z, F->mctx);
    fmpz_mpoly_q_clear(v, F->mctx);
    fmpz_mpoly_q_clear(expr, F->mctx);
    fmpz_mpoly_q_clear(x, F->mctx);
    _fmpz_vec_clear(coef, ngens);
    fmpz_clear(npow);
    return changed;
}

int
_gr_tower_polylog_round(gr_tower_flat_t F, slong limit, slong depth)
{
    gr_tower_struct * T = F->T;
    slong * gens, ngens, j, i, n, s;
    int changed = 0;

    (void) depth;

    gens = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));

    /* each order s >= 3 present (the latest generator of that order) */
    for (i = limit - 1; i >= 0 && !changed; i--)
    {
        int seen = 0;
        s = T->gens[i].def_param;
        if (!_pl_gen(T->gens + i, s) || s < 3)
            continue;
        for (j = i + 1; j < limit && !seen; j++)
            if (_pl_gen(T->gens + j, s))
                seen = 1;
        if (seen)
            continue;

        ngens = 0;
        for (j = 0; j < limit; j++)
            if (_pl_gen(T->gens + j, s))
                gens[ngens++] = j;
        if (ngens < 2)
            continue;

        for (j = ngens - 1; j >= 0 && !changed; j--)
            for (n = 2; n <= POLYLOG_MAX_N && !changed; n++)
                changed = _pl_relation(F, s, j, n, gens, ngens);
    }

    flint_free(gens);
    return changed;
}

POP_OPTIONS
