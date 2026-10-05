/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Relations between Hurwitz zeta (polygamma) generators at rational
    arguments, found by the zero test (called from the relation search in
    certify.c): the additive analogue of the gamma relations
    (gamma_relations.c).

    The values x_k = zeta(s, k/N) (0 < k <= N, x_N = zeta(s)) satisfy the
    reflection and distribution relations (special.c), which involve the
    elementary numbers E_k = pi^s (-1)^(s-1) / (s-1)! P_{s-1}(cot(pi k/N))
    and, for even s, zeta(s) (a rational multiple of pi^s). By a
    conjecture of the Chowla-Milnor type they are all the linear
    relations between these values over the algebraic numbers; they are
    theorems, so they are computed exactly. The generators involved are

        psi^(s-1)(k/N) = (-1)^s (s-1)! zeta(s, k/N)    (polygamma),
        zeta(s)                                         (odd s),
        Catalan's constant G = (zeta(2, 1/4) - pi^2) / 8.

    With columns for all values x_k first, then the generators (latest
    first, linked to their values by one row each), then the constants
    E_k, zeta(s) and pi^s, a reduced row echelon form gives the relations
    between generators and constants in the rows whose pivot is a
    generator; each such relation makes its pivot generator (the latest
    one involved) algebraic of degree one, all of them in one round
    (polylogarithms at roots of unity of large order bring in dozens of
    generators at once): a rational combination of earlier ones plus
    pi^s times an element of a cyclotomic field (cot(pi k/N) =
    i (w^k + 1) / (w^k - 1) with w = exp(2 pi i/N), the inverse written
    as a polynomial in w). pi and the roots of unity are adjoined and
    moved to the front as needed. The relation is checked numerically as
    a safeguard.
*/

#include "ulong_extras.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "bernoulli.h"
#include "arb_fmpz_poly.h"
#include "fmpq_mat.h"
#include "fmpz_mpoly_q.h"
#include "arb.h"
#include "acb.h"
#include "acb_dirichlet.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

/* the generator as a Hurwitz value: sets *s (weight), *p, *q (the value
   is zeta(s, p/q), p/q in (0, 1]; for zeta(s): p = q = 1) and the
   factors alpha, beta with x = alpha g + beta pi^s; returns 1, or 0 if
   the generator is not one of these */
static int
_hurwitz_gen(slong * s, slong * p, slong * q, fmpq_t alpha, fmpq_t beta, const gr_tower_gen_struct * g, slong level_limit, slong weight_limit)
{
    fmpz_t a, b;
    int ok = 0;

    fmpq_zero(beta);

    if (g->kind == GR_TOWER_CONSTANT && g->def_param == GR_TOWER_CONST_CATALAN)
    {
        /* zeta(2, 1/4) = pi^2 + 8 G */
        *s = 2;
        *p = 1;
        *q = 4;
        fmpq_set_si(alpha, 8, 1);
        fmpq_one(beta);
        return 1;
    }

    if ((g->kind != GR_TOWER_POLYGAMMA && g->kind != GR_TOWER_ZETA) || g->arg.mctx == NULL)
        return 0;
    if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(&g->arg.data), g->arg.mctx) ||
        !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(&g->arg.data), g->arg.mctx))
        return 0;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_mpoly_get_fmpz(a, fmpz_mpoly_q_numref(&g->arg.data), g->arg.mctx);
    fmpz_mpoly_get_fmpz(b, fmpz_mpoly_q_denref(&g->arg.data), g->arg.mctx);
    if (fmpz_sgn(b) < 0)
    {
        fmpz_neg(a, a);
        fmpz_neg(b, b);
    }

    if (g->kind == GR_TOWER_ZETA)
    {
        /* zeta(s) at an odd integer s >= 3 */
        if (fmpz_is_one(b) && fmpz_cmp_ui(a, 3) >= 0 && fmpz_cmp_ui(a, weight_limit) <= 0 && fmpz_is_odd(a))
        {
            *s = fmpz_get_si(a);
            *p = *q = 1;
            fmpq_one(alpha);
            ok = 1;
        }
    }
    else if (g->def_param >= 1 && g->def_param + 1 <= weight_limit &&
             fmpz_sgn(a) > 0 && fmpz_cmp(a, b) <= 0 && fmpz_cmp_ui(b, level_limit) <= 0)
    {
        /* psi^(m)(a/b) = (-1)^(m+1) m! zeta(m+1, a/b) */
        fmpz_t f;
        fmpz_init(f);
        *s = g->def_param + 1;
        *p = fmpz_get_si(a);
        *q = fmpz_get_si(b);
        fmpz_fac_ui(f, g->def_param);
        if (g->def_param % 2 == 0)
            fmpz_neg(f, f);
        fmpq_set_fmpz(alpha, f);
        fmpq_inv(alpha, alpha);
        fmpz_clear(f);
        ok = 1;
    }

    fmpz_clear(a);
    fmpz_clear(b);
    return ok;
}

/* a polynomial reduced modulo x^m - 1 (exponents folded), in place */
static void
_fold_cyclic(fmpq_poly_t a, ulong m)
{
    slong i, len = fmpq_poly_length(a);
    if (len <= (slong) m)
        return;
    for (i = m; i < len; i++)
        fmpz_add(a->coeffs + (i % m), a->coeffs + (i % m), a->coeffs + i);
    _fmpq_poly_set_length(a, m);
    _fmpq_poly_normalise(a);
    fmpq_poly_canonicalise(a);
}

/*
    E = sum_k c_k P(cot(pi k/N)) (1 <= k < N/2) in Q(zeta_m) = Q[x]/Phi_m
    (4 | m, N | m), univariate and cheap: cot(pi k/N) = i (u + 1) / (u - 1)
    with u = x^((m/N) k), i = x^(m/4), and 1/(u - 1) = (1/r) sum_{j<r} j u^j
    (u of order r). The computation is done modulo x^m - 1 (the powers of
    x are then monomials), with one reduction modulo Phi_m at the end.
*/
void
_gr_tower_cot_sum_cyclo(fmpq_poly_t E, const fmpq * ck, slong N, const fmpz_poly_t P, ulong m)
{
    fmpq_poly_t c, w, t;
    fmpz_poly_t phi;
    fmpq_t q;
    slong k, j, r, i;

    fmpq_poly_init(c);
    fmpq_poly_init(w);
    fmpq_poly_init(t);
    fmpq_init(q);
    fmpq_poly_zero(E);

    for (k = 1; 2 * k < N; k++)
    {
        ulong e;

        if (fmpq_is_zero(ck + k))
            continue;

        r = N / n_gcd(N, k);
        e = ((m / N) * k) % m;
        /* t = (1/r) sum_{j<r} j u^j, then c = i (u + 1) t: with i u^j =
           x^(m/4 + e j) and i u^(j+1) = x^(m/4 + e (j + 1)) */
        {
            fmpz * v = _fmpz_vec_init(m);
            for (j = 1; j < r; j++)
            {
                fmpz_add_si(v + (m / 4 + e * j) % m, v + (m / 4 + e * j) % m, j);
                fmpz_add_si(v + (m / 4 + e * (j + 1)) % m, v + (m / 4 + e * (j + 1)) % m, j);
            }
            fmpq_poly_fit_length(c, m);
            _fmpz_vec_set(c->coeffs, v, m);
            fmpz_one(c->den);
            _fmpq_poly_set_length(c, m);
            _fmpq_poly_normalise(c);
            _fmpz_vec_clear(v, m);
        }
        fmpq_poly_scalar_div_si(c, c, r);
        /* w = P(c) */
        fmpq_poly_zero(w);
        for (i = fmpz_poly_degree(P); i >= 0; i--)
        {
            fmpq_poly_mul(w, w, c);
            _fold_cyclic(w, m);
            fmpq_set_fmpz(q, P->coeffs + i);
            fmpq_poly_add_fmpq(w, w, q);
        }
        fmpq_poly_scalar_mul_fmpq(w, w, ck + k);
        fmpq_poly_add(E, E, w);
    }

    fmpz_poly_init(phi);
    fmpz_poly_cyclotomic(phi, m);
    fmpq_poly_set_fmpz_poly(t, phi);
    fmpq_poly_rem(E, E, t);
    fmpz_poly_clear(phi);

    fmpq_poly_clear(c);
    fmpq_poly_clear(w);
    fmpq_poly_clear(t);
    fmpq_clear(q);
}

/*
    E = sum_k c_k P(cot(pi k/N)) (1 <= k < N/2) as a flat element: the sum
    in Q(zeta_m), then written through zeta_m in the generators, so the
    flat element has no algebraic denominator. zeta_m must be present.
*/
static int
_elementary_flat(fmpz_mpoly_q_t res, const fmpq * ck, slong N, const fmpz_poly_t P, ulong m, gr_tower_flat_t F)
{
    fmpq_poly_t E;
    int ok;

    fmpq_poly_init(E);
    _gr_tower_cot_sum_cyclo(E, ck, N, P, m);
    ok = _gr_tower_special_flat_cyclotomic_eval(res, E, m, F);
    fmpq_poly_clear(E);
    return ok;
}

/* the relations at weight s and level N between the selected
   generators (definition orders sel_d, values x_{sel_k}, factors
   alpha, beta); returns 1 if the tower was changed */
static int
_hurwitz_round_level(gr_tower_flat_t F, slong s, slong N, const slong * sel_d, const slong * sel_k,
    const fmpq * alpha, const fmpq * beta, slong ns)
{
    gr_tower_struct * T = F->T;
    const gr_tower_hurwitz_lattice_struct * L = _gr_tower_hurwitz_lattice(s, N);
    slong nh = N / 2, nb, cg, ce, cz, cp, ncols;
    slong nrows, rank, r, j, i, row = -1, top = -1, nrel = 0, irel;
    slong * bidx, * rows, * tops;
    fmpq_mat_t A, B;
    int changed = 0;

    /* the values x_{k_j} in the normal form of the (cached) lattice at
       level N: the relations x_{k_j} - alpha_j g_j - beta_j pi^s = 0
       become relations between the basis values (columns first, to be
       eliminated), the generators, E, Z and pi^s */
    bidx = flint_malloc(sizeof(slong) * N);
    nb = 0;
    for (j = 0; j < N; j++)
        bidx[j] = (L->row_of_col[j] < 0) ? nb++ : -1;
    cg = nb;
    ce = cg + ns;
    cz = ce + nh;
    cp = cz + 1;
    ncols = cp + 1;

    nrows = ns;
    fmpq_mat_init(A, nrows, ncols);
    for (i = 0; i < ns; i++)
    {
        slong c = L->col_of_k[sel_k[i]], lr = L->row_of_col[c];
        if (lr < 0)
            fmpq_set_si(fmpq_mat_entry(A, i, bidx[c]), 1, 1);
        else
        {
            /* x_c = -sum_{j != c} B[lr, j] (basis values, E, Z) */
            for (j = 0; j < L->ncols; j++)
            {
                const fmpq * v = fmpq_mat_entry(L->B, lr, j);
                fmpq * e;
                if (j == c || fmpq_is_zero(v))
                    continue;
                if (j < N)
                    e = fmpq_mat_entry(A, i, bidx[j]);
                else if (j < GR_TOWER_HURWITZ_CZ(N))
                    e = fmpq_mat_entry(A, i, ce + j - GR_TOWER_HURWITZ_CE(N));
                else
                    e = fmpq_mat_entry(A, i, cz);
                fmpq_sub(e, e, v);
            }
        }
        fmpq_neg(fmpq_mat_entry(A, i, cg + i), alpha + i);
        fmpq_neg(fmpq_mat_entry(A, i, cp), beta + i);
    }
    _gr_tower_hurwitz_lattice_release(L);
    flint_free(bidx);

    fmpq_mat_init(B, nrows, ncols);
    rank = fmpq_mat_rref(B, A);

    /* the rows with their pivots among the generators: independent
       relations, each expressing its pivot generator through later
       columns (earlier generators) and constants; they are all
       eliminated in this round */
    rows = flint_malloc(sizeof(slong) * FLINT_MAX(rank, 1));
    tops = flint_malloc(sizeof(slong) * FLINT_MAX(rank, 1));
    for (r = 0; r < rank; r++)
    {
        for (j = 0; j < ncols; j++)
            if (!fmpq_is_zero(fmpq_mat_entry(B, r, j)))
                break;
        if (j >= cg && j < ce)
        {
            rows[nrel] = r;
            tops[nrel] = j - cg;
            nrel++;
        }
        else if (j >= ce)
            break;
    }

    /* the constants for all of them: pi and zeta_m (moved in front of
       the earliest generator to be eliminated) */
    if (nrel > 0)
    {
        slong dmin = WORD_MAX;
        int any_pi = 0, any_cot = 0, prep;
        for (i = 0; i < nrel; i++)
        {
            dmin = FLINT_MIN(dmin, sel_d[tops[i]]);
            if (!fmpq_is_zero(fmpq_mat_entry(B, rows[i], cp)) || !fmpq_is_zero(fmpq_mat_entry(B, rows[i], cz)))
                any_pi = 1;
            for (j = 1; j <= nh; j++)
            {
                if (!fmpq_is_zero(fmpq_mat_entry(B, rows[i], ce + j - 1)))
                {
                    any_pi = 1;
                    if (2 * j != N)
                        any_cot = 1;
                }
            }
        }
        prep = _gr_tower_special_prepare_constants(T, any_cot ? (N / n_gcd(N, 4)) * 4 : 1, any_pi, dmin);
        if (prep != 0)
            nrel = 0;
        changed = (prep > 0);
    }

    for (irel = 0; irel < nrel; irel++)
    {
        row = rows[irel];
        top = tops[irel];
        {
            /* g_top = sum_j c_j g_j + pi^s (rho + sum_k c_k P(cot(pi k/N))) */
            slong dtop = sel_d[top];
            fmpq * c = _fmpq_vec_init(ns);
            fmpq * ck = _fmpq_vec_init(nh + 1);
            fmpq_t rho, t, bt;
            fmpz_t f;
            int ok = 1, need_cot = 0, need_pi;
            ulong m = 1;

            fmpq_init(rho);
            fmpq_init(t);
            fmpq_init(bt);
            fmpz_init(f);

            /* g_top = -(1/B_top) [sum_{j != top} B_j g_j + B_E E + B_Z Z + B_P pi^s] */
            fmpq_inv(bt, fmpq_mat_entry(B, row, cg + top));
            fmpq_neg(bt, bt);
            for (j = 0; j < ns; j++)
                if (j != top)
                    fmpq_mul(c + j, fmpq_mat_entry(B, row, cg + j), bt);
            fmpq_mul(rho, fmpq_mat_entry(B, row, cp), bt);
            if (s % 2 == 0 && !fmpq_is_zero(fmpq_mat_entry(B, row, cz)))
            {
                _gr_tower_zeta_even_over_pi(t, s);
                fmpq_mul(t, t, fmpq_mat_entry(B, row, cz));
                fmpq_mul(t, t, bt);
                fmpq_add(rho, rho, t);
            }
            /* E_k = pi^s (-1)^(s-1) / (s-1)! P(cot) */
            fmpz_fac_ui(f, s - 1);
            if (s % 2 == 0)
                fmpz_neg(f, f);
            for (j = 1; j <= nh; j++)
            {
                fmpq_mul(ck + j, fmpq_mat_entry(B, row, ce + j - 1), bt);
                fmpq_div_fmpz(ck + j, ck + j, f);
                if (!fmpq_is_zero(ck + j) && 2 * j != N)
                    need_cot = 1;
            }

            /* safeguard: the relation holds numerically */
            {
                acb_t S, v, w, a;
                fmpz_poly_t P;
                slong prec = 128;
                acb_init(S);
                acb_init(v);
                acb_init(w);
                acb_init(a);
                fmpz_poly_init(P);
                _gr_tower_cot_derivative_poly(P, s - 1);

                acb_set_si(w, s);
                for (j = 0; j < ns; j++)
                {
                    /* g_j = (x_j - beta_j pi^s) / alpha_j */
                    fmpq_t cj;
                    fmpq_init(cj);
                    if (j == top)
                        fmpq_set_si(cj, -1, 1);
                    else
                        fmpq_set(cj, c + j);
                    if (!fmpq_is_zero(cj))
                    {
                        acb_set_si(a, sel_k[j]);
                        acb_div_si(a, a, N, prec);
                        acb_hurwitz_zeta(v, w, a, prec);
                        if (!fmpq_is_zero(beta + j))
                        {
                            acb_t pis;
                            acb_init(pis);
                            acb_const_pi(pis, prec);
                            acb_pow_ui(pis, pis, s, prec);
                            acb_mul_fmpz(pis, pis, fmpq_numref(beta + j), prec);
                            acb_div_fmpz(pis, pis, fmpq_denref(beta + j), prec);
                            acb_sub(v, v, pis, prec);
                            acb_clear(pis);
                        }
                        acb_mul_fmpz(v, v, fmpq_denref(alpha + j), prec);
                        acb_div_fmpz(v, v, fmpq_numref(alpha + j), prec);
                        acb_mul_fmpz(v, v, fmpq_numref(cj), prec);
                        acb_div_fmpz(v, v, fmpq_denref(cj), prec);
                        acb_add(S, S, v, prec);
                    }
                    fmpq_clear(cj);
                }
                /* + pi^s (rho + sum ck P(cot)) */
                {
                    arb_t ct, pv, x;
                    acb_t E;
                    arb_init(ct);
                    arb_init(pv);
                    arb_init(x);
                    acb_init(E);
                    acb_set_fmpq(E, rho, prec);
                    for (j = 1; j <= nh; j++)
                    {
                        if (fmpq_is_zero(ck + j))
                            continue;
                        arb_set_si(x, j);
                        arb_div_si(x, x, N, prec);
                        arb_cot_pi(ct, x, prec);
                        arb_fmpz_poly_evaluate_arb(pv, P, ct, prec);
                        arb_mul_fmpz(pv, pv, fmpq_numref(ck + j), prec);
                        arb_div_fmpz(pv, pv, fmpq_denref(ck + j), prec);
                        arb_add(acb_realref(E), acb_realref(E), pv, prec);
                    }
                    acb_const_pi(v, prec);
                    acb_pow_ui(v, v, s, prec);
                    acb_mul(E, E, v, prec);
                    acb_add(S, S, E, prec);
                    arb_clear(ct);
                    arb_clear(pv);
                    arb_clear(x);
                    acb_clear(E);
                }
                ok = acb_contains_zero(S);

                acb_clear(S);
                acb_clear(v);
                acb_clear(w);
                acb_clear(a);
                fmpz_poly_clear(P);
            }

            /* the constants: pi^s, and zeta_m for the cotangents */
            need_pi = !fmpq_is_zero(rho) || need_cot;
            if (need_cot)
                m = (N / n_gcd(N, 4)) * 4;


            if (ok)
            {
                fmpz_mpoly_q_t expr, x, y, E;
                fmpz_poly_t P;

                gr_tower_flat_ensure(F);
                fmpz_mpoly_q_init(expr, F->mctx);
                fmpz_mpoly_q_init(x, F->mctx);
                fmpz_mpoly_q_init(y, F->mctx);
                fmpz_mpoly_q_init(E, F->mctx);
                fmpz_poly_init(P);
                _gr_tower_cot_derivative_poly(P, s - 1);

                fmpz_mpoly_q_zero(expr, F->mctx);
                for (j = 0; j < ns && ok; j++)
                {
                    if (j == top || fmpq_is_zero(c + j))
                        continue;
                    if (sel_d[j] >= dtop)
                    {
                        ok = 0;    /* (not expected: the latest generator is on top) */
                        break;
                    }
                    fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, sel_d[j]), F->mctx);
                    fmpz_mpoly_q_mul_fmpq(x, x, c + j, F->mctx);
                    fmpz_mpoly_q_add(expr, expr, x, F->mctx);
                }

                if (ok && need_pi)
                {
                    /* rho, the term of cot(pi/2) = 0 (the constant term of P),
                       and the cotangents */
                    fmpz_mpoly_q_set_fmpq(E, rho, F->mctx);
                    if (N % 2 == 0 && !fmpq_is_zero(ck + N / 2))
                    {
                        fmpz_poly_get_coeff_fmpz(f, P, 0);
                        fmpq_mul_fmpz(t, ck + N / 2, f);
                        fmpz_mpoly_q_add_fmpq(E, E, t, F->mctx);
                    }
                    if (need_cot)
                    {
                        ok = _elementary_flat(y, ck, N, P, m, F);
                        fmpz_mpoly_q_add(E, E, y, F->mctx);
                    }
                    if (ok)
                    {
                        fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, _gr_tower_certify_find_pi(T, T->num_gens)), F->mctx);
                        _gr_tower_certify_mul_pow_si(E, x, s, F);
                        fmpz_mpoly_q_add(expr, expr, E, F->mctx);
                    }
                }

                if (ok && gr_tower_flat_reduce(expr, F) != GR_SUCCESS)
                    ok = 0;
                if (ok && _gr_tower_flat_eliminate_gen(dtop, expr, F))
                    changed = 1;

                fmpz_mpoly_q_clear(expr, F->mctx);
                fmpz_mpoly_q_clear(x, F->mctx);
                fmpz_mpoly_q_clear(y, F->mctx);
                fmpz_mpoly_q_clear(E, F->mctx);
                fmpz_poly_clear(P);
            }

            _fmpq_vec_clear(c, ns);
            _fmpq_vec_clear(ck, nh + 1);
            fmpq_clear(rho);
            fmpq_clear(t);
            fmpq_clear(bt);
            fmpz_clear(f);
        }
    }

    flint_free(rows);
    flint_free(tops);
    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
    return changed;
}

/*
    The relation round: the Hurwitz generators with definition order
    < limit, grouped by weight; for each generator (latest first), the
    level is the lcm of its denominator and those of the others of the
    same weight which keep it within the limit.
*/
int
_gr_tower_hurwitz_round(gr_tower_flat_t F, slong limit, slong depth)
{
    slong level_limit = GR_TOWER_OPTION(F->T, GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT);
    gr_tower_struct * T = F->T;
    slong * ad, * as, * ap, * aq, * bd, * bk, * tried_s, * tried_N;
    fmpq * aal, * abe, * bal, * bbe;
    slong na = 0, ntried = 0, t, j;
    int changed = 0;

    (void) depth;

    if (limit <= 0)
        return 0;

    ad = flint_malloc(sizeof(slong) * limit);
    as = flint_malloc(sizeof(slong) * limit);
    ap = flint_malloc(sizeof(slong) * limit);
    aq = flint_malloc(sizeof(slong) * limit);
    bd = flint_malloc(sizeof(slong) * limit);
    bk = flint_malloc(sizeof(slong) * limit);
    tried_s = flint_malloc(sizeof(slong) * limit);
    tried_N = flint_malloc(sizeof(slong) * limit);
    aal = _fmpq_vec_init(limit);
    abe = _fmpq_vec_init(limit);
    bal = _fmpq_vec_init(limit);
    bbe = _fmpq_vec_init(limit);

    for (j = limit - 1; j >= 0; j--)
    {
        if (_hurwitz_gen(as + na, ap + na, aq + na, aal + na, abe + na, T->gens + j, level_limit, GR_TOWER_OPTION(T, GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT)))
        {
            ad[na] = j;
            na++;
        }
    }

    for (t = 0; t < na && !changed; t++)
    {
        slong s = as[t], N = aq[t], nb = 0, i;
        int seen = 0, others = 0;

        for (j = 0; j < na; j++)
        {
            slong M;
            if (j == t || as[j] != s)
                continue;
            others = 1;
            M = (N / n_gcd(N, aq[j])) * aq[j];
            if (M <= level_limit)
                N = M;
        }

        /* (a single generator has relations only with the constants:
           at level 1 none; psi^(m)(1/2) and zeta(s) for even s are not
           generators) */
        if (!others && N <= 2)
            continue;

        for (i = 0; i < ntried; i++)
            if (tried_s[i] == s && tried_N[i] == N)
                seen = 1;
        if (seen)
            continue;
        tried_s[ntried] = s;
        tried_N[ntried] = N;
        ntried++;

        for (j = 0; j < na; j++)
        {
            if (as[j] == s && N % aq[j] == 0)
            {
                bd[nb] = ad[j];
                bk[nb] = ap[j] * (N / aq[j]);
                fmpq_set(bal + nb, aal + j);
                fmpq_set(bbe + nb, abe + j);
                nb++;
            }
        }

        changed = _hurwitz_round_level(F, s, N, bd, bk, bal, bbe, nb);
    }

    flint_free(ad); flint_free(as); flint_free(ap); flint_free(aq);
    flint_free(bd); flint_free(bk);
    flint_free(tried_s); flint_free(tried_N);
    _fmpq_vec_clear(aal, limit);
    _fmpq_vec_clear(abe, limit);
    _fmpq_vec_clear(bal, limit);
    _fmpq_vec_clear(bbe, limit);

    /* polygamma generators on rational affine lines (gamma_relations.c) */
    if (!changed)
        changed = _gr_tower_hurwitz_line_round(F, limit, depth);

    return changed;
}
