/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

/*
    Factorization of polynomials over the fields of a tower by Trager's
    method. The norm N(M) of a polynomial M over F_k is the iterated
    resultant of M with the moduli of the steps, a polynomial over the
    base field F_0 of degree deg(M) [F_k : F_0]. If the norm of
    M(x) = h(x - theta) is squarefree for a shift theta (a small integer
    combination of the generators), the irreducible factors of h over F_k
    are gcd(M, N_i)(x + theta) for the irreducible factors N_i of the
    norm over F_0. Polynomials over F_0 are factored by fmpz_poly_factor
    (F_0 = QQ) or fmpz_mpoly_factor (F_0 = QQ(t_1, ..., t_r), the norm
    with denominators cleared being a polynomial in r + 1 variables).

    The steps of the tower up to F_k must be proven irreducible, so that
    F_k is a field; the norms are computed exactly, so the factorization
    is exact (not a numerical heuristic). It fails with GR_UNABLE when no
    shift among the GR_TOWER_FACTOR_SHIFTS tried makes the norm
    squarefree.
*/

/* the norm of M over F_j down to the base field F_0 */
int
_gr_tower_poly_norm(gr_poly_t res, const gr_poly_t M, slong j, gr_tower_t T)
{
    gr_poly_t cur, next;
    slong i;
    int status = GR_SUCCESS;

    gr_poly_init(cur, gr_tower_field_at(T, j));
    status |= gr_poly_set(cur, M, gr_tower_field_at(T, j));

    /* invariant: cur is a polynomial over F_i = F_{i-1}[a_i] / (m_i) */
    for (i = j; i >= 1 && status == GR_SUCCESS; i--)
    {
        gr_poly_init(next, gr_tower_field_at(T, i - 1));
        status = gr_poly_quotient_norm_poly(next, cur, gr_tower_field_at(T, i));
        gr_poly_clear(cur, gr_tower_field_at(T, i));
        gr_poly_init(cur, gr_tower_field_at(T, i - 1));
        gr_poly_swap(cur, next, gr_tower_field_at(T, i - 1));
        gr_poly_clear(next, gr_tower_field_at(T, i - 1));
    }

    if (status == GR_SUCCESS)
        gr_poly_swap(res, cur, T->base);
    gr_poly_clear(cur, gr_tower_field_at(T, i));
    return status;
}

/* theta = sum_{j <= k} s_j a_j in F_k */
static int
_theta(gr_ptr theta, const slong * s, slong k, gr_tower_t T)
{
    gr_ctx_struct * F = gr_tower_field_at(T, k);
    gr_ptr t;
    slong j;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, F);
    status |= gr_zero(theta, F);
    for (j = 1; j <= k; j++)
    {
        if (s[j - 1] != 0)
        {
            gr_ctx_struct * Fj = gr_tower_field_at(T, j);
            gr_ptr a;
            GR_TMP_INIT(a, Fj);
            status |= gr_gen(a, Fj);
            status |= gr_tower_promote(t, a, j, k, T);
            status |= gr_mul_si(t, t, s[j - 1], F);
            status |= gr_add(theta, theta, t, F);
            GR_TMP_CLEAR(a, Fj);
        }
    }
    GR_TMP_CLEAR(t, F);
    return status;
}

/* the polynomial f over F_j as a polynomial over F_k (j <= k) */
static int
_promote_poly(gr_poly_t res, const gr_poly_t f, slong j, slong k, gr_tower_t T)
{
    gr_ctx_struct * Fk = gr_tower_field_at(T, k);
    slong i;
    int status = GR_SUCCESS;

    gr_poly_fit_length(res, f->length, Fk);
    for (i = 0; i < f->length && status == GR_SUCCESS; i++)
        status |= gr_tower_promote(gr_poly_coeff_ptr(res, i, Fk), gr_poly_coeff_srcptr(f, i, gr_tower_field_at(T, j)), j, k, T);
    _gr_poly_set_length(res, f->length, Fk);
    _gr_poly_normalise(res, Fk);
    return status;
}

static int
_append(gr_vec_t fac, fmpz_vec_t mult, gr_poly_t f, ulong e, gr_ctx_t pctx)
{
    slong n = fac->length;
    gr_vec_set_length(fac, n + 1, pctx);
    gr_swap(gr_vec_entry_ptr(fac, n, pctx), f, pctx);
    if (mult != NULL)
        fmpz_vec_append_ui(mult, e);
    return GR_SUCCESS;
}

/*
    The factorization of a nonzero polynomial N over the base field of
    the tower into monic irreducible factors (the constant factor is
    dropped), with multiplicities.
*/
int
_gr_tower_base_poly_factor(gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t N, gr_tower_t T)
{
    gr_ctx_struct * B = T->base;
    gr_ctx_t pctx;
    gr_poly_t f;
    slong i;
    int status = GR_SUCCESS;

    gr_ctx_init_gr_poly(pctx, B);
    gr_vec_set_length(fac, 0, pctx);
    fmpz_vec_set_length(mult, 0);
    gr_poly_init(f, B);

    if (N->length == 0)
        status = GR_DOMAIN;
    else if (N->length == 1)
    {
    }
    else if (B->which_ring == GR_CTX_FMPQ)
    {
        fmpq_poly_t q;
        fmpz_poly_t p;
        fmpz_poly_factor_t pf;

        fmpq_poly_init(q);
        fmpz_poly_init(p);
        fmpz_poly_factor_init(pf);

        for (i = 0; i < N->length; i++)
            fmpq_poly_set_coeff_fmpq(q, i, (const fmpq *) gr_poly_coeff_srcptr(N, i, B));
        fmpq_poly_get_numerator(p, q);
        fmpz_poly_factor(pf, p);

        for (i = 0; i < pf->num && status == GR_SUCCESS; i++)
        {
            status |= gr_poly_set_fmpz_poly(f, pf->p + i, B);
            status |= gr_poly_make_monic(f, f, B);
            status |= _append(fac, mult, f, pf->exp[i], pctx);
        }

        fmpz_poly_factor_clear(pf);
        fmpz_poly_clear(p);
        fmpq_poly_clear(q);
    }
    else if (B->which_ring == GR_CTX_FMPZ_MPOLY_Q)
    {
        fmpz_mpoly_ctx_struct * mctx = gr_ctx_fmpz_mpoly_q_mctx(B);
        slong r = mctx->minfo->nvars, v;
        fmpz_mpoly_ctx_t mctx2;
        fmpz_mpoly_t L, g, t, P, xi;
        fmpz_mpoly_factor_t pf;
        slong * cto, * cback;

        fmpz_mpoly_ctx_init(mctx2, r + 1, ORD_LEX);
        fmpz_mpoly_init(L, mctx);
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

        /* L = lcm of the denominators */
        fmpz_mpoly_one(L, mctx);
        for (i = 0; i < N->length; i++)
        {
            const fmpz_mpoly_q_struct * c = gr_poly_coeff_srcptr(N, i, B);
            if (fmpz_mpoly_q_is_zero(c, mctx))
                continue;
            fmpz_mpoly_gcd(g, L, fmpz_mpoly_q_denref(c), mctx);
            fmpz_mpoly_mul(t, L, fmpz_mpoly_q_denref(c), mctx);
            fmpz_mpoly_divides(L, t, g, mctx);
        }

        /* P = L N as a polynomial in t_1, ..., t_r, x */
        fmpz_mpoly_zero(P, mctx2);
        for (i = 0; i < N->length; i++)
        {
            const fmpz_mpoly_q_struct * c = gr_poly_coeff_srcptr(N, i, B);
            fmpz_mpoly_t u;
            if (fmpz_mpoly_q_is_zero(c, mctx))
                continue;
            fmpz_mpoly_init(u, mctx2);
            fmpz_mpoly_divides(t, L, fmpz_mpoly_q_denref(c), mctx);
            fmpz_mpoly_mul(t, t, fmpz_mpoly_q_numref(c), mctx);
            fmpz_mpoly_compose_fmpz_mpoly_gen(u, t, cto, mctx, mctx2);
            fmpz_mpoly_gen(xi, r, mctx2);
            fmpz_mpoly_pow_ui(xi, xi, i, mctx2);
            fmpz_mpoly_mul(u, u, xi, mctx2);
            fmpz_mpoly_add(P, P, u, mctx2);
            fmpz_mpoly_clear(u, mctx2);
        }

        if (!fmpz_mpoly_factor(pf, P, mctx2))
            status = GR_UNABLE;

        for (i = 0; i < pf->num && status == GR_SUCCESS; i++)
        {
            slong d = fmpz_mpoly_degree_si(pf->poly + i, r, mctx2), e;
            fmpz_mpoly_t ce;

            if (d < 1)
                continue;   /* (content in the t_j: a unit) */

            fmpz_mpoly_init(ce, mctx2);
            gr_poly_fit_length(f, d + 1, B);
            for (e = 0; e <= d; e++)
            {
                slong var = r;
                ulong ee = e;
                fmpz_mpoly_q_struct * fe = gr_poly_coeff_ptr(f, e, B);
                fmpz_mpoly_get_coeff_vars_ui(ce, pf->poly + i, &var, &ee, 1, mctx2);
                fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_numref(fe), ce, cback, mctx2, mctx);
                fmpz_mpoly_one(fmpz_mpoly_q_denref(fe), mctx);
            }
            _gr_poly_set_length(f, d + 1, B);
            _gr_poly_normalise(f, B);
            fmpz_mpoly_clear(ce, mctx2);

            status |= gr_poly_make_monic(f, f, B);
            status |= _append(fac, mult, f, fmpz_get_ui(pf->exp + i), pctx);
        }

        flint_free(cto);
        flint_free(cback);
        fmpz_mpoly_factor_clear(pf, mctx2);
        fmpz_mpoly_clear(xi, mctx2);
        fmpz_mpoly_clear(P, mctx2);
        fmpz_mpoly_clear(t, mctx);
        fmpz_mpoly_clear(g, mctx);
        fmpz_mpoly_clear(L, mctx);
        fmpz_mpoly_ctx_clear(mctx2);
    }
    else
        status = GR_UNABLE;

    gr_poly_clear(f, B);
    gr_ctx_clear(pctx);
    return status;
}

/* whether the norm (over the base field) is certainly not squarefree;
   a cheap test only over QQ, where the factorization is not needed */
static int
_base_not_squarefree(const gr_poly_t N, gr_tower_t T)
{
    int res = 0;

    if (GR_TOWER_BASE_IS_CONSTS(T) && GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL))
    {
        fmpq_poly_t q;
        fmpz_poly_t p;
        slong i;
        fmpq_poly_init(q);
        fmpz_poly_init(p);
        for (i = 0; i < N->length; i++)
            fmpq_poly_set_coeff_fmpq(q, i, (const fmpq *) gr_poly_coeff_srcptr(N, i, T->base));
        fmpq_poly_get_numerator(p, q);
        res = !fmpz_poly_is_squarefree(p);
        fmpz_poly_clear(p);
        fmpq_poly_clear(q);
    }

    return res;
}

/* the shifts tried: 0, then small combinations of the generators with
   growing coefficients (the norm is squarefree for all but finitely
   many shifts) */
static int
_shift(slong * s, slong k, slong attempt)
{
    slong j, c;
    int nonzero = 0;

    if (attempt == 0)
    {
        for (j = 0; j < k; j++)
            s[j] = 0;
        return 1;
    }

    if (k == 0)
        return 0;

    c = 1 + attempt / 4;
    for (j = 0; j < k; j++)
    {
        s[j] = ((attempt * (2 * j + 3) + j * j) % (2 * c + 1)) - c;
        nonzero |= (s[j] != 0);
    }
    if (!nonzero)
        s[k - 1] = 1;
    return 1;
}

/* monic gcd over the field F */
static int
_gcd_monic(gr_poly_t g, const gr_poly_t a, const gr_poly_t b, gr_ctx_t F)
{
    int status = gr_poly_gcd_subresultant(g, a, b, F);
    if (status == GR_SUCCESS && g->length > 0)
        status = gr_poly_make_monic(g, g, F);
    return status;
}

#define GR_TOWER_FACTOR_SHIFTS 12

/*
    The irreducible factors (monic) over F_k of the monic squarefree
    polynomial h over F_k, k >= 1, whose steps are proven. Returns
    GR_UNABLE if no shift making the norm squarefree was found.
*/
int
_gr_tower_poly_factor_squarefree_trager(gr_vec_t fac, const gr_poly_t h, slong k, gr_tower_t T)
{
    gr_ctx_struct * F = gr_tower_field_at(T, k);
    gr_ctx_t pctx, bctx;
    slong * s, D = 1, j, attempt, total;
    gr_poly_t M, N, Nf, G, back;
    gr_vec_t nfac;
    fmpz_vec_t nexp;
    gr_ptr theta;
    int status = GR_SUCCESS, done = 0;

    D = gr_tower_degree_at(T, k);

    gr_ctx_init_gr_poly(pctx, F);
    gr_ctx_init_gr_poly(bctx, T->base);
    gr_vec_set_length(fac, 0, pctx);

    if (h->length <= 2)
    {
        gr_vec_set_length(fac, 1, pctx);
        status = gr_poly_set(gr_vec_entry_ptr(fac, 0, pctx), h, F);
        gr_ctx_clear(pctx);
        gr_ctx_clear(bctx);
        return status;
    }

    s = flint_calloc(FLINT_MAX(k, 1), sizeof(slong));
    gr_poly_init(M, F);
    gr_poly_init(N, T->base);
    gr_poly_init(Nf, F);
    gr_poly_init(G, F);
    gr_poly_init(back, F);
    gr_vec_init(nfac, 0, bctx);
    fmpz_vec_init(nexp, 0);
    GR_TMP_INIT(theta, F);

    for (attempt = 0; attempt < GR_TOWER_FACTOR_SHIFTS && !done && status == GR_SUCCESS; attempt++)
    {
        int squarefree = 1;

        if (!_shift(s, k, attempt))
            break;

        /* M(x) = h(x - theta) */
        status = _theta(theta, s, k, T);
        status |= gr_neg(theta, theta, F);
        if (status == GR_SUCCESS)
            status = gr_poly_taylor_shift(M, h, theta, F);
        status |= gr_neg(theta, theta, F);
        if (status == GR_SUCCESS)
            status = _gr_tower_poly_norm(N, M, k, T);
        if (status != GR_SUCCESS)
            break;

        if (N->length - 1 != D * (h->length - 1) || _base_not_squarefree(N, T))
            continue;

        status = _gr_tower_base_poly_factor(nfac, nexp, N, T);
        if (status != GR_SUCCESS)
            break;

        for (j = 0; j < nfac->length; j++)
            if (!fmpz_is_one(nexp->entries + j))
                squarefree = 0;
        if (!squarefree)
            continue;

        done = 1;

        if (nfac->length == 1)
        {
            gr_vec_set_length(fac, 1, pctx);
            status = gr_poly_set(gr_vec_entry_ptr(fac, 0, pctx), h, F);
            break;
        }

        /* the factors gcd(M, N_i)(x + theta), the norm factors taken
           by increasing degree; each factor found is divided out of M,
           and the last one is the cofactor */
        total = 0;
        {
            slong * order = flint_malloc(sizeof(slong) * nfac->length), a, b;

            for (a = 0; a < nfac->length; a++)
                order[a] = a;
            for (a = 1; a < nfac->length; a++)
                for (b = a; b > 0 && ((gr_poly_struct *) gr_vec_entry_ptr(nfac, order[b], bctx))->length <
                                     ((gr_poly_struct *) gr_vec_entry_ptr(nfac, order[b - 1], bctx))->length; b--)
                    FLINT_SWAP(slong, order[b], order[b - 1]);

            for (a = 0; a < nfac->length && status == GR_SUCCESS; a++)
            {
                if (a == nfac->length - 1)
                {
                    status = gr_poly_set(G, M, F);
                }
                else
                {
                    status |= _promote_poly(Nf, gr_vec_entry_ptr(nfac, order[a], bctx), 0, k, T);
                    if (status == GR_SUCCESS)
                        status = _gcd_monic(G, M, Nf, F);
                    if (status == GR_SUCCESS && G->length >= 2)
                        status = gr_poly_div(M, M, G, F);
                }
                if (status == GR_SUCCESS && G->length < 2)
                    status = GR_UNABLE;   /* (should not happen) */
                if (status == GR_SUCCESS)
                    status = gr_poly_taylor_shift(back, G, theta, F);
                if (status == GR_SUCCESS)
                    status = gr_poly_make_monic(back, back, F);
                if (status == GR_SUCCESS)
                {
                    total += back->length - 1;
                    status = _append(fac, NULL, back, 0, pctx);
                }
            }

            flint_free(order);
        }

        if (status == GR_SUCCESS && total != h->length - 1)
            status = GR_UNABLE;
    }

    if (status == GR_SUCCESS && !done)
        status = GR_UNABLE;

    GR_TMP_CLEAR(theta, F);
    gr_vec_clear(nfac, bctx);
    fmpz_vec_clear(nexp);
    gr_poly_clear(M, F);
    gr_poly_clear(N, T->base);
    gr_poly_clear(Nf, F);
    gr_poly_clear(G, F);
    gr_poly_clear(back, F);
    flint_free(s);
    gr_ctx_clear(pctx);
    gr_ctx_clear(bctx);
    return status;
}

/* whether the steps 1..k are proven, attempting proofs (modular, then
   Trager within the degree limit) for those which are not */
static int
_prove_up_to(slong k, slong degree_limit, gr_tower_t T)
{
    slong j;

    for (j = 1; j <= k; j++)
    {
        gr_tower_gen_struct * g = GR_TOWER_STEP(T, j - 1);

        if (g->status == GR_TOWER_STATUS_PROVEN)
            continue;

        if (GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL))
            gr_tower_prove_step_modular(T, j, GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES));
        if (g->status != GR_TOWER_STATUS_PROVEN)
            gr_tower_prove_step_trager(T, j, degree_limit);
        if (g->status != GR_TOWER_STATUS_PROVEN)
            return 0;
    }

    return 1;
}

int
gr_tower_poly_factor_limit(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, slong k, slong degree_limit, gr_tower_t T)
{
    gr_ctx_struct * F;
    gr_ctx_t pctx;
    gr_poly_t g;
    gr_vec_t part;
    gr_poly_vec_t sqf;
    fmpz_vec_t sexp;
    gr_ptr u;
    slong i, j, D = 1;
    int status = GR_SUCCESS;

    if (f->length == 0)
        return GR_DOMAIN;

    if (!_prove_up_to(k, degree_limit, T))
        return GR_UNABLE;

    /* (proving may have refined the steps: F_k is read afterwards) */
    F = gr_tower_field_at(T, k);
    D = gr_tower_degree_at(T, k);

    gr_ctx_init_gr_poly(pctx, F);
    gr_vec_set_length(fac, 0, pctx);
    fmpz_vec_set_length(mult, 0);

    gr_poly_init(g, F);
    gr_poly_vec_init(sqf, 0, F);
    gr_vec_init(part, 0, pctx);
    fmpz_vec_init(sexp, 0);
    GR_TMP_INIT(u, F);

    /* c = lc(f), g = f / c (fresh, reduced coefficients) */
    gr_poly_fit_length(g, f->length, F);
    for (i = 0; i < f->length; i++)
        status |= gr_set(gr_poly_coeff_ptr(g, i, F), gr_poly_coeff_srcptr(f, i, F), F);
    _gr_poly_set_length(g, f->length, F);
    _gr_poly_normalise(g, F);
    if (g->length == 0)
        status = GR_DOMAIN;
    if (status == GR_SUCCESS)
        status = gr_set(c, gr_poly_coeff_srcptr(g, g->length - 1, F), F);
    if (status == GR_SUCCESS)
        status = gr_poly_make_monic(g, g, F);

    if (status != GR_SUCCESS || g->length <= 1)
        goto cleanup;

    if (k == 0)
    {
        status = _gr_tower_base_poly_factor(fac, mult, g, T);
        goto cleanup;
    }

    status = gr_poly_factor_squarefree(u, sqf, sexp, g, F);

    for (i = 0; i < sqf->length && status == GR_SUCCESS; i++)
    {
        gr_poly_struct * h = gr_poly_vec_entry_ptr(sqf, i, F);
        ulong e = fmpz_get_ui(sexp->entries + i);

        if (h->length <= 1)
            continue;

        status = gr_poly_make_monic(h, h, F);
        if (status != GR_SUCCESS)
            break;

        if (h->length == 2)
        {
            status = _append(fac, mult, h, e, pctx);
            continue;
        }

        if (D * (h->length - 1) > degree_limit)
        {
            status = GR_UNABLE;
            break;
        }

        status = _gr_tower_poly_factor_squarefree_trager(part, h, k, T);
        for (j = 0; j < part->length && status == GR_SUCCESS; j++)
            status = _append(fac, mult, gr_vec_entry_ptr(part, j, pctx), e, pctx);
    }

cleanup:
    GR_TMP_CLEAR(u, F);
    fmpz_vec_clear(sexp);
    gr_vec_clear(part, pctx);
    gr_poly_vec_clear(sqf, F);
    gr_poly_clear(g, F);
    gr_ctx_clear(pctx);
    return status;
}

int
gr_tower_poly_factor(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, slong k, gr_tower_t T)
{
    return gr_tower_poly_factor_limit(c, fac, mult, f, k, GR_TOWER_OPTION(T, GR_TOWER_OPT_FACTOR_DEGREE_LIMIT), T);
}

/* the roots of f in F_k (the linear factors), with multiplicities */
int
gr_tower_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t f, slong k, gr_tower_t T)
{
    gr_ctx_struct * F;
    gr_ctx_t pctx;
    gr_vec_t fac;
    fmpz_vec_t fexp;
    gr_ptr c;
    slong i;
    int status;

    if (f->length == 0)
        return GR_DOMAIN;

    /* (the field is read after the proofs, which may refine it) */
    if (!_prove_up_to(k, GR_TOWER_OPTION(T, GR_TOWER_OPT_FACTOR_DEGREE_LIMIT), T))
        return GR_UNABLE;

    F = gr_tower_field_at(T, k);
    gr_ctx_init_gr_poly(pctx, F);
    gr_vec_init(fac, 0, pctx);
    fmpz_vec_init(fexp, 0);
    GR_TMP_INIT(c, F);

    gr_vec_set_length(roots, 0, F);
    fmpz_vec_set_length(mult, 0);

    status = gr_tower_poly_factor(c, fac, fexp, f, k, T);

    for (i = 0; i < fac->length && status == GR_SUCCESS; i++)
    {
        gr_poly_struct * h = gr_vec_entry_ptr(fac, i, pctx);
        if (h->length == 2)
        {
            gr_vec_set_length(roots, roots->length + 1, F);
            status |= gr_neg(gr_vec_entry_ptr(roots, roots->length - 1, F), gr_poly_coeff_srcptr(h, 0, F), F);
            fmpz_vec_append(mult, fexp->entries + i);
        }
    }

    GR_TMP_CLEAR(c, F);
    fmpz_vec_clear(fexp);
    gr_vec_clear(fac, pctx);
    gr_ctx_clear(pctx);
    return status;
}
