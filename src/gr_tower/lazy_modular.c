/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Modular forms and functions in the lazy fields.

    A point tau of the upper half-plane is reduced exactly to the
    fundamental domain

        F = { -1/2 < Re(tau) <= 1/2, |tau| > 1 } u { |tau| = 1, 0 <= Re(tau) <= 1/2 }

    by g in PSL(2, Z) (found numerically, then verified and corrected with
    exact comparisons), tau0 = g tau. The values at tau are given by the
    transformation laws with exact multipliers: the 8th roots of unity of
    the theta functions (acb_modular_theta_transform), the 24th roots of
    unity of eta (acb_modular_epsilon_arg), principal square roots of
    i / (c tau + d), powers of c tau + d.

    At tau0, every value is expressed through the generator
    Lambda = lambda(tau0) (GR_TOWER_MODULAR_LAMBDA) and the complete
    elliptic integrals K(Lambda), E(Lambda) (themselves canonical: for tau0
    in F, |Lambda - 1| <= 1, with Im(Lambda) > 0 on the boundary
    Re(tau0) = 1/2, which matches the canonical arguments of K and E):

        theta_3^2 = 2 K / pi,
        theta_2 = theta_3 Lambda^(1/4),  theta_4 = theta_3 (1 - Lambda)^(1/4),
        eta = (theta_2 theta_3 theta_4 / 2)^(1/3),
        j = 256 (1 - Lambda + Lambda^2)^3 / (Lambda (1 - Lambda))^2,
        E_4 = theta_3^8 (1 - Lambda + Lambda^2),
        E_6 = theta_3^12 (1 + Lambda) (2 - Lambda) (1 - 2 Lambda) / 2,
        E_2 = theta_3^4 (3 E / K - 2 + Lambda),

    with principal roots, which are the right branches on F. Thus the
    relations between these values (Jacobi's identity theta_3^4 = theta_2^4
    + theta_4^4, theta_1' = 2 pi eta^3, j in terms of the theta functions,
    the Eisenstein series in terms of each other, ...) hold by
    construction, and the relations of K and E found by the zero test
    (Legendre's relation, which is the transformation tau -> -1/tau, and
    Landen's transformation, which is tau -> 2 tau) apply to the modular
    values as well.

    The theta functions of z are reduced in z as well, modulo the half
    lattice (Z + tau0 Z)/2 with the characteristics, by the parity, and
    by the stabilizer of tau0 at tau0 = i, rho. At the reduced point z0,
    theta_1(z0, tau0) and theta_4(z0, tau0) are generators, and
    theta_2(z0, tau0), theta_3(z0, tau0) algebraic generators of degree 2
    over them (Jacobi's relations); see _mod_theta_z0.
*/

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "fmpq.h"
#include "arith.h"
#include "bernoulli.h"
#include "acb.h"
#include <math.h>
#include "acb_modular.h"
#include "acb_elliptic.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include "gr_tower/impl.h"
#include "gr_tower/lazy_impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

#define FLAGS(ctx) gr_tower_lazy_ctx_field_flags(ctx)
#define ALG(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_ALGEBRAIC) && _gr_tower_lazy_outermost(ctx))
#define REAL(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_REAL) && _gr_tower_lazy_outermost(ctx))

#define CHECK_PREC 64
/* CM points: lambda as an algebraic number for class numbers up to this */
#define CM_CLASS_NUMBER_LIMIT 16
/* ... and discriminants up to this in absolute value */
#define CM_DISC_LIMIT 100000

/* -------------------------------------------------------------------- */
/* exact signs                                                           */
/* -------------------------------------------------------------------- */

/* sign of Re(x) - c */
static int
_mod_re_cmp(int * sgn, gr_srcptr x, slong cnum, slong cden, gr_ctx_t ctx)
{
    fmpq_t c;
    int status;
    fmpq_init(c);
    fmpq_set_si(c, cnum, cden);
    status = _gr_tower_lazy_re_cmp(sgn, x, c, ctx);
    fmpq_clear(c);
    return status;
}

/* sign of |x|^2 - 1 */
static int
_mod_abs_cmp_one(int * sgn, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t, u;
    int status;
    GR_TMP_INIT2(t, u, ctx);
    status = gr_conj(t, x, ctx);
    status |= gr_mul(t, t, x, ctx);
    status |= gr_re(t, t, ctx);
    status |= gr_sub_ui(t, t, 1, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_fast(sgn, t, ctx);
    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* reduction to the fundamental domain                                   */
/* -------------------------------------------------------------------- */

/* res = g tau = (a tau + b) / (c tau + d) */
static int
_mod_apply(gr_ptr res, const psl2z_t g, gr_srcptr tau, gr_ctx_t ctx)
{
    gr_ptr n, d;
    int status;
    GR_TMP_INIT2(n, d, ctx);
    status = gr_mul_fmpz(n, tau, &g->a, ctx);
    status |= gr_add_fmpz(n, n, &g->b, ctx);
    status |= gr_mul_fmpz(d, tau, &g->c, ctx);
    status |= gr_add_fmpz(d, d, &g->d, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(res, n, d, ctx);
    GR_TMP_CLEAR2(n, d, ctx);
    return status;
}

/* res = c tau + d */
static int
_mod_cd(gr_ptr res, const psl2z_t g, gr_srcptr tau, gr_ctx_t ctx)
{
    int status = gr_mul_fmpz(res, tau, &g->c, ctx);
    status |= gr_add_fmpz(res, res, &g->d, ctx);
    return status;
}

/* g := h g, normalized (c >= 0, d > 0 if c = 0) */
static void
_mod_left_mul(psl2z_t g, const psl2z_t h)
{
    psl2z_mul(g, h, g);
}

/*
    tau0 = g tau in the fundamental domain F (see the header), exactly.
    GR_DOMAIN if Im(tau) <= 0.
*/
int
_gr_tower_lazy_modular_reduce(gr_ptr tau0, psl2z_t g, gr_srcptr tau, gr_ctx_t ctx)
{
    acb_t z, w;
    arf_t one_minus_eps;
    psl2z_t h;
    gr_ptr t;
    slong iter;
    int status, sgn;

    /* Im(tau) > 0 */
    GR_TMP_INIT(t, ctx);
    status = gr_im(t, tau, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_fast(&sgn, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    if (status != GR_SUCCESS)
        return status;
    if (sgn <= 0)
        return GR_DOMAIN;

    acb_init(z);
    acb_init(w);
    arf_init(one_minus_eps);
    psl2z_init(h);

    /* a numerical candidate */
    psl2z_one(g);
    {
        slong prec;
        for (prec = CHECK_PREC; prec <= 4096; prec *= 2)
        {
            if (gr_tower_lazy_get_acb(z, tau, prec, ctx) != GR_SUCCESS)
                break;
            if (arb_is_positive(acb_imagref(z)) && acb_rel_accuracy_bits(z) > 20)
            {
                arf_set_d(one_minus_eps, 1.0 - 1.0 / 64);
                acb_modular_fundamental_domain_approx(w, g, z, one_minus_eps, prec);
                if (psl2z_is_correct(g))
                    break;
            }
            psl2z_one(g);
        }
    }

    status = _mod_apply(tau0, g, tau, ctx);

    /* exact corrections */
    for (iter = 0; iter < 64 && status == GR_SUCCESS; iter++)
    {
        int changed = 0;

        /* -1/2 < Re(tau0) <= 1/2: translations */
        status = _mod_re_cmp(&sgn, tau0, 1, 2, ctx);
        if (status == GR_SUCCESS && sgn > 0)
        {
            /* tau0 - 1 */
            psl2z_one(h);
            fmpz_set_si(&h->b, -1);
            _mod_left_mul(g, h);
            status = gr_sub_ui(tau0, tau0, 1, ctx);
            changed = 1;
        }
        else if (status == GR_SUCCESS)
        {
            status = _mod_re_cmp(&sgn, tau0, -1, 2, ctx);
            if (status == GR_SUCCESS && sgn <= 0)
            {
                psl2z_one(h);
                fmpz_set_si(&h->b, 1);
                _mod_left_mul(g, h);
                status = gr_add_ui(tau0, tau0, 1, ctx);
                changed = 1;
            }
        }

        if (status != GR_SUCCESS || changed)
            continue;

        /* |tau0| >= 1, with Re(tau0) >= 0 on the circle */
        status = _mod_abs_cmp_one(&sgn, tau0, ctx);
        if (status == GR_SUCCESS)
        {
            int s2 = 1;
            if (sgn == 0)
                status = _mod_re_cmp(&s2, tau0, 0, 1, ctx);
            if (status == GR_SUCCESS && (sgn < 0 || (sgn == 0 && s2 < 0)))
            {
                /* -1/tau0 */
                psl2z_one(h);
                fmpz_zero(&h->a);
                fmpz_set_si(&h->b, -1);
                fmpz_one(&h->c);
                fmpz_zero(&h->d);
                _mod_left_mul(g, h);
                status = gr_inv(tau0, tau0, ctx);
                status |= gr_neg(tau0, tau0, ctx);
                changed = 1;
            }
        }

        if (status == GR_SUCCESS && !changed)
            break;
    }

    if (status == GR_SUCCESS && iter == 64)
        status = GR_UNABLE;

    acb_clear(z);
    acb_clear(w);
    arf_clear(one_minus_eps);
    psl2z_clear(h);
    return status;
}

/* -------------------------------------------------------------------- */
/* values at a point of the fundamental domain                           */
/* -------------------------------------------------------------------- */

/*
    Whether tau0 (in F) is a CM point, a root of A tau^2 + B tau + C with
    integers A > 0, gcd(A, B, C) = 1 and D = B^2 - 4 A C < 0: then
    Re(tau0) = -B / (2A) and |tau0|^2 = C / A are rational (exactly).
    Returns 1 (setting D, A, B, C), 0, or -1 (unknown).
*/
int
_gr_tower_lazy_modular_cm_form(fmpz_t D, fmpz_t A, fmpz_t B, fmpz_t C, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr t, u;
    fmpq_t x, n;
    fmpz_t g;
    int r1, r2, res = 0;

    GR_TMP_INIT2(t, u, ctx);
    fmpq_init(x);
    fmpq_init(n);
    fmpz_init(g);

    if (gr_re(t, tau0, ctx) != GR_SUCCESS)
        res = -1;
    else
    {
        r1 = _gr_tower_lazy_rational_recognize(x, t, ctx);
        if (r1 == 1)
        {
            if (gr_conj(u, tau0, ctx) == GR_SUCCESS && gr_mul(u, u, tau0, ctx) == GR_SUCCESS)
            {
                r2 = _gr_tower_lazy_rational_recognize(n, u, ctx);
                if (r2 == 1)
                {
                    /* tau^2 - 2 x tau + n = 0 */
                    fmpz_lcm(A, fmpq_denref(x), fmpq_denref(n));
                    fmpz_mul_2exp(B, fmpq_numref(x), 1);
                    fmpz_mul(B, B, A);
                    fmpz_divexact(B, B, fmpq_denref(x));
                    fmpz_neg(B, B);
                    fmpz_mul(C, fmpq_numref(n), A);
                    fmpz_divexact(C, C, fmpq_denref(n));
                    fmpz_gcd3(g, A, B, C);
                    fmpz_divexact(A, A, g);
                    fmpz_divexact(B, B, g);
                    fmpz_divexact(C, C, g);
                    fmpz_mul(D, B, B);
                    fmpz_submul(D, A, C);
                    fmpz_submul(D, A, C);
                    fmpz_submul(D, A, C);
                    fmpz_submul(D, A, C);
                    res = (fmpz_sgn(D) < 0);
                }
                else if (r2 == -1)
                    res = -1;
            }
            else
                res = -1;
        }
        else if (r1 == -1)
            res = -1;
    }

    GR_TMP_CLEAR2(t, u, ctx);
    fmpq_clear(x);
    fmpq_clear(n);
    fmpz_clear(g);
    return res;
}

/*
    lambda(tau0) at a CM point of discriminant D: an algebraic number,
    the root near the enclosure of the factor of
    H_D(256 (1 - x + x^2)^3 / (x (1 - x))^2) (x (1 - x))^(2h)
    (H_D the Hilbert class polynomial, of degree h). Sets *found.
*/
/* the class number h(D) of the discriminant D < 0: the reduced
   primitive forms (a, b, c), |b| <= a <= c, b >= 0 when |b| = a or
   a = c */
static slong
_mod_class_number(slong D)
{
    slong a, b, c, h = 0;

    for (a = 1; 3 * a * a <= -D; a++)
    {
        for (b = -a + 1; b <= a; b++)
        {
            slong n = b * b - D;
            if (n % (4 * a) != 0)
                continue;
            c = n / (4 * a);
            if (c < a)
                continue;
            if (c == a && b < 0)
                continue;
            if (n_gcd(n_gcd(a, FLINT_ABS(b)), c) != 1)
                continue;
            h++;
        }
    }
    return h;
}

static int
_mod_cm_lambda(gr_ptr res, int * found, gr_srcptr tau0, gr_ctx_t ctx)
{
    fmpz_t D, A, B, C;
    int r, status = GR_SUCCESS;

    *found = 0;
    fmpz_init(D);
    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(C);

    r = _gr_tower_lazy_modular_cm_form(D, A, B, C, tau0, ctx);

    if (r == 1 && fmpz_bits(D) <= 30 && FLINT_ABS(fmpz_get_si(D)) <= CM_DISC_LIMIT)
    {
        fmpz_poly_t H, N, M, P, T, U;
        fmpz_poly_factor_t fac;
        acb_t z, t;
        slong h, k, prec;

        fmpz_poly_init(H);
        fmpz_poly_init(N);
        fmpz_poly_init(M);
        fmpz_poly_init(P);
        fmpz_poly_init(T);
        fmpz_poly_init(U);
        fmpz_poly_factor_init(fac);
        acb_init(z);
        acb_init(t);

        /* (the class number first: the class polynomial is costly when
           it is large) */
        h = _mod_class_number(fmpz_get_si(D));
        if (h >= 1 && h <= CM_CLASS_NUMBER_LIMIT)
        {
            acb_modular_hilbert_class_poly(H, fmpz_get_si(D));
            h = fmpz_poly_degree(H);
        }

        if (h >= 1 && h <= CM_CLASS_NUMBER_LIMIT)
        {
            /* N = 256 (1 - x + x^2)^3, M = x^2 (1 - x)^2 */
            fmpz_poly_set_coeff_si(T, 0, 1);
            fmpz_poly_set_coeff_si(T, 1, -1);
            fmpz_poly_set_coeff_si(T, 2, 1);
            fmpz_poly_pow(N, T, 3);
            fmpz_poly_scalar_mul_ui(N, N, 256);
            fmpz_poly_zero(T);
            fmpz_poly_set_coeff_si(T, 1, 1);
            fmpz_poly_set_coeff_si(T, 2, -1);
            fmpz_poly_pow(M, T, 2);

            /* P = sum_k H_k N^k M^(h-k) */
            fmpz_poly_zero(P);
            for (k = 0; k <= h; k++)
            {
                fmpz_poly_pow(T, N, k);
                fmpz_poly_pow(U, M, h - k);
                fmpz_poly_mul(T, T, U);
                fmpz_poly_scalar_mul_fmpz(T, T, H->coeffs + k);
                fmpz_poly_add(P, P, T);
            }

            fmpz_poly_factor(fac, P);

            /* the factor with the root at lambda(tau0) */
            for (prec = 64; prec <= 1024 && !*found && status == GR_SUCCESS; prec *= 2)
            {
                if (gr_tower_lazy_get_acb(t, tau0, prec, ctx) != GR_SUCCESS)
                {
                    status = GR_UNABLE;
                    break;
                }
                acb_modular_lambda(z, t, prec);

                /* the roots of P (all its factors) in z: exactly one */
                {
                    slong count = 0, best = -1, j, d;
                    for (k = 0; k < fac->num && count <= 1; k++)
                    {
                        fmpz_poly_struct * f = fac->p + k;
                        acb_ptr roots;
                        d = fmpz_poly_degree(f);
                        if (d < 1)
                            continue;
                        if (fmpz_sgn(f->coeffs + d) < 0)
                            fmpz_poly_neg(f, f);
                        roots = _acb_vec_init(d);
                        arb_fmpz_poly_complex_roots(roots, f, 0, prec);
                        for (j = 0; j < d; j++)
                            if (acb_overlaps(roots + j, z))
                            {
                                count++;
                                best = k;
                            }
                        _acb_vec_clear(roots, d);
                    }
                    if (count == 1)
                    {
                        qqbar_t q;
                        qqbar_init(q);
                        if (qqbar_set_fmpz_poly_root(q, fac->p + best, z, 4 * prec))
                        {
                            gr_ctx_t QQbar;
                            gr_ctx_init_complex_qqbar(QQbar);
                            status = gr_set_other(res, q, QQbar, ctx);
                            gr_ctx_clear(QQbar);
                            *found = (status == GR_SUCCESS);
                        }
                        qqbar_clear(q);
                    }
                }
            }
        }

        fmpz_poly_clear(H);
        fmpz_poly_clear(N);
        fmpz_poly_clear(M);
        fmpz_poly_clear(P);
        fmpz_poly_clear(T);
        fmpz_poly_clear(U);
        fmpz_poly_factor_clear(fac);
        acb_clear(z);
        acb_clear(t);
    }
    else if (r == -1)
        status = GR_UNABLE;

    fmpz_clear(D);
    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(C);
    return status;
}

/*
    The values at a point tau0 of F from Lambda = lambda(tau0) (given) and
    the elliptic integrals at Lambda: theta_2, theta_3, theta_4, E_2 (any
    of the outputs may be NULL).
*/
static int
_mod_base(gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_ptr e2, gr_srcptr L, gr_ctx_t ctx)
{
    gr_ptr th3, K, u, v;
    fmpq_t e;
    int status;

    GR_TMP_INIT4(th3, K, u, v, ctx);
    fmpq_init(e);

    /* theta_3 = sqrt(2 K(Lambda) / pi) */
    status = gr_elliptic_k(K, L, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_pi(u, ctx);
        status |= gr_div(th3, K, u, ctx);
        status |= gr_mul_ui(th3, th3, 2, ctx);
        status |= gr_sqrt(th3, th3, ctx);
    }

    fmpq_set_si(e, 1, 4);
    if (status == GR_SUCCESS && t2 != NULL)
    {
        status = gr_pow_fmpq(u, L, e, ctx);
        status |= gr_mul(t2, u, th3, ctx);
    }
    if (status == GR_SUCCESS && t4 != NULL)
    {
        status = gr_sub_ui(v, L, 1, ctx);
        status |= gr_neg(v, v, ctx);
        status |= gr_pow_fmpq(u, v, e, ctx);
        status |= gr_mul(t4, u, th3, ctx);
    }
    if (status == GR_SUCCESS && e2 != NULL)
    {
        /* E_2 = theta_3^4 (3 E / K - 2 + Lambda) */
        status = gr_elliptic_e(v, L, ctx);
        if (status == GR_SUCCESS)
            status = gr_div(u, v, K, ctx);
        status |= gr_mul_ui(u, u, 3, ctx);
        status |= gr_sub_ui(u, u, 2, ctx);
        status |= gr_add(u, u, L, ctx);
        status |= gr_pow_ui(v, th3, 4, ctx);
        status |= gr_mul(e2, u, v, ctx);
    }
    if (status == GR_SUCCESS && t3 != NULL)
        status = gr_set(t3, th3, ctx);

    fmpq_clear(e);
    GR_TMP_CLEAR4(th3, K, u, v, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* commensurable points                                                  */
/* -------------------------------------------------------------------- */

/*
    A point tau0 of F which is not a CM point is related to the points
    whose lambda is already a generator (the anchors) when tau0 = gamma
    tau1 with gamma in GL_2^+(Q) of determinant 2^e 3^f, up to
    COMMENSURABLE_LEVEL and COMMENSURABLE_LEVEL_3: its values are then
    computed from those at tau1 (gamma = g M, M = [[a, b], [0, d]], g in
    SL_2(Z)) by the duplication formulas of the theta functions,

        theta_3(2w)^2 = (theta_3^2 + theta_4^2) / 2,  theta_4(2w)^2 = theta_3 theta_4,
        theta_2(2w)^2 = (theta_3^2 - theta_4^2) / 2,
        theta_3(w/2)^2 = theta_3^2 + theta_2^2,  theta_4(w/2)^2 = theta_3^2 - theta_2^2,
        theta_2(w/2)^2 = 2 theta_2 theta_3

    (the roots chosen numerically), with 2 E_2(2w) - E_2(w) = theta_3(2w)^4
    + theta_2(2w)^4, by Jacobi's modular equation of degree 3 for the
    triplings and thirdings (see _mod_trip_step: a root of a quartic,
    adjoined as a new algebraic generator), the translations, and the
    transformation laws of g. Values at commensurable points are thus
    expressed through one generator lambda(tau1) and the elliptic
    integrals at it, and the modular equations of these levels hold by
    construction (lambda(2 tau) = ((1 - k') / (1 + k'))^2, j(2 tau) and
    j(tau) on the modular curve of level 2, the Hauptmoduln of Gamma_0(2)
    and Gamma_0(3), eta quotients, ...).
*/

#define COMMENSURABLE_LEVEL 4     /* 2^4 */
#define COMMENSURABLE_LEVEL_3 2   /* 3^2 */
#define COMMENSURABLE_LEVEL_23 1  /* 2^1 when 3 divides the determinant */

static int _mod_root_of_unity(gr_ptr res, slong r, slong n, gr_ctx_t ctx);
#define COMMENSURABLE_CANDIDATES 64

/* the newest (by definition id) of the generators lambda(tau) on which
   the generator d of F depends through the moduli: the root of the
   chain of steps which gave the value lambda at d (0 if none) */
static ulong
_mod_lambda_root(gr_tower_flat_struct * F, slong d)
{
    gr_tower_struct * T = F->T;
    char * need;
    int * used;
    slong e, v, nvars;
    ulong root = 0;

    gr_tower_flat_ensure(F);
    nvars = F->mctx->minfo->nvars;
    need = flint_calloc(d + 1, 1);
    used = flint_malloc(sizeof(int) * FLINT_MAX(nvars, 1));
    need[d] = 1;
    for (e = d; e >= 0; e--)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(T, e);
        if (!need[e])
            continue;
        if (g->kind == GR_TOWER_MODULAR_LAMBDA)
            root = FLINT_MAX(root, g->def_id);
        else if (g->kind == GR_TOWER_ALGEBRAIC)
        {
            fmpz_mpoly_used_vars(used, _gr_tower_flat_ideal_elem(F, g->index), F->mctx);
            for (v = 0; v < nvars; v++)
            {
                slong f = F->cap - 1 - v;
                if (used[v] && f >= 0 && f < e)
                    need[f] = 1;
            }
        }
    }
    flint_free(need);
    flint_free(used);
    return root;
}

/* the last use of the values coming from the generator lambda of
   definition root (the latest stamp of the points of the cache with
   this root) */
static ulong
_mod_root_last_use(ulong root, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    ulong res = 0;
    slong k;
    for (k = 0; k < L->num_mod_pts; k++)
        if (L->mod_pts[k].root == root)
            res = FLINT_MAX(res, L->mod_pts[k].stamp);
    return res;
}

/* the anchor tau1 and gamma = [[p, q], [r, s]] (gam[0..3]) with tau0 =
   gamma tau1, and the root of its values (see _mod_lambda_root); of the
   anchors found among the first candidates, an SL_2(Z)-equivalent one,
   else those of the most recently used root, and of these the one with
   the cheapest chain of steps (a duplication counts 1, a triplication
   3). (Two roots arise from points commensurable beyond the limits of
   the steps; the values in an evaluation in progress, j(tau) and
   eta(tau) with eta(3 tau), say, then come from the same root, rather
   than from the root of an earlier evaluation at (3 tau + 1)/4 with a
   cheaper chain to 3 tau.) */
static int
_mod_find_anchor(gr_ptr tau1, fmpz * gam, ulong * root_out, int * heavy_out, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * Lz = LAZY(ctx);
    acb_ptr vec;
    acb_t t0;
    fmpz * rel, * cg;
    fmpz_t det, g;
    slong i, d, tried = 0, prec = 256, best = WORD_MAX;
    ulong best_root = 0, best_use = 0, root, use;
    slong best_e3 = 0, cand_e3 = 0;
    ulong seen[COMMENSURABLE_CANDIDATES];
    gr_ptr u, v, c1;

    vec = _acb_vec_init(4);
    acb_init(t0);
    rel = _fmpz_vec_init(4);
    cg = _fmpz_vec_init(4);
    fmpz_init(det);
    fmpz_init(g);
    GR_TMP_INIT3(u, v, c1, ctx);

    if (gr_tower_lazy_get_acb(t0, tau0, prec, ctx) != GR_SUCCESS)
        goto cleanup;

    /* (the latest towers first, each definition once: the generators
       of a long session recur in many towers) */
    for (i = Lz->num_towers - 1; i >= 0 && best > 0 && tried < COMMENSURABLE_CANDIDATES; i--)
    {
        gr_tower_flat_struct * F = Lz->towers[i];
        gr_tower_struct * T = F->T;

        for (d = 0; d < T->num_gens && best > 0 && tried < COMMENSURABLE_CANDIDATES; d++)
        {
            const gr_tower_gen_struct * gen = GR_TOWER_GEN(T, d);
            slong cost, s;
            int ok;

            if (!(gen->kind == GR_TOWER_MODULAR_LAMBDA ||
                  (gen->kind == GR_TOWER_ALGEBRAIC && gen->def_kind == GR_TOWER_MODULAR_LAMBDA)) ||
                gen->arg.mctx == NULL)
                continue;

            for (s = 0; s < tried && seen[s] != gen->def_id; s++) ;
            if (s < tried && gen->def_id != 0)
                continue;
            seen[tried++] = gen->def_id;
            _gr_tower_lazy_set_flat(c1, F, &gen->arg.data, gen->arg.mctx, ctx);

            /* r0 + r1 tau1 + r2 tau0 + r3 tau0 tau1 = 0 */
            acb_one(vec + 0);
            if (gr_tower_lazy_get_acb(vec + 1, c1, prec, ctx) != GR_SUCCESS)
                continue;
            acb_set(vec + 2, t0);
            acb_mul(vec + 3, t0, vec + 1, prec);
            if (!_qqbar_acb_lindep(rel, vec, 4, 1, prec))
                continue;

            /* tau0 = (-r1 tau1 - r0) / (r3 tau1 + r2) */
            fmpz_neg(cg + 0, rel + 1);
            fmpz_neg(cg + 1, rel + 0);
            fmpz_set(cg + 2, rel + 3);
            fmpz_set(cg + 3, rel + 2);
            fmpz_gcd(g, cg + 0, cg + 1);
            fmpz_gcd(g, g, cg + 2);
            fmpz_gcd(g, g, cg + 3);
            if (fmpz_is_zero(g))
                continue;
            _fmpz_vec_scalar_divexact_fmpz(cg, cg, 4, g);
            fmpz_mul(det, cg + 0, cg + 3);
            fmpz_submul(det, cg + 1, cg + 2);
            if (fmpz_sgn(det) < 0)
            {
                _fmpz_vec_neg(cg, cg, 4);
                fmpz_neg(det, det);
            }
            /* (a negative determinant would map the half-plane to the lower one) */
            if (fmpz_sgn(det) <= 0)
                continue;
            {
                fmpz_t r;
                fmpz three = 3;
                slong e2 = fmpz_val2(det), e3;
                fmpz_init(r);
                fmpz_tdiv_q_2exp(r, det, e2);
                e3 = fmpz_remove(r, r, &three);
                ok = fmpz_is_one(r) && e2 <= COMMENSURABLE_LEVEL && e3 <= COMMENSURABLE_LEVEL_3;
                /* (with a tripling or thirding, at most one duplication
                   or halving: the values lie in towers of large degree
                   whose generators are radicals over the roots of the
                   modular equation, and the expressions of j, say, swell
                   beyond that) */
                ok = ok && (e3 == 0 || e2 <= COMMENSURABLE_LEVEL_23);
                /* (nor in the conjugations, lambda(-conj(tau)) for a
                   realness test, say: the new generator is cheaper there) */
                ok = ok && (e3 == 0 || Lz->conj_depth == 0);
                /* (not with triplings followed by thirdings, gamma =
                   [[3 a', b], [0, 3 d']] after reduction, tau1 + 1/3 say:
                   the second modular equation is solved over the field of
                   the first, and the expressions swell) */
                if (ok && e3 >= 2)
                {
                    fmpz_gcd(r, cg + 0, cg + 2);
                    if (fmpz_divisible_si(r, 3))
                    {
                        fmpz_divexact(r, det, r);
                        ok = !fmpz_divisible_si(r, 3);
                    }
                }
                fmpz_clear(r);
                cost = e2 + 3 * e3;
                cand_e3 = e3;
                if (!ok)
                    continue;
                root = _mod_lambda_root(F, d);
                use = (cost == 0) ? 0 : _mod_root_last_use(root, ctx);
                if (!(cost == 0 || best == WORD_MAX || use > best_use ||
                      (use == best_use && cost < best)))
                    continue;
            }

            /* exactly */
            ok = (gr_mul_fmpz(u, c1, cg + 2, ctx) == GR_SUCCESS);
            ok = ok && (gr_add_fmpz(u, u, cg + 3, ctx) == GR_SUCCESS);
            ok = ok && (gr_mul(u, u, tau0, ctx) == GR_SUCCESS);
            ok = ok && (gr_mul_fmpz(v, c1, cg + 0, ctx) == GR_SUCCESS);
            ok = ok && (gr_add_fmpz(v, v, cg + 1, ctx) == GR_SUCCESS);
            if (ok && gr_equal(u, v, ctx) == T_TRUE)
            {
                best = cost;
                best_root = root;
                best_use = use;
                best_e3 = cand_e3;
                gr_swap(tau1, c1, ctx);
                _fmpz_vec_set(gam, cg, 4);
            }

            /* (the towers may have changed) */
            if (i >= Lz->num_towers || Lz->towers[i] != F)
                break;
        }
    }

cleanup:
    _acb_vec_clear(vec, 4);
    acb_clear(t0);
    _fmpz_vec_clear(rel, 4);
    _fmpz_vec_clear(cg, 4);
    fmpz_clear(det);
    fmpz_clear(g);
    GR_TMP_CLEAR3(u, v, c1, ctx);
    *root_out = best_root;
    *heavy_out = (best != WORD_MAX && best_e3 >= 1 && best >= 4);
    return best != WORD_MAX;
}

/* res = +- sqrt(x), either sign: as i sqrt(-x) when x is closer to the
   negative real axis (the principal root of x would have to decide on
   which side of the branch cut x lies, by a realness test, which costs
   a conjugation) */
static int
_mod_sqrt_any(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    acb_t z;
    int status;

    acb_init(z);
    status = gr_tower_lazy_get_acb(z, x, 64, ctx);
    if (status == GR_SUCCESS && arf_sgn(arb_midref(acb_realref(z))) < 0)
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status = gr_neg(t, x, ctx);
        status |= gr_sqrt(t, t, ctx);
        status |= gr_i(res, ctx);
        status |= gr_mul(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);
    }
    else if (status == GR_SUCCESS)
        status = gr_sqrt(res, x, ctx);
    acb_clear(z);
    return status;
}

/*
    res = +-res, the sign for which res overlaps the enclosure ref of the
    value while -res does not (more precision for res if both do, as when
    the enclosure of res contains 0): an exact choice between the two
    exact candidates; GR_UNABLE if undecided.
*/
static int
_mod_fix_sign(gr_ptr res, const acb_t ref, gr_ctx_t ctx)
{
    acb_t z, mz;
    slong prec;
    int status = GR_UNABLE, o1, o2;

    acb_init(z);
    acb_init(mz);
    for (prec = 128; prec <= 1024; prec *= 2)
    {
        if (gr_tower_lazy_get_acb(z, res, prec, ctx) != GR_SUCCESS)
            break;
        acb_neg(mz, z);
        o1 = acb_overlaps(z, ref);
        o2 = acb_overlaps(mz, ref);
        if (o1 && !o2)
        {
            status = GR_SUCCESS;
            break;
        }
        if (o2 && !o1)
        {
            status = gr_neg(res, res, ctx);
            break;
        }
        if (!o1 && !o2)
            break;
    }
    acb_clear(z);
    acb_clear(mz);
    return status;
}

/* res = +- sqrt(x), the sign agreeing with ref */
static int
_mod_sqrt_near(gr_ptr res, gr_srcptr x, const acb_t ref, gr_ctx_t ctx)
{
    int status = _mod_sqrt_any(res, x, ctx);
    if (status == GR_SUCCESS)
        status = _mod_fix_sign(res, ref, ctx);
    return status;
}

/*
    th[0..2] = theta_2, theta_3, theta_4 (and *e2 = E_2, if not NULL) at
    w, times scale, are replaced by the values at 2w (dir = 1) or w/2 (dir
    = -1), times scale; w is updated.
*/
static int
_mod_dup_step(gr_ptr * th, gr_ptr e2, gr_ptr w, int dir, const acb_t scale, gr_ctx_t ctx)
{
    gr_ptr s2, s3, s4, n2, n3, n4, u;
    acb_t aw, r1, r2, r3, r4;
    int status;

    GR_TMP_INIT3(s2, s3, s4, ctx);
    GR_TMP_INIT4(n2, n3, n4, u, ctx);
    acb_init(aw); acb_init(r1); acb_init(r2); acb_init(r3); acb_init(r4);

    status = gr_sqr(s2, th[0], ctx);
    status |= gr_sqr(s3, th[1], ctx);
    status |= gr_sqr(s4, th[2], ctx);

    if (dir > 0)
    {
        status |= gr_add(n3, s3, s4, ctx);
        status |= gr_div_ui(n3, n3, 2, ctx);
        status |= gr_mul(n4, th[1], th[2], ctx);
        status |= gr_sub(n2, s3, s4, ctx);
        status |= gr_div_ui(n2, n2, 2, ctx);
        status |= gr_mul_ui(w, w, 2, ctx);
    }
    else
    {
        status |= gr_add(n3, s3, s2, ctx);
        status |= gr_sub(n4, s3, s2, ctx);
        status |= gr_mul(n2, th[0], th[1], ctx);
        status |= gr_mul_ui(n2, n2, 2, ctx);
        status |= gr_div_ui(w, w, 2, ctx);
    }

    if (status == GR_SUCCESS)
    {
        status = gr_tower_lazy_get_acb(aw, w, 128, ctx);
        if (status == GR_SUCCESS)
        {
            acb_zero(r1);
            acb_modular_theta(r1, r2, r3, r4, r1, aw, 128);
            acb_mul(r2, r2, scale, 128);
            acb_mul(r3, r3, scale, 128);
            acb_mul(r4, r4, scale, 128);
            status = _mod_sqrt_near(th[0], n2, r2, ctx);
            if (status == GR_SUCCESS)
                status = _mod_sqrt_near(th[1], n3, r3, ctx);
            if (status == GR_SUCCESS)
                status = _mod_sqrt_near(th[2], n4, r4, ctx);
        }
    }

    if (status == GR_SUCCESS && e2 != NULL)
    {
        if (dir > 0)
        {
            /* E_2(2w) = (E_2(w) + theta_3(2w)^4 + theta_2(2w)^4) / 2 */
            status = gr_sqr(u, n3, ctx);
            status |= gr_add(e2, e2, u, ctx);
            status |= gr_sqr(u, n2, ctx);
            status |= gr_add(e2, e2, u, ctx);
            status |= gr_div_ui(e2, e2, 2, ctx);
        }
        else
        {
            /* E_2(w/2) = 2 E_2(w) - theta_3(w)^4 - theta_2(w)^4 */
            status = gr_mul_ui(e2, e2, 2, ctx);
            status |= gr_sqr(u, s3, ctx);
            status |= gr_sub(e2, e2, u, ctx);
            status |= gr_sqr(u, s2, ctx);
            status |= gr_sub(e2, e2, u, ctx);
        }
    }

    GR_TMP_CLEAR3(s2, s3, s4, ctx);
    GR_TMP_CLEAR4(n2, n3, n4, u, ctx);
    acb_clear(aw); acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(r4);
    return status;
}

/*
    th[0..2] = theta_2, theta_3, theta_4 (and *e2 = E_2, if not NULL) at
    w, times scale, are replaced by the values at 3w (dir = 1) or w/3 (dir
    = -1), times scale; w is updated. With u = sqrt(theta_2/theta_3)
    at w and v at 3w (with the branches of the acb conventions, up to a
    common sign), Jacobi's modular equation of degree 3 reads
    u^4 - v^4 - 2uv(1 - u^2 v^2) = 0, the multiplier is
    theta_3(w)^2 / theta_3(3w)^2 = 1 + 2 v^3 / u, and
    3 E_2(3w) - E_2(w) = sum_{k=2,3,4} theta_k(w)^2 theta_k(3w)^2.
    The root of the equation is chosen numerically (+-v is known
    numerically, and exactly one root of the quartic matches).
*/
static int
_mod_trip_step(gr_ptr * th, gr_ptr e2, gr_ptr w, int dir, const acb_t scale, gr_ctx_t ctx)
{
    gr_ptr u, v, M, t, s, n2, n3, n4, old[3];
    gr_poly_t p;
    acb_t aw, r1, r2, r3, r4, ref;
    slong k;
    int status;

    GR_TMP_INIT4(u, v, M, t, ctx);
    GR_TMP_INIT4(s, n2, n3, n4, ctx);
    GR_TMP_INIT3(old[0], old[1], old[2], ctx);
    gr_poly_init(p, ctx);
    acb_init(aw); acb_init(r1); acb_init(r2); acb_init(r3); acb_init(r4); acb_init(ref);

    for (k = 0; k < 3; k++)
        GR_MUST_SUCCEED(gr_set(old[k], th[k], ctx));

    /* the known root: at w (dir = 1) or at w = 3 (w/3) (dir = -1) */
    status = gr_div(t, th[0], th[1], ctx);
    if (status == GR_SUCCESS)
        status = _mod_sqrt_any((dir > 0) ? u : v, t, ctx);

    /* the new point, numerically */
    if (dir > 0)
        status |= gr_mul_ui(w, w, 3, ctx);
    else
        status |= gr_div_ui(w, w, 3, ctx);
    if (status == GR_SUCCESS)
        status = gr_tower_lazy_get_acb(aw, w, 128, ctx);
    if (status == GR_SUCCESS)
    {
        acb_zero(r1);
        acb_modular_theta(r1, r2, r3, r4, r1, aw, 128);
        acb_div(ref, r2, r3, 128);
        acb_sqrt(ref, ref, 128);
        acb_div(r4, r4, r3, 128);   /* theta_4 / theta_3 at the new point */
    }

    /* the other root: v from u (X^4 - 2u^3 X^3 + 2u X - u^4), or u from
       v (X^4 + 2v^3 X^3 - 2v X - v^4) */
    if (status == GR_SUCCESS)
    {
        gr_srcptr a = (dir > 0) ? u : v;
        int sg = (dir > 0) ? 1 : -1;
        gr_poly_fit_length(p, 5, ctx);
        status = gr_pow_ui(t, a, 4, ctx);
        status |= gr_neg(gr_poly_coeff_ptr(p, 0, ctx), t, ctx);
        status |= gr_mul_si(gr_poly_coeff_ptr(p, 1, ctx), a, 2 * sg, ctx);
        status |= gr_zero(gr_poly_coeff_ptr(p, 2, ctx), ctx);
        status |= gr_pow_ui(t, a, 3, ctx);
        status |= gr_mul_si(gr_poly_coeff_ptr(p, 3, ctx), t, -2 * sg, ctx);
        status |= gr_one(gr_poly_coeff_ptr(p, 4, ctx), ctx);
        _gr_poly_set_length(p, 5, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_poly_root_near((dir > 0) ? v : u, p, ref, 1, ctx);
    }

    /* M = 1 + 2 v^3 / u; theta_3(3w)^2 = theta_3(w)^2 / M */
    if (status == GR_SUCCESS)
    {
        status = gr_pow_ui(M, v, 3, ctx);
        status |= gr_mul_ui(M, M, 2, ctx);
        status |= gr_div(M, M, u, ctx);
        status |= gr_add_ui(M, M, 1, ctx);
        status |= gr_sqr(n3, th[1], ctx);
        if (dir > 0)
            status |= gr_div(n3, n3, M, ctx);
        else
            status |= gr_mul(n3, n3, M, ctx);
    }
    if (status == GR_SUCCESS)
    {
        acb_mul(r3, r3, scale, 128);
        status = _mod_sqrt_near(th[1], n3, r3, ctx);
    }
    /* theta_2 = (u or v)^2 theta_3 */
    if (status == GR_SUCCESS)
    {
        status = gr_sqr(t, (dir > 0) ? v : u, ctx);
        status |= gr_mul(th[0], t, th[1], ctx);
    }
    /* theta_4 = theta_3 (1 - (theta_2/theta_3)^4)^(1/4) (Jacobi's identity) */
    if (status == GR_SUCCESS)
    {
        acb_t q;
        acb_init(q);
        status = gr_pow_ui(s, t, 4, ctx);
        status |= gr_sub_ui(s, s, 1, ctx);
        status |= gr_neg(s, s, ctx);
        acb_sqr(q, r4, 128);
        if (status == GR_SUCCESS)
            status = _mod_sqrt_near(n4, s, q, ctx);
        if (status == GR_SUCCESS)
            status = _mod_sqrt_near(n2, n4, r4, ctx);
        if (status == GR_SUCCESS)
            status = gr_mul(th[2], n2, th[1], ctx);
        acb_clear(q);
    }

    /* E_2: 3 E_2(3w) - E_2(w) = sum theta_k(w)^2 theta_k(3w)^2 */
    if (status == GR_SUCCESS && e2 != NULL)
    {
        status = gr_zero(s, ctx);
        for (k = 0; k < 3; k++)
        {
            status |= gr_mul(t, old[k], th[k], ctx);
            status |= gr_sqr(t, t, ctx);
            status |= gr_add(s, s, t, ctx);
        }
        if (dir > 0)
        {
            status |= gr_add(e2, e2, s, ctx);
            status |= gr_div_ui(e2, e2, 3, ctx);
        }
        else
        {
            status |= gr_mul_ui(e2, e2, 3, ctx);
            status |= gr_sub(e2, e2, s, ctx);
        }
    }

    GR_TMP_CLEAR4(u, v, M, t, ctx);
    GR_TMP_CLEAR4(s, n2, n3, n4, ctx);
    GR_TMP_CLEAR3(old[0], old[1], old[2], ctx);
    gr_poly_clear(p, ctx);
    acb_clear(aw); acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(r4); acb_clear(ref);
    return status;
}

/* the values at tau0 = gamma tau1 from those at the anchor tau1 */
static int
_mod_from_anchor(gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_ptr e2, const fmpz * gam, gr_srcptr tau1, int ratios, gr_ctx_t ctx)
{
    gr_ptr th[3], E, w, L1, A, u;
    acb_t scale;
    fmpz_t x, y, e, k, a, b, dd, N;
    fmpz three = 3;
    slong k3;
    psl2z_t g;
    slong i;
    int R[4], S[4], C;
    int status;

    GR_TMP_INIT3(th[0], th[1], th[2], ctx);
    GR_TMP_INIT5(E, w, L1, A, u, ctx);
    acb_init(scale);
    acb_one(scale);
    fmpz_init(x); fmpz_init(y); fmpz_init(e); fmpz_init(k);
    fmpz_init(a); fmpz_init(b); fmpz_init(dd); fmpz_init(N);
    psl2z_init(g);

    /* the values at the anchor (with ratios, divided by theta_3(tau1):
       the duplication formulas and the transformation laws are
       homogeneous, and K is not needed) */
    status = _gr_tower_lazy_special_gen_locked(L1, tau1, GR_TOWER_MODULAR_LAMBDA, 0, ctx);
    if (status == GR_SUCCESS && ratios)
    {
        fmpq_t q;
        fmpq_init(q);
        fmpq_set_si(q, 1, 4);
        status = gr_pow_fmpq(th[0], L1, q, ctx);
        status |= gr_one(th[1], ctx);
        status |= gr_sub_ui(th[2], L1, 1, ctx);
        status |= gr_neg(th[2], th[2], ctx);
        status |= gr_pow_fmpq(th[2], th[2], q, ctx);
        fmpq_clear(q);
        e2 = NULL;
        /* scale = 1 / theta_3(tau1) */
        if (status == GR_SUCCESS)
        {
            acb_t z0, r1, r2, r4;
            acb_init(z0); acb_init(r1); acb_init(r2); acb_init(r4);
            status = gr_tower_lazy_get_acb(z0, tau1, 128, ctx);
            if (status == GR_SUCCESS)
            {
                acb_zero(r1);
                acb_modular_theta(r1, r2, scale, r4, r1, z0, 128);
                acb_inv(scale, scale, 128);
            }
            acb_clear(z0); acb_clear(r1); acb_clear(r2); acb_clear(r4);
        }
    }
    else if (status == GR_SUCCESS)
        status = _mod_base(th[0], th[1], th[2], (e2 != NULL) ? E : NULL, L1, ctx);

    /* gamma = g M: g^-1 = [[x, y], [-r/e, p/e]] (x p + y r = e), then
       M = [[e, b], [0, N/e]] with b reduced modulo N/e */
    fmpz_mul(N, gam + 0, gam + 3);
    fmpz_submul(N, gam + 1, gam + 2);
    fmpz_xgcd(e, x, y, gam + 0, gam + 2);
    if (fmpz_sgn(e) < 0)
    {
        fmpz_neg(e, e); fmpz_neg(x, x); fmpz_neg(y, y);
    }
    fmpz_set(a, e);
    fmpz_divexact(dd, N, e);
    fmpz_mul(b, x, gam + 1);
    fmpz_addmul(b, y, gam + 3);
    fmpz_fdiv_qr(k, b, b, dd);
    {
        /* g^-1 = [[1, -k], [0, 1]] [[x, y], [-r/e, p/e]] */
        fmpz_t gi00, gi01, gi10, gi11;
        fmpz_init(gi00); fmpz_init(gi01); fmpz_init(gi10); fmpz_init(gi11);
        fmpz_divexact(gi10, gam + 2, e);
        fmpz_neg(gi10, gi10);
        fmpz_divexact(gi11, gam + 0, e);
        fmpz_set(gi00, x);
        fmpz_submul(gi00, k, gi10);
        fmpz_set(gi01, y);
        fmpz_submul(gi01, k, gi11);
        /* g = (g^-1)^-1 = [[gi11, -gi01], [-gi10, gi00]] */
        fmpz_set(&g->a, gi11);
        fmpz_neg(&g->b, gi01);
        fmpz_neg(&g->c, gi10);
        fmpz_set(&g->d, gi00);
        if (fmpz_sgn(&g->c) < 0 || (fmpz_is_zero(&g->c) && fmpz_sgn(&g->d) < 0))
        {
            fmpz_neg(&g->a, &g->a); fmpz_neg(&g->b, &g->b);
            fmpz_neg(&g->c, &g->c); fmpz_neg(&g->d, &g->d);
        }
        fmpz_clear(gi00); fmpz_clear(gi01); fmpz_clear(gi10); fmpz_clear(gi11);
    }

    /* w = M tau1: doublings and triplings, the translation, halvings and
       thirdings */
    status |= gr_set(w, tau1, ctx);
    for (i = 0; i < (slong) fmpz_val2(a) && status == GR_SUCCESS; i++)
        status = _mod_dup_step(th, (e2 != NULL) ? E : NULL, w, 1, scale, ctx);
    {
        fmpz_t r;
        fmpz_init(r);
        k3 = fmpz_remove(r, a, (const fmpz *) &three);
        for (i = 0; i < k3 && status == GR_SUCCESS; i++)
            status = _mod_trip_step(th, (e2 != NULL) ? E : NULL, w, 1, scale, ctx);
        fmpz_clear(r);
    }
    if (status == GR_SUCCESS && !fmpz_is_zero(b))
    {
        /* theta_2(w + b) = exp(pi i b / 4) theta_2(w); theta_3, theta_4
           exchanged for odd b */
        slong bb = fmpz_fdiv_ui(b, 8);
        status = _mod_root_of_unity(u, bb, 4, ctx);
        status |= gr_mul(th[0], th[0], u, ctx);
        if (bb % 2 == 1)
            gr_swap(th[1], th[2], ctx);
        status |= gr_add_fmpz(w, w, b, ctx);
    }
    for (i = 0; i < (slong) fmpz_val2(dd) && status == GR_SUCCESS; i++)
        status = _mod_dup_step(th, (e2 != NULL) ? E : NULL, w, -1, scale, ctx);
    {
        fmpz_t r;
        fmpz_init(r);
        k3 = fmpz_remove(r, dd, (const fmpz *) &three);
        for (i = 0; i < k3 && status == GR_SUCCESS; i++)
            status = _mod_trip_step(th, (e2 != NULL) ? E : NULL, w, -1, scale, ctx);
        fmpz_clear(r);
    }

    /* tau0 = g w: theta_{1+i}(w) = exp(pi i R_i / 4) A theta_{1+S_i}(tau0),
       A = sqrt(i / (c w + d)) */
    if (status == GR_SUCCESS)
    {
        gr_ptr out[4];
        acb_modular_theta_transform(R, S, &C, g);
        if (C)
        {
            status = _mod_cd(u, g, w, ctx);
            status |= gr_inv(A, u, ctx);
            status |= gr_i(u, ctx);
            status |= gr_mul(A, A, u, ctx);
            status |= gr_sqrt(A, A, ctx);
        }
        else
            status = gr_one(A, ctx);

        out[0] = NULL; out[1] = t2; out[2] = t3; out[3] = t4;
        for (i = 1; i < 4 && status == GR_SUCCESS; i++)
        {
            if (out[S[i]] == NULL)
                continue;
            status = _mod_root_of_unity(u, R[i], 4, ctx);
            status |= gr_mul(u, u, A, ctx);
            status |= gr_div(out[S[i]], th[i - 1], u, ctx);
        }

        /* E_2(g w) = (c w + d)^2 E_2(w) + 6 c (c w + d) / (pi i) */
        if (status == GR_SUCCESS && e2 != NULL)
        {
            status = _mod_cd(u, g, w, ctx);
            status |= gr_sqr(A, u, ctx);
            status |= gr_mul(e2, A, E, ctx);
            if (!fmpz_is_zero(&g->c))
            {
                status |= gr_mul_fmpz(u, u, &g->c, ctx);
                status |= gr_mul_ui(u, u, 6, ctx);
                status |= gr_pi(A, ctx);
                status |= gr_div(u, u, A, ctx);
                status |= gr_i(A, ctx);
                status |= gr_div(u, u, A, ctx);
                status |= gr_add(e2, e2, u, ctx);
            }
        }
    }

    GR_TMP_CLEAR3(th[0], th[1], th[2], ctx);
    GR_TMP_CLEAR5(E, w, L1, A, u, ctx);
    acb_clear(scale);
    fmpz_clear(x); fmpz_clear(y); fmpz_clear(e); fmpz_clear(k);
    fmpz_clear(a); fmpz_clear(b); fmpz_clear(dd); fmpz_clear(N);
    psl2z_clear(g);
    return status;
}

/*
    Linked generators. The chains of steps with a tripling and a further
    step (a halving, another tripling) give values in towers of large
    degree, where everything at tau0 costs (the theta functions of z at
    tau0 = 4i/pi through two thirdings from an earlier point 1/3 + pi i/4:
    3.5 seconds rather than 0.01; j(tau), eta(3 tau) after lambda((2 tau +
    1)/3): minutes). There lambda(tau0) is a new generator, and so are
    the values at tau0 from it (_mod_base: the fourth roots of lambda and
    1 - lambda, K(lambda), theta_3 = sqrt(2 K / pi), E(lambda)), and the
    context records their values through the chain from the anchor
    (rebase records triggered by the root of the anchor): a tower which
    meets values from both makes them the values of the chain, so that
    the identities between the two hold as with the chain; elsewhere the
    new generators cost what they cost in a new context.
*/

/* the chain of a link, shared by its records: computed when a record
   needs it */
typedef struct
{
    slong refs;
    slong visits;       /* (_mod_link_rule_refs: the records counted so far) */
    gr_tower_lazy_elem_struct tau1;
    fmpz gam[4];
    int have;           /* 1: a2, a3, a4, e2 computed; -1: failed */
    gr_tower_lazy_elem_struct a2, a3, a4, e2;
}
mod_link_struct;

#define LINK_LAMBDA 0
#define LINK_R2 1
#define LINK_R4 2
#define LINK_K 3
#define LINK_TH3 4
#define LINK_E 5

typedef struct
{
    mod_link_struct * S;
    int which;
    gr_tower_lazy_elem_struct scale;    /* the generator = scale * the value of x */
}
mod_link_rule_struct;

static void
_mod_link_rule_clear(void * data, gr_tower_lazy_ctx_struct * L)
{
    mod_link_rule_struct * R = data;
    mod_link_struct * S = R->S;
    slong k;
    _gr_tower_lazy_clear_data(&R->scale);
    _gr_tower_lazy_unref(L, R->scale.F);
    if (--S->refs == 0)
    {
        _gr_tower_lazy_clear_data(&S->tau1);
        _gr_tower_lazy_unref(L, S->tau1.F);
        for (k = 0; k < 4; k++)
            fmpz_clear(S->gam + k);
        if (S->have == 1)
        {
            _gr_tower_lazy_clear_data(&S->a2); _gr_tower_lazy_unref(L, S->a2.F);
            _gr_tower_lazy_clear_data(&S->a3); _gr_tower_lazy_unref(L, S->a3.F);
            _gr_tower_lazy_clear_data(&S->a4); _gr_tower_lazy_unref(L, S->a4.F);
            _gr_tower_lazy_clear_data(&S->e2); _gr_tower_lazy_unref(L, S->e2.F);
        }
        flint_free(S);
    }
    flint_free(R);
}

/* the references of the data of a record (the collection of the records,
   lazy.c: all the records are counted with delta = -1, then with +1; the
   chain shared by several records is counted once) */
static void
_mod_link_rule_refs(void * data, slong delta)
{
    mod_link_rule_struct * R = data;
    mod_link_struct * S = R->S;
    _gr_tower_lazy_ref_adjust(R->scale.F, delta);
    if ((delta < 0) ? (S->visits++ == 0) : (--S->visits == 0))
    {
        _gr_tower_lazy_ref_adjust(S->tau1.F, delta);
        if (S->have == 1)
        {
            _gr_tower_lazy_ref_adjust(S->a2.F, delta);
            _gr_tower_lazy_ref_adjust(S->a3.F, delta);
            _gr_tower_lazy_ref_adjust(S->a4.F, delta);
            _gr_tower_lazy_ref_adjust(S->e2.F, delta);
        }
    }
}

/* the value of a linked generator through the chain */
static int
_mod_link_value(gr_ptr res, void * data, gr_ctx_t ctx)
{
    mod_link_rule_struct * R = data;
    mod_link_struct * S = R->S;
    gr_ptr u, v, L, K;
    int status = GR_SUCCESS;

    if (S->have == 0)
    {
        gr_tower_lazy_elem_struct * t[4] = {&S->a2, &S->a3, &S->a4, &S->e2};
        slong k;
        for (k = 0; k < 4; k++)
            _gr_tower_lazy_init(t[k], ctx);
        status = _mod_from_anchor(&S->a2, &S->a3, &S->a4, &S->e2, S->gam, &S->tau1, 0, ctx);
        S->have = (status == GR_SUCCESS) ? 1 : -1;
        if (status != GR_SUCCESS)
            for (k = 0; k < 4; k++)
                _gr_tower_lazy_clear(t[k], ctx);
    }
    if (S->have != 1)
        return GR_UNABLE;

    GR_TMP_INIT4(u, v, L, K, ctx);
    switch (R->which)
    {
        case LINK_R2:
            status = gr_div(res, &S->a2, &S->a3, ctx);
            break;
        case LINK_R4:
            status = gr_div(res, &S->a4, &S->a3, ctx);
            break;
        case LINK_TH3:
            status = gr_set(res, &S->a3, ctx);
            break;
        case LINK_K:
            /* K = pi theta_3^2 / 2 */
            status = gr_pi(u, ctx);
            status |= gr_sqr(res, &S->a3, ctx);
            status |= gr_mul(res, res, u, ctx);
            status |= gr_div_ui(res, res, 2, ctx);
            break;
        default:
            /* lambda = (theta_2 / theta_3)^4, E = K (E_2 / theta_3^4 + 2 - lambda) / 3 */
            status = gr_div(u, &S->a2, &S->a3, ctx);
            status |= gr_pow_ui(L, u, 4, ctx);
            if (R->which == LINK_LAMBDA)
                status |= gr_set(res, L, ctx);
            else
            {
                status |= gr_pi(u, ctx);
                status |= gr_sqr(K, &S->a3, ctx);
                status |= gr_mul(K, K, u, ctx);
                status |= gr_div_ui(K, K, 2, ctx);
                status |= gr_sqr(v, &S->a3, ctx);
                status |= gr_sqr(v, v, ctx);
                status |= gr_div(res, &S->e2, v, ctx);
                status |= gr_add_ui(res, res, 2, ctx);
                status |= gr_sub(res, res, L, ctx);
                status |= gr_mul(res, res, K, ctx);
                status |= gr_div_ui(res, res, 3, ctx);
            }
    }
    status |= gr_mul(res, res, &R->scale, ctx);
    GR_TMP_CLEAR4(u, v, L, K, ctx);
    return status;
}

/* the record of a linked generator: x (a new value at tau0) is a
   monomial whose newest generator g has exponent 1; the generator g =
   (g / x) * (the value of x through the chain) */
static int
_mod_link_rule(gr_srcptr x_in, mod_link_struct * S, int which, ulong trigger, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * x = _gr_tower_lazy_flat_view(x_in);
    gr_tower_flat_struct * F = x->F;
    const fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(&x->elem.flat.data);
    fmpz_mpoly_ctx_struct * mctx = x->elem.flat.mctx;
    mod_link_rule_struct * R;
    ulong * exp;
    slong nvars, v, d = -1;
    ulong def;
    gr_ptr g;
    int status = GR_SUCCESS;

    if (F == LAZY(ctx)->trivial || num->length != 1 || mctx != F->mctx)
        return GR_UNABLE;

    /* the newest generator of the monomial */
    nvars = mctx->minfo->nvars;
    exp = flint_malloc(sizeof(ulong) * nvars);
    fmpz_mpoly_get_term_exp_ui(exp, num, 0, mctx);
    for (v = 0; v < nvars; v++)
        if (exp[v] != 0 && F->cap - 1 - v < F->T->num_gens)
            d = FLINT_MAX(d, F->cap - 1 - v);
    if (d < 0 || exp[F->cap - 1 - d] != 1 ||
        fmpz_mpoly_degree_si(fmpz_mpoly_q_denref(&x->elem.flat.data), F->cap - 1 - d, mctx) > 0)
        status = GR_UNABLE;
    flint_free(exp);
    if (status != GR_SUCCESS)
        return status;

    _gr_tower_lazy_assign_def_ids(F, ctx);
    def = GR_TOWER_GEN(F->T, d)->def_id;
    if (def == 0)
        return GR_UNABLE;

    R = flint_malloc(sizeof(mod_link_rule_struct));
    R->S = S;
    R->which = which;
    _gr_tower_lazy_init(&R->scale, ctx);
    GR_TMP_INIT(g, ctx);
    _gr_tower_lazy_set_gen_d(g, F, d, ctx);
    status = gr_div(&R->scale, g, x_in, ctx);
    GR_TMP_CLEAR(g, ctx);
    S->refs++;
    if (status != GR_SUCCESS)
    {
        _mod_link_rule_clear(R, LAZY(ctx));
        return status;
    }
    _gr_tower_lazy_rebase_add_lazy(def, trigger, _mod_link_value, R, _mod_link_rule_clear, _mod_link_rule_refs, ctx);
    return GR_SUCCESS;
}

/* the link of the new generator L0 = lambda(tau0) to the anchor of the
   entry ent, whose root has the definition root: the new generators of
   the values at tau0 are created, with their records */
static int
_mod_link(gr_srcptr L0, slong ent, ulong root, gr_ctx_t ctx)
{
    mod_link_struct * S;
    gr_ptr x, u, v;
    fmpq_t q;
    int status;

    if (root == 0)
        return GR_UNABLE;

    S = flint_malloc(sizeof(mod_link_struct));
    S->refs = 1;    /* (held here until the records are made) */
    S->visits = 0;
    S->have = 0;
    _gr_tower_lazy_init(&S->tau1, ctx);
    GR_MUST_SUCCEED(gr_set(&S->tau1, &LAZY(ctx)->mod_pts[ent].tau1, ctx));
    S->gam[0] = S->gam[1] = S->gam[2] = S->gam[3] = 0;
    _fmpz_vec_set(S->gam, LAZY(ctx)->mod_pts[ent].gam, 4);

    GR_TMP_INIT3(x, u, v, ctx);
    fmpq_init(q);
    fmpq_set_si(q, 1, 4);

    status = _mod_link_rule(L0, S, LINK_LAMBDA, root, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_pow_fmpq(x, L0, q, ctx);
        if (status == GR_SUCCESS)
            status = _mod_link_rule(x, S, LINK_R2, root, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = gr_sub_ui(v, L0, 1, ctx);
        status |= gr_neg(v, v, ctx);
        status |= gr_pow_fmpq(x, v, q, ctx);
        if (status == GR_SUCCESS)
            status = _mod_link_rule(x, S, LINK_R4, root, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = gr_elliptic_k(x, L0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_link_rule(x, S, LINK_K, root, ctx);
    }
    if (status == GR_SUCCESS)
    {
        /* theta_3 = sqrt(2 K / pi), as in _mod_base */
        status = gr_elliptic_k(x, L0, ctx);
        status |= gr_pi(u, ctx);
        status |= gr_div(x, x, u, ctx);
        status |= gr_mul_ui(x, x, 2, ctx);
        status |= gr_sqrt(x, x, ctx);
        if (status == GR_SUCCESS)
            status = _mod_link_rule(x, S, LINK_TH3, root, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = gr_elliptic_e(x, L0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_link_rule(x, S, LINK_E, root, ctx);
    }

    {
        /* (the reference held here) */
        mod_link_rule_struct * R = flint_malloc(sizeof(mod_link_rule_struct));
        R->S = S;
        _gr_tower_lazy_init(&R->scale, ctx);
        _mod_link_rule_clear(R, LAZY(ctx));
    }

    GR_TMP_CLEAR3(x, u, v, ctx);
    fmpq_clear(q);
    return status;
}

/*
    The values at tau0 in F (any of the outputs may be NULL): theta_2,
    theta_3, theta_4, Lambda = lambda(tau0), E_2; with ratios, t2 and t4
    are set to theta_2 / theta_3 and theta_4 / theta_3 (t3 and e2 are not
    used), which need no elliptic integral. Lambda is algebraic at CM
    points; at a point commensurable with an anchor, everything comes
    from the anchor; otherwise Lambda is a generator (and tau0 becomes an
    anchor).
*/
static int
_mod_point_compute(gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_ptr lam, gr_ptr e2, gr_srcptr tau0, int ratios, slong ent, gr_ctx_t ctx)
{
    gr_ptr L, tau1, a2, a3, a4;
    fmpz * gam;
    int found, status, anchored = 0;

    GR_TMP_INIT5(L, tau1, a2, a3, a4, ctx);
    gam = _fmpz_vec_init(4);

    if (ratios)
    {
        t3 = NULL;
        e2 = NULL;
    }

    status = _mod_cm_lambda(L, &found, tau0, ctx);

    /* the anchor, chosen once per point (the cache entry ent): values
       computed later at tau0 (other outputs) come from the same anchor,
       whatever anchors appeared meanwhile */
    if (status == GR_SUCCESS && !found)
    {
        if (LAZY(ctx)->mod_pts[ent].anchor == 0)
        {
            ulong root = 0;
            int heavy = 0;
            int a = _mod_find_anchor(tau1, gam, &root, &heavy, tau0, ctx) &&
                !(fmpz_is_one(gam + 0) && fmpz_is_zero(gam + 1) && fmpz_is_zero(gam + 2) && fmpz_is_one(gam + 3));
            /* (the entry may have moved: by its index) */
            gr_tower_lazy_mod_point_struct * e = LAZY(ctx)->mod_pts + ent;
            if (a)
            {
                _gr_tower_lazy_init(&e->tau1, ctx);
                GR_MUST_SUCCEED(gr_set(&e->tau1, tau1, ctx));
                e->gam[0] = e->gam[1] = e->gam[2] = e->gam[3] = 0;   /* (fmpz init) */
                _fmpz_vec_set(e->gam, gam, 4);
                e->anchor = heavy ? 2 : 1;
                e->root = root;
            }
            else
                e->anchor = -1;
        }
        if (LAZY(ctx)->mod_pts[ent].anchor == 1)
        {
            anchored = 1;
            status = gr_set(tau1, &LAZY(ctx)->mod_pts[ent].tau1, ctx);
            _fmpz_vec_set(gam, LAZY(ctx)->mod_pts[ent].gam, 4);
        }
    }

    if (status == GR_SUCCESS && !found && anchored)
    {
        int need = (t2 != NULL || t3 != NULL || t4 != NULL || e2 != NULL);
        status = _mod_from_anchor(a2, a3, a4, e2, gam, tau1, ratios || !need, ctx);
        if (status == GR_SUCCESS && lam != NULL)
        {
            /* Lambda = theta_2^4 / theta_3^4 */
            status = gr_div(L, a2, a3, ctx);
            status |= gr_pow_ui(lam, L, 4, ctx);
        }
        if (status == GR_SUCCESS && t2 != NULL)
            status = ratios ? gr_div(t2, a2, a3, ctx) : gr_set(t2, a2, ctx);
        if (status == GR_SUCCESS && t3 != NULL)
            status = gr_set(t3, a3, ctx);
        if (status == GR_SUCCESS && t4 != NULL)
            status = ratios ? gr_div(t4, a4, a3, ctx) : gr_set(t4, a4, ctx);
    }
    else if (status == GR_SUCCESS)
    {
        if (!found)
        {
            status = _gr_tower_lazy_special_gen_locked(L, tau0, GR_TOWER_MODULAR_LAMBDA, 0, ctx);
            /* a new generator linked to the anchor (once): if the link
               fails, the values come from the anchor after all */
            if (status == GR_SUCCESS && LAZY(ctx)->mod_pts[ent].anchor == 2)
            {
                if (_mod_link(L, ent, LAZY(ctx)->mod_pts[ent].root, ctx) != GR_SUCCESS)
                {
                    LAZY(ctx)->mod_pts[ent].anchor = 1;
                    GR_TMP_CLEAR5(L, tau1, a2, a3, a4, ctx);
                    _fmpz_vec_clear(gam, 4);
                    return _mod_point_compute(t2, t3, t4, lam, e2, tau0, ratios, ent, ctx);
                }
                LAZY(ctx)->mod_pts[ent].anchor = 3;
            }
            if (status == GR_SUCCESS)
                LAZY(ctx)->mod_pts[ent].root = _gr_tower_lazy_gen_def_of(L, ctx);
        }
        if (status == GR_SUCCESS && ratios)
        {
            fmpq_t q;
            fmpq_init(q);
            fmpq_set_si(q, 1, 4);
            if (t2 != NULL)
                status = gr_pow_fmpq(t2, L, q, ctx);
            if (t4 != NULL)
            {
                status |= gr_sub_ui(t4, L, 1, ctx);
                status |= gr_neg(t4, t4, ctx);
                status |= gr_pow_fmpq(t4, t4, q, ctx);
            }
            fmpq_clear(q);
        }
        else if (status == GR_SUCCESS && (t2 != NULL || t3 != NULL || t4 != NULL || e2 != NULL))
            status = _mod_base(t2, t3, t4, e2, L, ctx);
        if (status == GR_SUCCESS && lam != NULL)
            status = gr_set(lam, L, ctx);
    }

    GR_TMP_CLEAR5(L, tau1, a2, a3, a4, ctx);
    _fmpz_vec_clear(gam, 4);
    return status;
}

/* the entry of the reduced point tau0 in the cache of the context
   (created if absent) */
static slong
_mod_point_entry(gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * Lz = LAZY(ctx);
    acb_t a, b;
    slong k, res = -1;

    acb_init(a);
    acb_init(b);
    if (gr_tower_lazy_get_acb(a, tau0, 128, ctx) == GR_SUCCESS)
    {
        for (k = 0; k < Lz->num_mod_pts && res < 0; k++)
        {
            if (gr_tower_lazy_get_acb(b, &Lz->mod_pts[k].tau0, 128, ctx) != GR_SUCCESS || !acb_overlaps(a, b))
                continue;
            if (gr_equal(tau0, &Lz->mod_pts[k].tau0, ctx) == T_TRUE)
                res = k;
        }
    }
    acb_clear(a);
    acb_clear(b);

    if (res < 0)
    {
        gr_tower_lazy_mod_point_struct * e;
        slong m;
        if (Lz->num_mod_pts == Lz->alloc_mod_pts)
        {
            slong alloc = FLINT_MAX(8, 2 * Lz->alloc_mod_pts);
            gr_tower_lazy_mod_point_struct * p = flint_malloc(sizeof(gr_tower_lazy_mod_point_struct) * alloc);
            for (k = 0; k < Lz->num_mod_pts; k++)
                p[k] = Lz->mod_pts[k];
            flint_free(Lz->mod_pts);
            Lz->mod_pts = p;
            Lz->alloc_mod_pts = alloc;
        }
        e = Lz->mod_pts + Lz->num_mod_pts;
        _gr_tower_lazy_init(&e->tau0, ctx);
        GR_MUST_SUCCEED(gr_set(&e->tau0, tau0, ctx));
        for (m = 0; m < 2; m++)
            for (k = 0; k < 5; k++)
                e->have[m][k] = 0;
        e->anchor = 0;
        e->root = 0;
        e->stamp = 0;
        res = Lz->num_mod_pts++;
    }
    return res;
}

/*
    A cached element, moved to a fresh copy of its prefix when its tower
    has grown since (towers grow in place: the cached values would
    otherwise drag every generator adjoined to that tower later into the
    computations which use them).
*/
static void
_gr_tower_lazy_cache_refresh(gr_tower_lazy_elem_struct * x, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * v = _gr_tower_lazy_flat_view(x);
    if (v->level > 0 && v->F->T->num_gens - v->level >= LAZY_FORK_SUFFIX)
    {
        gr_tower_lazy_elem_struct y;
        _gr_tower_lazy_init(&y, ctx);
        _gr_tower_lazy_prefix_copy(&y, v, ctx);
        _gr_tower_lazy_swap(x, &y, ctx);
        _gr_tower_lazy_clear(&y, ctx);
    }
}

/*
    The caches keep the values of the most recently used points only:
    cached elements keep their towers alive (the garbage collection
    cannot reclaim them), and with them generators which the searches of
    later operations visit (the arguments of the hypergeometric
    generators, the anchors).
*/
#define MOD_POINT_CACHE 8
#define MOD_POINT_VALUES_CACHE 1
#define THETA_POINT_CACHE 32

static void
_mod_point_evict(slong keep, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    slong i, n, best;
    int m, k;

    /* (the values, m = 0, involve the elliptic integrals at lambda: the
       generators K(lambda), E(lambda) kept alive in a long session make
       the later hypergeometric normal forms compare with them, and costly
       exact comparisons follow: the values of one point only) */
    for (m = 0; m < 2; m++)
    {
        slong limit = (m == 0) ? MOD_POINT_VALUES_CACHE : MOD_POINT_CACHE;
        for (;;)
        {
            n = 0;
            best = -1;
            for (i = 0; i < L->num_mod_pts; i++)
            {
                int any = 0;
                for (k = 0; k < 5; k++)
                    any |= L->mod_pts[i].have[m][k];
                if (!any)
                    continue;
                n++;
                if (i != keep && (best < 0 || L->mod_pts[i].stamp < L->mod_pts[best].stamp))
                    best = i;
            }
            if (n <= limit || best < 0)
                break;
            for (k = 0; k < 5; k++)
                if (L->mod_pts[best].have[m][k])
                {
                    _gr_tower_lazy_clear(&L->mod_pts[best].v[m][k], ctx);
                    L->mod_pts[best].have[m][k] = 0;
                }
        }
    }
}

static void
_mod_tp_evict(slong keep, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    slong i, n, best, q;

    for (;;)
    {
        n = 0;
        best = -1;
        for (i = 0; i < L->num_theta_pts; i++)
        {
            if (!L->theta_pts[i].have_vals)
                continue;
            n++;
            if (i != keep && !L->theta_pts[i].pinned &&
                (best < 0 || L->theta_pts[i].stamp < L->theta_pts[best].stamp))
                best = i;
        }
        if (n <= THETA_POINT_CACHE || best < 0)
            break;
        for (q = 0; q < 4; q++)
            _gr_tower_lazy_clear(&L->theta_pts[best].vals[q], ctx);
        L->theta_pts[best].have_vals = 0;
    }
}

/*
    theta_2, theta_3, theta_4 (or with ratios theta_2 / theta_3, theta_4 /
    theta_3), lambda and E_2 at the reduced point tau0 (any output may be
    NULL), cached in the context: repeated evaluations give the same
    elements (and the merges of their towers are found again).
*/
static int
_mod_point(gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_ptr lam, gr_ptr e2, gr_srcptr tau0, int ratios, gr_ctx_t ctx)
{
    gr_ptr out[5], tmp[5];
    slong k, m, ent;
    int status = GR_SUCCESS, missing = 0;

    m = (ratios != 0);
    out[0] = t2; out[1] = ratios ? NULL : t3; out[2] = t4; out[3] = lam; out[4] = ratios ? NULL : e2;


    ent = _mod_point_entry(tau0, ctx);

    /* the cached outputs first: the computation of the others may
       evaluate at other points, whose caching may evict this entry */
    for (k = 0; k < 5; k++)
    {
        tmp[k] = NULL;
        if (out[k] == NULL)
            continue;
        if (LAZY(ctx)->mod_pts[ent].have[m][k])
        {
            _gr_tower_lazy_cache_refresh(&LAZY(ctx)->mod_pts[ent].v[m][k], ctx);
            if (status == GR_SUCCESS)
                status = gr_set(out[k], &LAZY(ctx)->mod_pts[ent].v[m][k], ctx);
        }
        else
        {
            missing = 1;
            tmp[k] = out[k];
        }
    }

    if (missing && status == GR_SUCCESS)
    {
        status = _mod_point_compute(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tau0, ratios, ent, ctx);
        for (k = 0; k < 5 && status == GR_SUCCESS; k++)
        {
            /* (the entry may have moved: by its index) */
            gr_tower_lazy_mod_point_struct * e = LAZY(ctx)->mod_pts + ent;
            if (tmp[k] == NULL || e->have[m][k])
                continue;
            _gr_tower_lazy_init(&e->v[m][k], ctx);
            GR_MUST_SUCCEED(gr_set(&e->v[m][k], tmp[k], ctx));
            e->have[m][k] = 1;
        }
    }

    LAZY(ctx)->mod_pts[ent].stamp = ++LAZY(ctx)->cache_clock;
    if (missing)
        _mod_point_evict(ent, ctx);

    return status;
}

/* Lambda = lambda(tau0): algebraic at CM points, otherwise a generator
   (or through an anchor) */
static int
_mod_lambda0(gr_ptr res, gr_srcptr tau0, gr_ctx_t ctx)
{
    return _mod_point(NULL, NULL, NULL, res, NULL, tau0, 0, ctx);
}

/*
    The theta constants theta_2, theta_3, theta_4 at tau0 in F (any of
    the outputs may be NULL), and Lambda.
*/
static int
_mod_theta0(gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_ptr lam, gr_srcptr tau0, gr_ctx_t ctx)
{
    return _mod_point(t2, t3, t4, lam, NULL, tau0, 0, ctx);
}

/* eta(tau0) = (theta_2 theta_3 theta_4 / 2)^(1/3) */
static int
_mod_eta0(gr_ptr res, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr t2, t3, t4;
    fmpq_t e;
    int status;

    GR_TMP_INIT3(t2, t3, t4, ctx);
    fmpq_init(e);
    status = _mod_theta0(t2, t3, t4, NULL, tau0, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_mul(res, t2, t3, ctx);
        status |= gr_mul(res, res, t4, ctx);
        status |= gr_div_ui(res, res, 2, ctx);
        fmpq_set_si(e, 1, 3);
        status |= gr_pow_fmpq(res, res, e, ctx);
    }
    fmpq_clear(e);
    GR_TMP_CLEAR3(t2, t3, t4, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* transformation laws                                                   */
/* -------------------------------------------------------------------- */

/* res = exp(pi i r / n) (a root of unity) */
static int
_mod_root_of_unity(gr_ptr res, slong r, slong n, gr_ctx_t ctx)
{
    fmpq_t q;
    int status;
    fmpq_init(q);
    /* exp(pi i r / n) = exp(2 pi i (r / (2n))) */
    fmpq_set_si(q, r, 2 * n);
    status = gr_set_fmpq(res, q, ctx);
    if (status == GR_SUCCESS)
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status = gr_pi(t, ctx);
        status |= gr_mul(res, res, t, ctx);
        status |= gr_i(t, ctx);
        status |= gr_mul(res, res, t, ctx);
        status |= gr_mul_ui(res, res, 2, ctx);
        status |= gr_exp(res, res, ctx);
        GR_TMP_CLEAR(t, ctx);
    }
    fmpq_clear(q);
    return status;
}

/*
    res = exp(pi i alpha), with the root of unity exp(pi i Re(alpha))
    split off when Re(alpha) is rational: the remaining factor exp(-pi
    Im(alpha)) is real, so that the factors of the transformation laws at
    rational points are roots of unity times powers of exp(pi) rather
    than exponentials mixing both, whose multiplicative relations have
    large degrees.
*/
static int
_mod_exp_pi_i(gr_ptr res, gr_srcptr alpha, gr_ctx_t ctx)
{
    gr_ptr t, u;
    fmpq_t a;
    int status, split = 0;

    GR_TMP_INIT2(t, u, ctx);
    fmpq_init(a);

    status = gr_re(t, alpha, ctx);
    if (status == GR_SUCCESS && _gr_tower_lazy_rational_recognize(a, t, ctx) == 1 && !fmpq_is_zero(a) &&
        fmpz_fits_si(fmpq_numref(a)) && fmpz_fits_si(fmpq_denref(a)))
    {
        /* exp(pi i a) exp(-pi Im(alpha)) */
        status = _mod_root_of_unity(u, fmpz_get_si(fmpq_numref(a)), fmpz_get_si(fmpq_denref(a)), ctx);
        status |= gr_im(t, alpha, ctx);
        status |= gr_neg(t, t, ctx);
        status |= gr_pi(res, ctx);
        status |= gr_mul(t, t, res, ctx);
        status |= gr_exp(t, t, ctx);
        status |= gr_mul(res, t, u, ctx);
        split = 1;
    }

    if (!split)
    {
        status = gr_pi(u, ctx);
        status |= gr_mul(t, alpha, u, ctx);
        status |= gr_i(u, ctx);
        status |= gr_mul(t, t, u, ctx);
        status |= gr_exp(res, t, ctx);
    }

    GR_TMP_CLEAR2(t, u, ctx);
    fmpq_clear(a);
    return status;
}

/*
    The theta constants at tau (any output may be NULL): with tau0 = g tau,
    theta_{1+i}(0, tau) = exp(pi i R_i / 4) A theta_{1+S_i}(0, tau0),
    A = sqrt(i / (c tau + d)) (1 if C = 0).
*/
static int
_mod_theta_consts(gr_ptr t2, gr_ptr t3, gr_ptr t4, gr_srcptr tau, gr_ctx_t ctx)
{
    psl2z_t g;
    gr_ptr tau0, th[4], A, u;
    int R[4], S[4], C, i;
    int status;

    psl2z_init(g);
    GR_TMP_INIT3(tau0, A, u, ctx);
    GR_TMP_INIT4(th[0], th[1], th[2], th[3], ctx);

    status = _gr_tower_lazy_modular_reduce(tau0, g, tau, ctx);
    if (status == GR_SUCCESS)
        status = _mod_theta0(th[1], th[2], th[3], NULL, tau0, ctx);
    status |= gr_zero(th[0], ctx);

    if (status == GR_SUCCESS)
    {
        acb_modular_theta_transform(R, S, &C, g);

        if (C)
        {
            /* A = sqrt(i / (c tau + d)) */
            status = _mod_cd(u, g, tau, ctx);
            status |= gr_inv(A, u, ctx);
            status |= gr_i(u, ctx);
            status |= gr_mul(A, A, u, ctx);
            status |= gr_sqrt(A, A, ctx);
        }
        else
            status = gr_one(A, ctx);

        for (i = 1; i < 4 && status == GR_SUCCESS; i++)
        {
            gr_ptr out = (i == 1) ? t2 : (i == 2) ? t3 : t4;
            if (out == NULL)
                continue;
            status = _mod_root_of_unity(u, R[i], 4, ctx);
            status |= gr_mul(u, u, A, ctx);
            status |= gr_mul(out, u, th[S[i]], ctx);
        }
    }

    psl2z_clear(g);
    GR_TMP_CLEAR3(tau0, A, u, ctx);
    GR_TMP_CLEAR4(th[0], th[1], th[2], th[3], ctx);
    return status;
}

/* lambda(tau) = theta_2^4 / theta_3^4: a rational function of Lambda */
static int
_mod_lambda(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx)
{
    psl2z_t g;
    gr_ptr tau0, L, q[4], u;
    int R[4], S[4], C;
    int status;

    psl2z_init(g);
    GR_TMP_INIT3(tau0, L, u, ctx);
    GR_TMP_INIT4(q[0], q[1], q[2], q[3], ctx);

    status = _gr_tower_lazy_modular_reduce(tau0, g, tau, ctx);
    if (status == GR_SUCCESS)
        status = _mod_lambda0(L, tau0, ctx);

    if (status == GR_SUCCESS)
    {
        /* theta_j(tau0)^4 / theta_3(tau0)^4: 0, Lambda, 1, 1 - Lambda */
        status = gr_zero(q[0], ctx);
        status |= gr_set(q[1], L, ctx);
        status |= gr_one(q[2], ctx);
        status |= gr_sub_ui(q[3], L, 1, ctx);
        status |= gr_neg(q[3], q[3], ctx);

        acb_modular_theta_transform(R, S, &C, g);

        /* lambda(tau) = exp(pi i (R_1 - R_2)) q[S_1] / q[S_2] */
        status |= gr_div(res, q[S[1]], q[S[2]], ctx);
        if ((R[1] - R[2]) % 2 != 0)
            status |= gr_neg(res, res, ctx);
    }

    psl2z_clear(g);
    GR_TMP_CLEAR3(tau0, L, u, ctx);
    GR_TMP_CLEAR4(q[0], q[1], q[2], q[3], ctx);
    return status;
}

/* eta(tau) = eta(tau0) / (eps(g) sqrt(c tau + d)) */
static int
_mod_eta(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx)
{
    psl2z_t g;
    gr_ptr tau0, u;
    int status;

    psl2z_init(g);
    GR_TMP_INIT2(tau0, u, ctx);

    status = _gr_tower_lazy_modular_reduce(tau0, g, tau, ctx);
    if (status == GR_SUCCESS)
        status = _mod_eta0(res, tau0, ctx);
    if (status == GR_SUCCESS && !psl2z_is_one(g))
    {
        status = _mod_root_of_unity(u, -acb_modular_epsilon_arg(g), 12, ctx);
        status |= gr_mul(res, res, u, ctx);
        status |= _mod_cd(u, g, tau, ctx);
        status |= gr_rsqrt(u, u, ctx);
        status |= gr_mul(res, res, u, ctx);
    }

    psl2z_clear(g);
    GR_TMP_CLEAR2(tau0, u, ctx);
    return status;
}

/* j(tau) = 256 (1 - L + L^2)^3 / (L (1 - L))^2, L = lambda(tau0) */
static int
_mod_j(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx)
{
    psl2z_t g;
    gr_ptr tau0, L, u, v;
    int status;

    psl2z_init(g);
    GR_TMP_INIT4(tau0, L, u, v, ctx);

    status = _gr_tower_lazy_modular_reduce(tau0, g, tau, ctx);
    if (status == GR_SUCCESS)
        status = _mod_lambda0(L, tau0, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_sqr(u, L, ctx);
        status |= gr_sub(u, u, L, ctx);
        status |= gr_add_ui(u, u, 1, ctx);
        status |= gr_pow_ui(u, u, 3, ctx);
        status |= gr_mul_ui(u, u, 256, ctx);
        status |= gr_sub_ui(v, L, 1, ctx);
        status |= gr_mul(v, v, L, ctx);
        status |= gr_sqr(v, v, ctx);
        if (status == GR_SUCCESS)
            status = gr_div(res, u, v, ctx);
    }

    psl2z_clear(g);
    GR_TMP_CLEAR4(tau0, L, u, v, ctx);
    return status;
}

/*
    The Eisenstein series E_k(tau) (normalized: constant term 1), k even:
    E_2, E_4, E_6 at tau0 from the theta constants, E_k for k >= 8 by the
    recurrence of the coefficients of the Weierstrass function, and the
    transformation laws E_k(g tau) = (c tau + d)^k E_k(tau) (k >= 4),
    E_2(g tau) = (c tau + d)^2 E_2(tau) + 6 c (c tau + d) / (pi i).
*/
static int
_mod_eisenstein0(gr_ptr res, ulong k, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr t3, L, u, v;
    int status;

    GR_TMP_INIT2(t3, L, ctx);
    GR_TMP_INIT2(u, v, ctx);

    if (k == 2)
        status = _mod_point(NULL, NULL, NULL, NULL, res, tau0, 0, ctx);
    else
        status = _mod_theta0(NULL, t3, NULL, L, tau0, ctx);

    if (k == 2)
        ;
    else if (status == GR_SUCCESS && k == 4)
    {
        /* theta_3^8 (1 - Lambda + Lambda^2) */
        status = gr_sqr(u, L, ctx);
        status |= gr_sub(u, u, L, ctx);
        status |= gr_add_ui(u, u, 1, ctx);
        status |= gr_pow_ui(v, t3, 8, ctx);
        status |= gr_mul(res, u, v, ctx);
    }
    else if (status == GR_SUCCESS && k == 6)
    {
        /* theta_3^12 (1 + Lambda) (2 - Lambda) (1 - 2 Lambda) / 2 */
        status = gr_add_ui(u, L, 1, ctx);
        status |= gr_sub_ui(v, L, 2, ctx);
        status |= gr_neg(v, v, ctx);
        status |= gr_mul(u, u, v, ctx);
        status |= gr_mul_ui(v, L, 2, ctx);
        status |= gr_sub_ui(v, v, 1, ctx);
        status |= gr_neg(v, v, ctx);
        status |= gr_mul(u, u, v, ctx);
        status |= gr_div_ui(u, u, 2, ctx);
        status |= gr_pow_ui(v, t3, 12, ctx);
        status |= gr_mul(res, u, v, ctx);
    }
    else if (status == GR_SUCCESS)
    {
        /*
            c_n = (2n + 1) zeta(2n + 2) E_{2n+2} / (2 pi^(2n+2))-free form:
            with G_{2k} = 2 zeta(2k) E_{2k} and the Weierstrass coefficients
            c_n = (2n + 1) G_{2n+2}, c_n = 3 / ((2n + 3)(n - 2)) sum_{m=1}^{n-2} c_m c_{n-1-m}
            for n >= 3 (c_1 = 3 G_4, c_2 = 5 G_6).
        */
        slong n, m, N = k / 2 - 1;
        gr_ptr c;
        fmpq_t z2;

        c = gr_heap_init_vec(N + 1, ctx);
        fmpq_init(z2);

        for (n = 1; n <= N && status == GR_SUCCESS; n++)
        {
            gr_ptr cn = GR_ENTRY(c, n, ctx->sizeof_elem);
            if (n <= 2)
            {
                /* G_{2n+2} = 2 zeta(2n+2) E_{2n+2} */
                status = _mod_eisenstein0(cn, 2 * n + 2, tau0, ctx);
                _gr_tower_zeta_even_over_pi(z2, 2 * n + 2);
                status |= gr_mul_fmpq(cn, cn, z2, ctx);
                status |= gr_mul_ui(cn, cn, 2 * (2 * n + 1), ctx);
                status |= gr_pi(u, ctx);
                status |= gr_pow_ui(u, u, 2 * n + 2, ctx);
                status |= gr_mul(cn, cn, u, ctx);
            }
            else
            {
                status |= gr_zero(cn, ctx);
                for (m = 1; m <= n - 2; m++)
                {
                    status |= gr_mul(u, GR_ENTRY(c, m, ctx->sizeof_elem), GR_ENTRY(c, n - 1 - m, ctx->sizeof_elem), ctx);
                    status |= gr_add(cn, cn, u, ctx);
                }
                status |= gr_mul_ui(cn, cn, 3, ctx);
                status |= gr_div_ui(cn, cn, (2 * n + 3) * (n - 2), ctx);
            }
        }

        /* E_k = c_N / ((2N + 1) 2 zeta(k)) */
        if (status == GR_SUCCESS)
        {
            status = gr_div_ui(res, GR_ENTRY(c, N, ctx->sizeof_elem), 2 * (2 * N + 1), ctx);
            _gr_tower_zeta_even_over_pi(z2, k);
            fmpq_inv(z2, z2);
            status |= gr_mul_fmpq(res, res, z2, ctx);
            status |= gr_pi(u, ctx);
            status |= gr_pow_ui(u, u, k, ctx);
            status |= gr_div(res, res, u, ctx);
        }

        gr_heap_clear_vec(c, N + 1, ctx);
        fmpq_clear(z2);
    }

    GR_TMP_CLEAR2(t3, L, ctx);
    GR_TMP_CLEAR2(u, v, ctx);
    return status;
}

static int
_mod_eisenstein(gr_ptr res, ulong k, gr_srcptr tau, gr_ctx_t ctx)
{
    psl2z_t g;
    gr_ptr tau0, u, v;
    int status;

    /* (as for acb: even k >= 2) */
    if (k == 0 || k % 2 == 1)
        return GR_DOMAIN;
    if (k > 1000)
        return GR_UNABLE;

    psl2z_init(g);
    GR_TMP_INIT3(tau0, u, v, ctx);

    status = _gr_tower_lazy_modular_reduce(tau0, g, tau, ctx);
    if (status == GR_SUCCESS)
        status = _mod_eisenstein0(res, k, tau0, ctx);

    if (status == GR_SUCCESS && !psl2z_is_one(g))
    {
        /* E_k(tau) = E_k(tau0) / (c tau + d)^k, minus the correction for k = 2 */
        status = _mod_cd(u, g, tau, ctx);
        if (k == 2 && !fmpz_is_zero(&g->c))
        {
            /* E_2(tau) = (E_2(tau0) - 6 c (c tau + d) / (pi i)) / (c tau + d)^2 */
            status |= gr_mul_fmpz(v, u, &g->c, ctx);
            status |= gr_mul_ui(v, v, 6, ctx);
            {
                gr_ptr p;
                GR_TMP_INIT(p, ctx);
                status |= gr_pi(p, ctx);
                status |= gr_div(v, v, p, ctx);
                status |= gr_i(p, ctx);
                status |= gr_div(v, v, p, ctx);
                GR_TMP_CLEAR(p, ctx);
            }
            status |= gr_sub(res, res, v, ctx);
        }
        status |= gr_pow_ui(u, u, k, ctx);
        if (status == GR_SUCCESS)
            status = gr_div(res, res, u, ctx);
    }

    psl2z_clear(g);
    GR_TMP_CLEAR3(tau0, u, v, ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* interface                                                             */
/* (the inputs are copied: res may alias them, and the evaluation reads
   them after writing res) */
#define MOD_WRAP1(name, call) \
int name(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t tau_in, gr_ctx_t ctx) \
{ \
    int status, real, alg; \
    gr_ptr tau; \
    _gr_tower_lazy_lock(ctx); \
    real = REAL(ctx); \
    alg = ALG(ctx); \
    GR_TMP_INIT(tau, ctx); \
    status = gr_set(tau, tau_in, ctx); \
    if (status == GR_SUCCESS) \
        status = _gr_tower_lazy_view_finish(call, res, real, alg, ctx); \
    GR_TMP_CLEAR(tau, ctx); \
    _gr_tower_lazy_unlock(ctx); \
    return status; \
}

MOD_WRAP1(gr_tower_lazy_modular_lambda, _mod_lambda(res, tau, ctx))
MOD_WRAP1(gr_tower_lazy_modular_j, _mod_j(res, tau, ctx))
MOD_WRAP1(gr_tower_lazy_dedekind_eta, _mod_eta(res, tau, ctx))

static int
_mod_delta(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx)
{
    int status = _mod_eta(res, tau, ctx);
    if (status == GR_SUCCESS)
        status = gr_pow_ui(res, res, 24, ctx);
    return status;
}

MOD_WRAP1(gr_tower_lazy_modular_delta, _mod_delta(res, tau, ctx))

int
gr_tower_lazy_eisenstein_e(gr_tower_lazy_elem_t res, ulong k, const gr_tower_lazy_elem_t tau_in, gr_ctx_t ctx)
{
    int status, real, alg;
    gr_ptr tau;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    GR_TMP_INIT(tau, ctx);
    status = gr_set(tau, tau_in, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_view_finish(_mod_eisenstein(res, k, tau, ctx), res, real, alg, ctx);
    GR_TMP_CLEAR(tau, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* G_k = 2 zeta(k) E_k */
int
gr_tower_lazy_eisenstein_g(gr_tower_lazy_elem_t res, ulong k, const gr_tower_lazy_elem_t tau_in, gr_ctx_t ctx)
{
    int status, real, alg;
    gr_ptr tau;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    GR_TMP_INIT(tau, ctx);
    status = gr_set(tau, tau_in, ctx);
    if (status == GR_SUCCESS)
        status = _mod_eisenstein(res, k, tau, ctx);
    if (status == GR_SUCCESS)
    {
        fmpq_t z;
        gr_ptr p;
        fmpq_init(z);
        GR_TMP_INIT(p, ctx);
        _gr_tower_zeta_even_over_pi(z, k);
        fmpq_mul_2exp(z, z, 1);
        status = gr_mul_fmpq(res, res, z, ctx);
        status |= gr_pi(p, ctx);
        status |= gr_pow_ui(p, p, k, ctx);
        status |= gr_mul(res, res, p, ctx);
        GR_TMP_CLEAR(p, ctx);
        fmpq_clear(z);
    }
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    GR_TMP_CLEAR(tau, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* theta functions of z                                                  */
/* -------------------------------------------------------------------- */

/* n = floor(x) for a real element x, exactly */
static int
_mod_floor(slong * n, gr_srcptr x, gr_ctx_t ctx)
{
    acb_t z;
    arb_t f;
    fmpz_t a;
    slong prec;
    int status = GR_UNABLE;

    acb_init(z);
    arb_init(f);
    fmpz_init(a);

    for (prec = CHECK_PREC; prec <= 4096; prec *= 4)
    {
        if (gr_tower_lazy_get_acb(z, x, prec, ctx) != GR_SUCCESS || !arb_is_finite(acb_realref(z)))
            break;
        if (!arb_contains_int(acb_realref(z)))
        {
            arb_floor(f, acb_realref(z), prec);
            if (arb_is_exact(f) && arf_is_int(arb_midref(f)))
            {
                arf_get_fmpz(a, arb_midref(f), ARF_RND_NEAR);
                if (fmpz_fits_si(a))
                {
                    *n = fmpz_get_si(a);
                    status = GR_SUCCESS;
                }
            }
            break;
        }
        /* an integer within the enclosure: is x that integer? */
        arf_get_fmpz(a, arb_midref(acb_realref(z)), ARF_RND_NEAR);
        if (fmpz_fits_si(a) && mag_cmp_2exp_si(arb_radref(acb_realref(z)), -2) < 0)
        {
            gr_ptr t;
            truth_t eq;
            GR_TMP_INIT(t, ctx);
            eq = (gr_set_fmpz(t, a, ctx) == GR_SUCCESS) ? gr_equal(t, x, ctx) : T_UNKNOWN;
            GR_TMP_CLEAR(t, ctx);
            if (eq == T_TRUE)
            {
                *n = fmpz_get_si(a);
                status = GR_SUCCESS;
                break;
            }
            if (eq == T_UNKNOWN)
                break;
            /* x != a: more precision separates them */
        }
    }

    acb_clear(z);
    arb_clear(f);
    fmpz_clear(a);
    return status;
}

/*
    The values theta_1, ..., theta_4 at (z0, tau0) for z0 != 0 reduced:
    theta_1(z0, tau0) and theta_4(z0, tau0) are generators
    (GR_TOWER_JACOBI_THETA), and theta_2(z0, tau0), theta_3(z0, tau0) are
    algebraic generators over them, by Jacobi's relations
    theta_2(z)^2 theta_4^2 = theta_4(z)^2 theta_2^2 - theta_1(z)^2 theta_3^2,
    theta_3(z)^2 theta_4^2 = theta_4(z)^2 theta_3^2 - theta_1(z)^2 theta_2^2
    (theta_k = theta_k(0, tau0)), divided by theta_3^2:
    theta_2(z)^2 sqrt(1 - L) = theta_4(z)^2 sqrt(L) - theta_1(z)^2,
    theta_3(z)^2 sqrt(1 - L) = theta_4(z)^2 - theta_1(z)^2 sqrt(L), with
    L = lambda(tau0). (As generators with the definition theta_j(z0, tau0)
    rather than square roots of these expressions, whose form depends on
    the representation of their ingredients: equal arguments give the
    same generators.)
*/
static int
_mod_theta_z0_gens(gr_ptr * th, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr lam, s, r, args, u, v, w;
    slong k;
    int status;

    GR_TMP_INIT3(lam, s, r, ctx);
    GR_TMP_INIT3(u, v, w, ctx);
    args = gr_heap_init_vec(2, ctx);

    /* sqrt(L) = (theta_2 / theta_3)^2 and sqrt(1 - L) = (theta_4 /
       theta_3)^2 with the theta constants of _mod_theta0, so that the
       relations involve the same generators */
    status = _mod_point(s, NULL, r, NULL, NULL, tau0, 1, ctx);
    status |= gr_sqr(s, s, ctx);
    status |= gr_sqr(r, r, ctx);
    status |= gr_set(args, z0, ctx);
    status |= gr_set(GR_ENTRY(args, 1, ctx->sizeof_elem), tau0, ctx);

    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_special_gen_multi_locked(th[0], args, 2, GR_TOWER_JACOBI_THETA, 1, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_special_gen_multi_locked(th[3], args, 2, GR_TOWER_JACOBI_THETA, 4, ctx);

    for (k = 2; k <= 3 && status == GR_SUCCESS; k++)
    {
        status = gr_sqr(u, th[3], ctx);
        status |= gr_sqr(v, th[0], ctx);
        if (k == 2)
            status |= gr_mul(u, u, s, ctx);
        else
            status |= gr_mul(v, v, s, ctx);
        status |= gr_sub(w, u, v, ctx);
        status |= gr_div(w, w, r, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_special_alg_gen_multi_locked(th[k - 1], args, 2, GR_TOWER_JACOBI_THETA, k, w, 2, ctx);
    }

    GR_TMP_CLEAR3(lam, s, r, ctx);
    GR_TMP_CLEAR3(u, v, w, ctx);
    gr_heap_clear_vec(args, 2, ctx);
    return status;
}

/* theta[a, b](z, tau) = sum exp(pi i (n + a/2)^2 tau + 2 pi i (n + a/2)(z + b/2)):
   theta_1 = -theta[1, 1], theta_2 = theta[1, 0], theta_3 = theta[0, 0],
   theta_4 = theta[0, 1] */
static void
_mod_char_of(int * a, int * b, int * sign, int k)
{
    *a = (k == 1 || k == 2);
    *b = (k == 1 || k == 4);
    *sign = (k == 1) ? -1 : 1;
}

static int
_mod_k_of_char(int a, int b)
{
    return (a && b) ? 1 : a ? 2 : b ? 4 : 3;
}

/* the sign of x - c/2 (quarter = 0) or of x - c/2 + 1/4 (quarter = 1),
   x real (exact) */
static int
_mod_half_sign(int * sgn, gr_srcptr x, slong c, int quarter, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_INIT(t, ctx);
    /* 4x - 2c + quarter */
    status = gr_mul_ui(t, x, 4, ctx);
    status |= gr_sub_si(t, t, 2 * c, ctx);
    status |= gr_add_si(t, t, quarter, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_fast(sgn, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/*
    theta_k(z, tau0) for k = 1, ..., 4 (outputs th[k-1]; NULL to skip),
    tau0 in F: z = z1 + (p + q tau0)/2 with z1 = x + y tau0, x, y in
    [-1/4, 1/4) (exact floors), then theta[a, b](z1 + (p + q tau0)/2) =
    exp(-pi i q^2 tau0 / 4 - pi i q (z1 + (b + p)/2)) theta[a + q, b + p](z1),
    with theta[a + 2, b] = theta[a, b], theta[a, b + 2] = (-1)^a theta[a, b],
    and the parity theta[a, b](-z) = (-1)^(ab) theta[a, b](z); the
    reduced point z0 = +-z1 is canonical (see below: p, q adjusted on the
    boundary of the box).
*/
static int
_mod_z_reduce(gr_ptr z0, gr_ptr z1, slong * pp, slong * qq, int * negp, gr_srcptr z, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr x, y, t;
    slong p = 0, q = 0;
    int neg = 0, status = GR_SUCCESS;

    GR_TMP_INIT3(x, y, t, ctx);

    /* y = Im(z) / Im(tau0), x = Re(z) - y Re(tau0) */
    status |= gr_im(y, z, ctx);
    status |= gr_im(t, tau0, ctx);
    status |= gr_div(y, y, t, ctx);
    status |= gr_re(t, tau0, ctx);
    status |= gr_mul(t, t, y, ctx);
    status |= gr_re(x, z, ctx);
    status |= gr_sub(x, x, t, ctx);

    /* p = floor(2x + 1/2), q = floor(2y + 1/2) */
    if (status == GR_SUCCESS)
    {
        status = gr_mul_ui(t, x, 4, ctx);
        status |= gr_add_ui(t, t, 1, ctx);
        status |= gr_div_ui(t, t, 2, ctx);
        if (status == GR_SUCCESS)
            status = _mod_floor(&p, t, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = gr_mul_ui(t, y, 4, ctx);
        status |= gr_add_ui(t, t, 1, ctx);
        status |= gr_div_ui(t, t, 2, ctx);
        if (status == GR_SUCCESS)
            status = _mod_floor(&q, t, ctx);
    }

    /* the representative: negation and the shifts by 1/2, tau0/2 to the
       canonical set y1 in (0, 1/4), x1 in [-1/4, 1/4), or y1 in {0, 1/4},
       x1 in [0, 1/4] (a fundamental domain for z -> +-z + (Z + tau0 Z)/2:
       negation maps the box [-1/4, 1/4)^2 to (-1/4, 1/4]^2, whose
       boundary points must be brought back consistently, or the reduced
       points of z and of -z, of P and of the reduced point of P would
       differ) */
    if (status == GR_SUCCESS)
    {
        int sy, sy4, sx, sx4;
        /* the signs of y1 = y - q/2, y1 + 1/4, x1 = x - p/2, x1 + 1/4 */
        status = _mod_half_sign(&sy, y, q, 0, ctx);
        status |= _mod_half_sign(&sy4, y, q, 1, ctx);
        status |= _mod_half_sign(&sx, x, p, 0, ctx);
        status |= _mod_half_sign(&sx4, x, p, 1, ctx);
        if (status == GR_SUCCESS)
        {
            if (sy > 0)
                neg = 0;
            else if (sy == 0)
                neg = (sx < 0);
            else if (sy4 == 0)
            {
                /* y1 = -1/4: (x1, -1/4) ~ (-x1, 1/4) ~ (x1, 1/4) */
                if (sx <= 0)
                    neg = 1;
                else
                {
                    neg = 0;
                    q -= 1;
                }
            }
            else
            {
                neg = 1;
                /* x1 = -1/4: -z1 has x = 1/4, shifted by -1/2 */
                if (sx4 == 0)
                    p -= 1;
            }
        }
    }

    /* z1 = z - (p + q tau0) / 2 */
    if (status == GR_SUCCESS)
    {
        status = gr_mul_si(t, tau0, q, ctx);
        status |= gr_add_si(t, t, p, ctx);
        status |= gr_div_ui(t, t, 2, ctx);
        status |= gr_sub(z1, z, t, ctx);
    }

    /* z0 = +-z1 */
    if (status == GR_SUCCESS)
        status = neg ? gr_neg(z0, z1, ctx) : gr_set(z0, z1, ctx);

    *pp = p;
    *qq = q;
    *negp = neg;
    GR_TMP_CLEAR3(x, y, t, ctx);
    return status;
}

/*
    The map of the theta functions at P to those at its reduced point z0
    (see _mod_z_reduce): theta_k(P) = sg[k-1] f[k-1] theta_{kk[k-1]}(z0),
    with kk a permutation of 1, ..., 4.
*/
static int
_mod_reduce_map(gr_ptr z0, gr_ptr * f, int * kk, int * sg, gr_srcptr P, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr t, u, z1;
    slong p, q, k;
    int neg = 0, status;

    GR_TMP_INIT3(t, u, z1, ctx);

    status = _mod_z_reduce(z0, z1, &p, &q, &neg, P, tau0, ctx);

    /* theta_k(P) = sign_k theta[a, b](z1 + (p + q tau0)/2) = sign_k f theta[a', b'](z1) */
    for (k = 1; k <= 4 && status == GR_SUCCESS; k++)
    {
        int a, b, sign, a2, b2;
        slong bp;

        _mod_char_of(&a, &b, &sign, k);

        /* f = exp(-pi i q^2 tau0 / 4 - pi i q (z1 + (b + p)/2)) */
        status = gr_set(t, z1, ctx);
        status |= gr_set_si(u, b + p, ctx);
        status |= gr_div_ui(u, u, 2, ctx);
        status |= gr_add(t, t, u, ctx);
        status |= gr_mul_si(t, t, q, ctx);
        status |= gr_mul_si(u, tau0, q, ctx);
        status |= gr_mul_si(u, u, q, ctx);
        status |= gr_div_ui(u, u, 4, ctx);
        status |= gr_add(t, t, u, ctx);
        status |= gr_neg(t, t, ctx);
        status |= _mod_exp_pi_i(f[k - 1], t, ctx);

        /* theta[a + q, b + p] reduced */
        a2 = (int) (((a + q) % 2 + 2) % 2);
        bp = b + p;
        b2 = (int) ((bp % 2 + 2) % 2);
        /* theta[a2, b2 + 2m] = (-1)^(a2 m) theta[a2, b2] */
        if (a2 && (((bp - b2) / 2) % 2 != 0))
            sign = -sign;

        /* theta[a2, b2](z1) = (-1)^(a2 b2 neg) theta[a2, b2](z0) */
        if (neg && a2 && b2)
            sign = -sign;

        kk[k - 1] = _mod_k_of_char(a2, b2);
        /* theta[a2, b2] = (kk == 1 ? -theta_1 : theta_kk) */
        if (kk[k - 1] == 1)
            sign = -sign;
        sg[k - 1] = sign;
    }

    GR_TMP_CLEAR3(t, u, z1, ctx);
    return status;
}

static int _mod_theta_z0(gr_ptr * th, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx);

static int
_mod_theta_z_reduced(gr_ptr * out, gr_srcptr z, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr z0, f[4], base[4];
    int kk[4], sg[4], k, status;
    truth_t zero;

    GR_TMP_INIT(z0, ctx);
    GR_TMP_INIT4(f[0], f[1], f[2], f[3], ctx);
    GR_TMP_INIT4(base[0], base[1], base[2], base[3], ctx);

    status = _mod_reduce_map(z0, f, kk, sg, z, tau0, ctx);

    /* the values at z0 */
    if (status == GR_SUCCESS)
    {
        zero = gr_is_zero(z0, ctx);
        if (zero == T_UNKNOWN)
            status = GR_UNABLE;
        else if (zero == T_TRUE)
        {
            status = _mod_theta0(base[1], base[2], base[3], NULL, tau0, ctx);
            status |= gr_zero(base[0], ctx);
        }
        else
            status = _mod_theta_z0(base, z0, tau0, ctx);
    }

    for (k = 0; k < 4 && status == GR_SUCCESS; k++)
    {
        if (out[k] == NULL)
            continue;
        status = gr_mul(out[k], f[k], base[kk[k] - 1], ctx);
        if (sg[k] < 0)
            status |= gr_neg(out[k], out[k], ctx);
    }

    GR_TMP_CLEAR(z0, ctx);
    GR_TMP_CLEAR4(f[0], f[1], f[2], f[3], ctx);
    GR_TMP_CLEAR4(base[0], base[1], base[2], base[3], ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* theta functions of z: multiplication, addition, torsion points        */
/* -------------------------------------------------------------------- */

/*
    With rho_k(z) = theta_k(z) / theta_3, c_k = theta_k(0) / theta_3
    (c_1 = 0, c_3 = 1) and r_k = theta_k(z) / theta_1(z), the addition
    formulas (y = n z, w = z)

        theta_1(y+w) theta_1(y-w) theta_4^2 = theta_1(y)^2 theta_4(w)^2 - theta_4(y)^2 theta_1(w)^2,
        theta_2(y+w) theta_2(y-w) theta_4^2 = theta_2(y)^2 theta_4(w)^2 - theta_3(y)^2 theta_1(w)^2,
        theta_3(y+w) theta_3(y-w) theta_4^2 = theta_3(y)^2 theta_4(w)^2 - theta_2(y)^2 theta_1(w)^2,
        theta_4(y+w) theta_4(y-w) theta_4^2 = theta_4(y)^2 theta_4(w)^2 - theta_1(y)^2 theta_1(w)^2

    and the duplication theta_1(2z) theta_2 theta_3 theta_4 = 2 theta_1(z)
    theta_2(z) theta_3(z) theta_4(z) give rho_k(n z) = rho_1(z)^(n^2)
    S_k(n) with S_k(0) = c_k, S_k(1) = r_k and

        S_1(n+1) S_1(n-1) c_4^2 = S_1(n)^2 r_4^2 - S_4(n)^2  (n >= 2; S_1(2) = 2 r_2 r_3 r_4 / (c_2 c_4)),
        S_4(n+1) S_4(n-1) c_4^2 = S_4(n)^2 r_4^2 - S_1(n)^2,
        S_2(n+1) S_2(n-1) c_4^2 = S_2(n)^2 r_4^2 - S_3(n)^2,
        S_3(n+1) S_3(n-1) c_4^2 = S_3(n)^2 r_4^2 - S_2(n)^2.

    The pair (1, 4) divides by theta_1, theta_4 at (n-1) z, the pair
    (2, 3) by theta_2, theta_3 there. Sets S[k-1] for the k of the pair
    (pair 0: k = 1, 4; pair 1: k = 2, 3).
*/
static int
_mod_mult_ratios(gr_ptr * S, int pair, gr_srcptr * r, gr_srcptr c2, gr_srcptr c4, slong n, gr_ctx_t ctx)
{
    gr_ptr A0, A1, A2, B0, B1, B2, r42, c42, t;
    slong m;
    int status = GR_SUCCESS;

    GR_TMP_INIT3(A0, A1, A2, ctx);
    GR_TMP_INIT3(B0, B1, B2, ctx);
    GR_TMP_INIT3(r42, c42, t, ctx);

    status |= gr_sqr(r42, r[3], ctx);
    status |= gr_sqr(c42, c4, ctx);

    if (pair == 0)
    {
        status |= gr_zero(A0, ctx);
        status |= gr_one(A1, ctx);
        status |= gr_set(B0, c4, ctx);
        status |= gr_set(B1, r[3], ctx);
    }
    else
    {
        status |= gr_set(A0, c2, ctx);
        status |= gr_set(A1, r[1], ctx);
        status |= gr_one(B0, ctx);
        status |= gr_set(B1, r[2], ctx);
    }

    for (m = 1; m < n && status == GR_SUCCESS; m++)
    {
        /* A2 = A(m+1), B2 = B(m+1) */
        if (pair == 0 && m == 1)
        {
            status = gr_mul(A2, r[1], r[2], ctx);
            status |= gr_mul(A2, A2, r[3], ctx);
            status |= gr_mul_ui(A2, A2, 2, ctx);
            status |= gr_mul(t, c2, c4, ctx);
            status |= gr_div(A2, A2, t, ctx);
        }
        else
        {
            status = gr_sqr(A2, A1, ctx);
            status |= gr_mul(A2, A2, r42, ctx);
            status |= gr_sqr(t, B1, ctx);
            status |= gr_sub(A2, A2, t, ctx);
            status |= gr_mul(t, A0, c42, ctx);
            status |= gr_div(A2, A2, t, ctx);
        }
        status |= gr_sqr(B2, B1, ctx);
        status |= gr_mul(B2, B2, r42, ctx);
        status |= gr_sqr(t, A1, ctx);
        status |= gr_sub(B2, B2, t, ctx);
        status |= gr_mul(t, B0, c42, ctx);
        status |= gr_div(B2, B2, t, ctx);

        gr_swap(A0, A1, ctx); gr_swap(A1, A2, ctx);
        gr_swap(B0, B1, ctx); gr_swap(B1, B2, ctx);
    }

    if (status == GR_SUCCESS)
    {
        if (pair == 0)
        {
            status = gr_set(S[0], A1, ctx);
            status |= gr_set(S[3], B1, ctx);
        }
        else
        {
            status = gr_set(S[1], A1, ctx);
            status |= gr_set(S[2], B1, ctx);
        }
    }

    GR_TMP_CLEAR3(A0, A1, A2, ctx);
    GR_TMP_CLEAR3(B0, B1, B2, ctx);
    GR_TMP_CLEAR3(r42, c42, t, ctx);
    return status;
}

/*
    The division polynomial psi_n of the lattice, psi_n(z) = sigma(n z) /
    sigma(z)^(n^2), as f_n(x) with psi_n = y^(n mod 2 == 0) f_n(x), x = wp(z),
    y = wp'(z), y^2 = 4 x^3 - g2 x - g3:

        psi_2 = -y, psi_3 = 3 x^4 - 3/2 g2 x^2 - 3 g3 x - g2^2/16,
        psi_4 = -y (2 x^6 - 5/2 g2 x^4 - 10 g3 x^3 - 5/8 g2^2 x^2 - 1/2 g2 g3 x - g3^2 + g2^3/32),
        psi_(2m+1) = psi_(m+2) psi_m^3 - psi_(m-1) psi_(m+1)^3,
        psi_(2m) = psi_m (psi_(m+2) psi_(m-1)^2 - psi_(m-2) psi_(m+1)^2) / psi_2.
*/
#define TORSION_LEVEL 6

static int
_mod_division_polys(gr_poly_struct * res, slong N, gr_srcptr g2, gr_srcptr g3, gr_ctx_t ctx)
{
    gr_poly_struct f[TORSION_LEVEL + 3];
    gr_poly_t F, t, u;
    gr_ptr c;
    slong n, i;
    int status = GR_SUCCESS;

    if (N < 1 || N > TORSION_LEVEL + 1)
        return GR_UNABLE;

    for (i = 0; i < TORSION_LEVEL + 3; i++)
        gr_poly_init(f + i, ctx);
    gr_poly_init(F, ctx);
    gr_poly_init(t, ctx);
    gr_poly_init(u, ctx);
    GR_TMP_INIT(c, ctx);

    /* F = 4 x^3 - g2 x - g3 */
    status |= gr_poly_set_coeff_si(F, 3, 4, ctx);
    status |= gr_neg(c, g2, ctx);
    status |= gr_poly_set_coeff_scalar(F, 1, c, ctx);
    status |= gr_neg(c, g3, ctx);
    status |= gr_poly_set_coeff_scalar(F, 0, c, ctx);

    status |= gr_poly_zero(f + 0, ctx);
    status |= gr_poly_one(f + 1, ctx);
    status |= gr_poly_set_si(f + 2, -1, ctx);

    /* f_3 */
    status |= gr_poly_set_coeff_si(f + 3, 4, 3, ctx);
    status |= gr_mul_si(c, g2, -3, ctx);
    status |= gr_div_ui(c, c, 2, ctx);
    status |= gr_poly_set_coeff_scalar(f + 3, 2, c, ctx);
    status |= gr_mul_si(c, g3, -3, ctx);
    status |= gr_poly_set_coeff_scalar(f + 3, 1, c, ctx);
    status |= gr_sqr(c, g2, ctx);
    status |= gr_div_si(c, c, -16, ctx);
    status |= gr_poly_set_coeff_scalar(f + 3, 0, c, ctx);

    /* f_4 */
    {
        gr_ptr d;
        GR_TMP_INIT(d, ctx);
        status |= gr_poly_set_coeff_si(f + 4, 6, -2, ctx);
        status |= gr_mul_ui(c, g2, 5, ctx);
        status |= gr_div_ui(c, c, 2, ctx);
        status |= gr_poly_set_coeff_scalar(f + 4, 4, c, ctx);
        status |= gr_mul_ui(c, g3, 10, ctx);
        status |= gr_poly_set_coeff_scalar(f + 4, 3, c, ctx);
        status |= gr_sqr(c, g2, ctx);
        status |= gr_mul_ui(c, c, 5, ctx);
        status |= gr_div_ui(c, c, 8, ctx);
        status |= gr_poly_set_coeff_scalar(f + 4, 2, c, ctx);
        status |= gr_mul(c, g2, g3, ctx);
        status |= gr_div_ui(c, c, 2, ctx);
        status |= gr_poly_set_coeff_scalar(f + 4, 1, c, ctx);
        status |= gr_sqr(c, g3, ctx);
        status |= gr_pow_ui(d, g2, 3, ctx);
        status |= gr_div_ui(d, d, 32, ctx);
        status |= gr_sub(c, c, d, ctx);
        status |= gr_poly_set_coeff_scalar(f + 4, 0, c, ctx);
        GR_TMP_CLEAR(d, ctx);
    }

    for (n = 5; n <= N && status == GR_SUCCESS; n++)
    {
        slong m = n / 2;
        if (n % 2 == 1)
        {
            /* y^4 = F^2 appears in the term whose factors are even */
            status = gr_poly_pow_ui(t, f + m, 3, ctx);
            status |= gr_poly_mul(t, t, f + m + 2, ctx);
            status |= gr_poly_pow_ui(u, f + m + 1, 3, ctx);
            status |= gr_poly_mul(u, u, f + m - 1, ctx);
            if (m % 2 == 0)
            {
                /* psi_(m+2), psi_m even: y^4 in the first term */
                status |= gr_poly_mul(t, t, F, ctx);
                status |= gr_poly_mul(t, t, F, ctx);
            }
            else
            {
                status |= gr_poly_mul(u, u, F, ctx);
                status |= gr_poly_mul(u, u, F, ctx);
            }
            status |= gr_poly_sub(f + n, t, u, ctx);
        }
        else
        {
            /* psi_(m+2) psi_(m-1)^2 - psi_(m-2) psi_(m+1)^2: the y-powers
               of the terms are (m even) 1 + 0 and 1 + 0, (m odd) 0 + 2
               and 0 + 2; then times psi_m, divided by psi_2 = -y */
            status = gr_poly_mul(t, f + m - 1, f + m - 1, ctx);
            status |= gr_poly_mul(t, t, f + m + 2, ctx);
            status |= gr_poly_mul(u, f + m + 1, f + m + 1, ctx);
            status |= gr_poly_mul(u, u, f + m - 2, ctx);
            status |= gr_poly_sub(t, t, u, ctx);
            status |= gr_poly_mul(t, t, f + m, ctx);
            /* m even: y (inner) * y (psi_m) / (-y) = -y; m odd: y^2
               (inner) / (-y) = -y, psi_m odd: in both cases psi_n =
               -y f_m (inner without its y-powers) */
            status |= gr_poly_neg(f + n, t, ctx);
        }
    }

    for (i = 0; i <= N && status == GR_SUCCESS; i++)
        status = gr_poly_set(res + i, f + i, ctx);

    for (i = 0; i < TORSION_LEVEL + 3; i++)
        gr_poly_clear(f + i, ctx);
    gr_poly_clear(F, ctx);
    gr_poly_clear(t, ctx);
    gr_poly_clear(u, ctx);
    GR_TMP_CLEAR(c, ctx);
    return status;
}

/* z0 = (a + b tau0) / N with N <= TORSION_LEVEL (exactly): 1 */
static int
_mod_torsion_point(slong * N, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr x, y, t;
    fmpq_t qx, qy;
    int res = 0;

    GR_TMP_INIT3(x, y, t, ctx);
    fmpq_init(qx);
    fmpq_init(qy);

    if (gr_im(y, z0, ctx) == GR_SUCCESS && gr_im(t, tau0, ctx) == GR_SUCCESS &&
        gr_div(y, y, t, ctx) == GR_SUCCESS && gr_re(t, tau0, ctx) == GR_SUCCESS &&
        gr_mul(t, t, y, ctx) == GR_SUCCESS && gr_re(x, z0, ctx) == GR_SUCCESS &&
        gr_sub(x, x, t, ctx) == GR_SUCCESS &&
        _gr_tower_lazy_rational_recognize(qx, x, ctx) == 1 &&
        _gr_tower_lazy_rational_recognize(qy, y, ctx) == 1)
    {
        fmpz_t l;
        fmpz_init(l);
        fmpz_lcm(l, fmpq_denref(qx), fmpq_denref(qy));
        if (fmpz_cmp_ui(l, TORSION_LEVEL) <= 0 && fmpz_cmp_ui(l, 3) >= 0)
        {
            *N = fmpz_get_si(l);
            res = 1;
        }
        fmpz_clear(l);
    }

    GR_TMP_CLEAR3(x, y, t, ctx);
    fmpq_clear(qx);
    fmpq_clear(qy);
    return res;
}

/*
    The values at w from those at P = n w (n <= TORSION_LEVEL), or (thP =
    NULL) for P = n w a point of the lattice (w a torsion point): algebraic
    over the values at P and the theta constants. In the units of theta_3
    (rho_k = theta_k / theta_3, c_k = theta_k(0) / theta_3, wp~ = wp /
    (pi^2 theta_3^4), with e1 = (1 + c4^4)/3, e2 = (c2^4 - c4^4)/3, e3 =
    -(c2^4 + 1)/3, g2 = 2 (e1^2 + e2^2 + e3^2), g3 = 4 e1 e2 e3 in these
    units), X = wp~(w) is the root of the division polynomial f_n (P in
    the lattice) or of (x - wp~(P)) psi_n^2 - psi_(n-1) psi_(n+1) (from
    wp(n w) = wp(w) - psi_(n-1) psi_(n+1) / psi_n^2, of degree n^2), with
    wp~(P) = e1 + (c4 theta_2(P) / theta_1(P))^2; then the ratios
    r_k = theta_k(w) / theta_1(w),

        r_4^2 = (X - e3) / c2^2,  r_2^2 = (X - e1) / c4^2,  r_3^2 = (X - e2) / (c2 c4)^2,

    and with S_k(n) of _mod_mult_ratios, rho_1(w)^(n^2) = rho_k(P) / S_k(n)
    (all the roots chosen numerically, against acb).
*/
static int
_mod_theta_divide(gr_ptr * th, slong n, gr_srcptr w, gr_ptr * thP, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr c2, c4, e1, e2, e3, g2, g3, X, rv[4], S[4], t, u, W, rho1, t3, P, f[4], z0p, pw;
    gr_poly_struct D[TORSION_LEVEL + 2];
    gr_poly_t E, B;
    acb_t az, at, w1, w2, w3, w4, ref, q1, q2, q3, q4, zz;
    int kk[4], sg[4], pair, kref, k, status = GR_SUCCESS;
    slong i;

    if (n < 2 || n > TORSION_LEVEL)
        return GR_UNABLE;

    GR_TMP_INIT5(c2, c4, e1, e2, e3, ctx);
    GR_TMP_INIT5(g2, g3, X, t, u, ctx);
    GR_TMP_INIT5(W, rho1, t3, P, z0p, ctx);
    GR_TMP_INIT(pw, ctx);
    GR_TMP_INIT4(rv[0], rv[1], rv[2], rv[3], ctx);
    GR_TMP_INIT4(S[0], S[1], S[2], S[3], ctx);
    GR_TMP_INIT4(f[0], f[1], f[2], f[3], ctx);
    for (i = 0; i < TORSION_LEVEL + 2; i++)
        gr_poly_init(D + i, ctx);
    gr_poly_init(E, ctx);
    gr_poly_init(B, ctx);
    acb_init(az); acb_init(at); acb_init(w1); acb_init(w2); acb_init(w3); acb_init(w4);
    acb_init(ref); acb_init(q1); acb_init(q2); acb_init(q3); acb_init(q4); acb_init(zz);

    /* numerically: the ratios at w, theta_3 */
    status = gr_tower_lazy_get_acb(az, w, 128, ctx);
    status |= gr_tower_lazy_get_acb(at, tau0, 128, ctx);
    if (status == GR_SUCCESS)
    {
        acb_modular_theta(w1, w2, w3, w4, az, at, 128);
        acb_div(w2, w2, w1, 128);
        acb_div(w3, w3, w1, 128);
        acb_div(w4, w4, w1, 128);
        acb_zero(zz);
        acb_modular_theta(q1, q2, q3, q4, zz, at, 128);
    }

    /* the constants, in the units of theta_3 */
    if (status == GR_SUCCESS)
        status = _mod_point(c2, NULL, c4, NULL, NULL, tau0, 1, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_pow_ui(t, c2, 4, ctx);
        status |= gr_pow_ui(u, c4, 4, ctx);
        status |= gr_add_ui(e1, u, 1, ctx);
        status |= gr_div_ui(e1, e1, 3, ctx);
        status |= gr_sub(e2, t, u, ctx);
        status |= gr_div_ui(e2, e2, 3, ctx);
        status |= gr_add_ui(e3, t, 1, ctx);
        status |= gr_div_si(e3, e3, -3, ctx);
        status |= gr_sqr(g2, e1, ctx);
        status |= gr_sqr(t, e2, ctx);
        status |= gr_add(g2, g2, t, ctx);
        status |= gr_sqr(t, e3, ctx);
        status |= gr_add(g2, g2, t, ctx);
        status |= gr_mul_ui(g2, g2, 2, ctx);
        status |= gr_mul(g3, e1, e2, ctx);
        status |= gr_mul(g3, g3, e3, ctx);
        status |= gr_mul_ui(g3, g3, 4, ctx);
    }

    /* X = wp~(w) */
    if (status == GR_SUCCESS)
        status = _mod_division_polys(D, n + 1, g2, g3, ctx);
    if (status == GR_SUCCESS)
    {
        if (thP == NULL)
            status = gr_poly_set(E, D + n, ctx);
        else
        {
            /* (x - p) F^(e_n) f_n^2 - F^(e_(n+1)) f_(n-1) f_(n+1), F = 4x^3 - g2 x - g3 */
            gr_poly_t F;
            gr_poly_init(F, ctx);
            status = gr_poly_set_coeff_si(F, 3, 4, ctx);
            status |= gr_neg(t, g2, ctx);
            status |= gr_poly_set_coeff_scalar(F, 1, t, ctx);
            status |= gr_neg(t, g3, ctx);
            status |= gr_poly_set_coeff_scalar(F, 0, t, ctx);
            /* p = e1 + (c4 theta_2(P) / theta_1(P))^2 */
            status |= gr_div(pw, thP[1], thP[0], ctx);
            status |= gr_mul(pw, pw, c4, ctx);
            status |= gr_sqr(pw, pw, ctx);
            status |= gr_add(pw, pw, e1, ctx);
            status |= gr_poly_mul(E, D + n, D + n, ctx);
            if (n % 2 == 0)
                status |= gr_poly_mul(E, E, F, ctx);
            status |= gr_poly_zero(B, ctx);
            status |= gr_poly_set_coeff_si(B, 1, 1, ctx);
            status |= gr_neg(t, pw, ctx);
            status |= gr_poly_set_coeff_scalar(B, 0, t, ctx);
            status |= gr_poly_mul(E, E, B, ctx);
            status |= gr_poly_mul(B, D + n - 1, D + n + 1, ctx);
            if (n % 2 == 1)
                status |= gr_poly_mul(B, B, F, ctx);
            status |= gr_poly_sub(E, E, B, ctx);
            gr_poly_clear(F, ctx);
        }
    }
    if (status == GR_SUCCESS)
    {
        acb_elliptic_p(ref, az, at, 128);
        acb_pow_ui(q1, q3, 4, 128);
        acb_const_pi(q2, 128);
        acb_mul(q2, q2, q2, 128);
        acb_mul(q1, q1, q2, 128);
        acb_div(ref, ref, q1, 128);

        status = _gr_tower_lazy_poly_root_near(X, E, ref, 0, ctx);
    }

    /* the ratios r_k = theta_k(w) / theta_1(w) */
    if (status == GR_SUCCESS)
    {
        status = gr_one(rv[0], ctx);
        status |= gr_sub(t, X, e1, ctx);
        status |= gr_sqr(u, c4, ctx);
        status |= gr_div(t, t, u, ctx);
        if (status == GR_SUCCESS)
            status = _mod_sqrt_near(rv[1], t, w2, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = gr_sub(t, X, e2, ctx);
        status |= gr_mul(u, c2, c4, ctx);
        status |= gr_sqr(u, u, ctx);
        status |= gr_div(t, t, u, ctx);
        if (status == GR_SUCCESS)
            status = _mod_sqrt_near(rv[2], t, w3, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = gr_sub(t, X, e3, ctx);
        status |= gr_sqr(u, c2, ctx);
        status |= gr_div(t, t, u, ctx);
        if (status == GR_SUCCESS)
            status = _mod_sqrt_near(rv[3], t, w4, ctx);
    }

    /* rho_k(P): from the lattice point (sg f c_kk), or theta_k(P) / theta_3 */
    if (status == GR_SUCCESS)
        status = _mod_theta0(NULL, t3, NULL, NULL, tau0, ctx);
    pair = 0;
    kref = 1;
    if (status == GR_SUCCESS && thP == NULL)
    {
        status = gr_mul_si(P, w, n, ctx);
        if (status == GR_SUCCESS)
            status = _mod_reduce_map(z0p, f, kk, sg, P, tau0, ctx);
        if (status == GR_SUCCESS && n % 2 == 0)
        {
            /* for even n, the pair avoiding the half period (n/2) w:
               theta_4 vanishes there when it maps to theta_1 */
            gr_ptr hf[4], hz;
            int hkk[4], hsg[4];
            GR_TMP_INIT4(hf[0], hf[1], hf[2], hf[3], ctx);
            GR_TMP_INIT(hz, ctx);
            status = gr_mul_si(t, w, n / 2, ctx);
            if (status == GR_SUCCESS)
                status = _mod_reduce_map(hz, hf, hkk, hsg, t, tau0, ctx);
            if (status == GR_SUCCESS && hkk[3] == 1)
                pair = 1;
            GR_TMP_CLEAR4(hf[0], hf[1], hf[2], hf[3], ctx);
            GR_TMP_CLEAR(hz, ctx);
        }
        kref = (pair == 0) ? 4 : 2;
        if (status == GR_SUCCESS)
        {
            k = kk[kref - 1];
            if (k == 1)
                status = GR_UNABLE;
            else if (k == 2)
                status = gr_set(t, c2, ctx);
            else if (k == 4)
                status = gr_set(t, c4, ctx);
            else
                status = gr_one(t, ctx);
            status |= gr_mul(t, t, f[kref - 1], ctx);
            if (sg[kref - 1] < 0)
                status |= gr_neg(t, t, ctx);
        }
    }
    else if (status == GR_SUCCESS)
    {
        status = gr_div(t, thP[0], t3, ctx);
    }

    if (status == GR_SUCCESS)
        status = _mod_mult_ratios(S, pair, (gr_srcptr *) rv, c2, c4, n, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(W, t, S[kref - 1], ctx);

    /* rho_1(w): the root of x^(n^2) - W near the value */
    if (status == GR_SUCCESS)
    {
        acb_div(ref, w1, q3, 128);
        status |= gr_poly_zero(B, ctx);
        status = gr_poly_set_coeff_si(B, n * n, 1, ctx);
        status |= gr_neg(t, W, ctx);
        status |= gr_poly_set_coeff_scalar(B, 0, t, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_poly_root_near(rho1, B, ref, 0, ctx);
    }

    /* theta_k(w) = theta_3 rho_1(w) r_k */
    for (k = 0; k < 4 && status == GR_SUCCESS; k++)
    {
        status = gr_mul(th[k], t3, rho1, ctx);
        status |= gr_mul(th[k], th[k], rv[k], ctx);
    }

    GR_TMP_CLEAR5(c2, c4, e1, e2, e3, ctx);
    GR_TMP_CLEAR5(g2, g3, X, t, u, ctx);
    GR_TMP_CLEAR5(W, rho1, t3, P, z0p, ctx);
    GR_TMP_CLEAR(pw, ctx);
    GR_TMP_CLEAR4(rv[0], rv[1], rv[2], rv[3], ctx);
    GR_TMP_CLEAR4(S[0], S[1], S[2], S[3], ctx);
    GR_TMP_CLEAR4(f[0], f[1], f[2], f[3], ctx);
    for (i = 0; i < TORSION_LEVEL + 2; i++)
        gr_poly_clear(D + i, ctx);
    gr_poly_clear(E, ctx);
    gr_poly_clear(B, ctx);
    acb_clear(az); acb_clear(at); acb_clear(w1); acb_clear(w2); acb_clear(w3); acb_clear(w4);
    acb_clear(ref); acb_clear(q1); acb_clear(q2); acb_clear(q3); acb_clear(q4); acb_clear(zz);
    return status;
}

#define THETA_MULT 6
#define THETA_DIV 3


/* whether 2 (z - w) is in Z + tau0 Z numerically (z = w modulo the half
   lattice), with x, y the coordinates of 2 (z - w) */
static int
_mod_half_lattice_num(const acb_t z, const acb_t w, const acb_t tau)
{
    acb_t d;
    arb_t y, x, t;
    int ok = 0;

    acb_init(d);
    arb_init(y);
    arb_init(x);
    arb_init(t);
    acb_sub(d, z, w, 128);
    acb_mul_2exp_si(d, d, 1);
    arb_div(y, acb_imagref(d), acb_imagref(tau), 128);
    arb_mul(t, y, acb_realref(tau), 128);
    arb_sub(x, acb_realref(d), t, 128);
    {
        fmpz_t n;
        fmpz_init(n);
        if (arf_cmpabs_2exp_si(arb_midref(x), 30) < 0 && arf_cmpabs_2exp_si(arb_midref(y), 30) < 0)
        {
            arf_get_fmpz(n, arb_midref(x), ARF_RND_NEAR);
            ok = arb_contains_fmpz(x, n);
            arf_get_fmpz(n, arb_midref(y), ARF_RND_NEAR);
            ok = ok && arb_contains_fmpz(y, n);
        }
        fmpz_clear(n);
    }
    acb_clear(d);
    arb_clear(y);
    arb_clear(x);
    arb_clear(t);
    return ok;
}

/* the same, exactly */
static int
_mod_half_lattice(gr_srcptr z, gr_srcptr w, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr d, x, y, t;
    slong n;
    int ok = 0;

    GR_TMP_INIT4(d, x, y, t, ctx);
    if (gr_sub(d, z, w, ctx) == GR_SUCCESS && gr_mul_ui(d, d, 2, ctx) == GR_SUCCESS &&
        gr_im(y, d, ctx) == GR_SUCCESS && gr_im(t, tau0, ctx) == GR_SUCCESS &&
        gr_div(y, y, t, ctx) == GR_SUCCESS && gr_re(t, tau0, ctx) == GR_SUCCESS &&
        gr_mul(t, t, y, ctx) == GR_SUCCESS && gr_re(x, d, ctx) == GR_SUCCESS &&
        gr_sub(x, x, t, ctx) == GR_SUCCESS)
    {
        {
            fmpq_t qa, qb;
            fmpq_init(qa);
            fmpq_init(qb);
            ok = _gr_tower_lazy_rational_recognize(qa, x, ctx) == 1 && fmpz_is_one(fmpq_denref(qa)) &&
                 _gr_tower_lazy_rational_recognize(qb, y, ctx) == 1 && fmpz_is_one(fmpq_denref(qb));
            fmpq_clear(qa);
            fmpq_clear(qb);
            (void) n;
        }
    }
    GR_TMP_CLEAR4(d, x, y, t, ctx);
    return ok;
}

/* val_k = theta_k(n w) from tw_k = theta_k(w): theta_3 (theta_1(w) / theta_3)^(n^2) S_k(n) */
static int
_mod_theta_mult_vals(gr_ptr * val, gr_ptr * tw, slong n, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr rv[4], S[4], c2, c4, t3, t;
    int k, status;

    GR_TMP_INIT4(rv[0], rv[1], rv[2], rv[3], ctx);
    GR_TMP_INIT4(S[0], S[1], S[2], S[3], ctx);
    GR_TMP_INIT4(c2, c4, t3, t, ctx);

    status = _mod_point(c2, NULL, c4, NULL, NULL, tau0, 1, ctx);
    status |= _mod_theta0(NULL, t3, NULL, NULL, tau0, ctx);
    for (k = 0; k < 4 && status == GR_SUCCESS; k++)
        status = gr_div(rv[k], tw[k], tw[0], ctx);
    if (status == GR_SUCCESS)
        status = _mod_mult_ratios(S, 0, (gr_srcptr *) rv, c2, c4, n, ctx);
    if (status == GR_SUCCESS)
        status = _mod_mult_ratios(S, 1, (gr_srcptr *) rv, c2, c4, n, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_div(t, tw[0], t3, ctx);
        status |= gr_pow_ui(t, t, n * n, ctx);
        status |= gr_mul(t, t, t3, ctx);
        for (k = 0; k < 4; k++)
            status |= gr_mul(val[k], t, S[k], ctx);
    }

    GR_TMP_CLEAR4(rv[0], rv[1], rv[2], rv[3], ctx);
    GR_TMP_CLEAR4(S[0], S[1], S[2], S[3], ctx);
    GR_TMP_CLEAR4(c2, c4, t3, t, ctx);
    return status;
}


/*
    The points of the theta functions of z. Every reduced point z0 != 0
    at which the theta functions are evaluated is recorded in the context
    with the method of its evaluation, chosen when it is first met (from
    the points recorded then), so that later evaluations at z0 give the
    same representation:

    TP_STANDARD   theta_1(z0), theta_4(z0) generators, theta_2(z0),
                  theta_3(z0) algebraic over them (_mod_theta_z0_gens);
    TP_TORSION    z0 = (a + b tau0) / N: algebraic over the theta
                  constants (_mod_theta_divide with the lattice point N z0);
    TP_MULT       dv z0 = s n u modulo the half lattice, u a recorded point:
                  the multiplication formulas from u to s n u, then (dv > 1)
                  the division by dv: algebraic over the values at u;
    TP_ADD        z0 = s (u + s2 v) with u - s2 v recorded too: the addition
                  formulas theta_k(u + v) theta_k(u - v) theta_4^2 = A_k(u, v),
                  rational;
    TP_ADD_HALF   z0 with u, u + z0 and u - z0 recorded: theta_1(z0)^2 and
                  theta_4(z0)^2 from the addition formulas for theta_1(u +
                  z0) theta_1(u - z0) and theta_4(u + z0) theta_4(u - z0),
                  linear in them, then theta_2(z0)^2, theta_3(z0)^2 (four
                  square roots); when u is of type TP_HALF, u = (P + Q) / 2,
                  and z0 = +-(P - Q) / 2: monomials in the generators of u
                  (_mod_tp_add_half_sym);
    TP_JACOBI     z0 = s (n u + s2 v), n = 1, 2, 3 (for n = 1, without the
                  partner of TP_ADD): the ratios of the theta functions at
                  z0 rational in the values at u and v (Jacobi's addition
                  formulas for sn, cn, dn, from the multiple n u),
                  theta_1(z0) a generator. (With n = 2: the values at z,
                  z + w, z - w are not independent, so that those at w
                  then come from TP_ADD_HALF consistently with all the
                  addition formulas.)
    TP_HALF       2 z0 = s (u + s2 v) otherwise, z0 = (P + Q) / 2 for P, Q
                  the recorded points modulo the half lattice: the ratios
                  monomials in the square roots of the forms S_ab =
                  theta_a(P) theta_b(Q) + theta_b(P) theta_a(Q)
                  (_mod_half_sum_ratios; else by the halving of P + Q),
                  theta_1(z0) a generator. (So with z + w and z - w
                  evaluated first, z is a half of their sum and w then of
                  type TP_ADD_HALF: the values at the four points satisfy
                  the addition formulas as rational identities in the S_ab
                  and D_ab = theta_a(P) theta_b(Q) - theta_b(P) theta_a(Q).)

    (Equalities modulo the half lattice and up to sign are those of the
    reduction: _mod_reduce_map relates the values.)
*/
#define TP_STANDARD 0
#define TP_TORSION 1
#define TP_MULT 2
#define TP_ADD 3
#define TP_ADD_HALF 4
#define TP_JACOBI 5
#define TP_HALF 6

/* at most this many recorded points at one tau0 are searched */
#define THETA_POINTS 32

/* the recorded point equal to (z0, tau0), or -1 */
static slong
_mod_tp_find(gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    acb_t a, b, c, d;
    slong k, res = -1;

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(d);
    if (gr_tower_lazy_get_acb(a, z0, 128, ctx) == GR_SUCCESS &&
        gr_tower_lazy_get_acb(c, tau0, 128, ctx) == GR_SUCCESS)
    {
        for (k = 0; k < L->num_theta_pts && res < 0; k++)
        {
            gr_tower_lazy_theta_point_struct * e = L->theta_pts + k;
            if (gr_tower_lazy_get_acb(b, &e->z0, 128, ctx) != GR_SUCCESS || !acb_overlaps(a, b))
                continue;
            if (gr_tower_lazy_get_acb(d, &e->tau0, 128, ctx) != GR_SUCCESS || !acb_overlaps(c, d))
                continue;
            if (gr_equal(z0, &e->z0, ctx) == T_TRUE && gr_equal(tau0, &e->tau0, ctx) == T_TRUE)
                res = k;
        }
    }
    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(d);
    return res;
}

/* the recorded points at tau0: their indices (at most THETA_POINTS) */
static slong
_mod_tp_points(slong * idx, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    acb_t c, d;
    slong k, n = 0;

    acb_init(c);
    acb_init(d);
    if (gr_tower_lazy_get_acb(c, tau0, 128, ctx) == GR_SUCCESS)
    {
        for (k = 0; k < L->num_theta_pts && n < THETA_POINTS; k++)
        {
            gr_tower_lazy_theta_point_struct * e = L->theta_pts + k;
            if (gr_tower_lazy_get_acb(d, &e->tau0, 128, ctx) != GR_SUCCESS || !acb_overlaps(c, d))
                continue;
            if (gr_equal(tau0, &e->tau0, ctx) == T_TRUE)
                idx[n++] = k;
        }
    }
    acb_clear(c);
    acb_clear(d);
    return n;
}

/* records (z0, tau0) with the method of e (copied); returns its index */
static slong
_mod_tp_record(gr_srcptr z0, gr_srcptr tau0, const gr_tower_lazy_theta_point_struct * e, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_lazy_theta_point_struct * f;

    if (L->num_theta_pts == L->alloc_theta_pts)
    {
        slong k, alloc = FLINT_MAX(8, 2 * L->alloc_theta_pts);
        gr_tower_lazy_theta_point_struct * p = flint_malloc(sizeof(gr_tower_lazy_theta_point_struct) * alloc);
        /* (the elements are moved: their bits are copied) */
        for (k = 0; k < L->num_theta_pts; k++)
            p[k] = L->theta_pts[k];
        flint_free(L->theta_pts);
        L->theta_pts = p;
        L->alloc_theta_pts = alloc;
    }

    f = L->theta_pts + L->num_theta_pts;
    *f = *e;
    f->have_vals = 0;
    f->stamp = 0;
    f->pinned = 0;
    _gr_tower_lazy_init(&f->z0, ctx);
    _gr_tower_lazy_init(&f->tau0, ctx);
    GR_MUST_SUCCEED(gr_set(&f->z0, z0, ctx));
    GR_MUST_SUCCEED(gr_set(&f->tau0, tau0, ctx));
    return L->num_theta_pts++;
}

/* stores the values of the record k */
static void
_mod_tp_store(slong k, gr_ptr * th, gr_ctx_t ctx)
{
    gr_tower_lazy_theta_point_struct * e = LAZY(ctx)->theta_pts + k;
    slong q;
    if (e->have_vals)
        return;
    for (q = 0; q < 4; q++)
    {
        _gr_tower_lazy_init(&e->vals[q], ctx);
        GR_MUST_SUCCEED(gr_set(&e->vals[q], th[q], ctx));
    }
    e->have_vals = 1;
    e->stamp = ++LAZY(ctx)->cache_clock;
    _mod_tp_evict(k, ctx);
}

/* the point of the record k (a copy) */
static int
_mod_tp_z(gr_ptr res, slong k, gr_ctx_t ctx)
{
    return gr_set(res, &LAZY(ctx)->theta_pts[k].z0, ctx);
}

/*
    The method for the new point z0 (not recorded) from the recorded
    points at tau0 (the search by enclosures, each candidate verified
    exactly): sets e (kind TP_STANDARD if none applies).
*/
static int _mod_torsion_field_ok(gr_srcptr tau0, gr_ctx_t ctx);

static int
_mod_tp_search(gr_tower_lazy_theta_point_struct * e, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    slong idx[THETA_POINTS], np, i, j, n, l, dv, maxdv, maxn;
    gr_ptr pts, P, Q, t;
    acb_ptr ap;
    acb_t az, at, w, dz, pp, qq, m;
    int status = GR_SUCCESS, found = 0;
    slong k5i = -1, k5j = -1, k5m = -1, k5n = 1, sz = ctx->sizeof_elem;
    int k5s = 1, k5s2 = 1;

    e->kind = TP_STANDARD;
    e->i = e->j = e->l = -1;
    e->n = e->dv = 0;
    e->s = e->s2 = 1;

    np = _mod_tp_points(idx, tau0, ctx);
    if (np == 0)
        return GR_SUCCESS;

    pts = gr_heap_init_vec(np, ctx);
    GR_TMP_INIT3(P, Q, t, ctx);
    ap = _acb_vec_init(np);
    acb_init(az); acb_init(at); acb_init(w); acb_init(dz);
    acb_init(pp); acb_init(qq); acb_init(m);

    status = gr_tower_lazy_get_acb(az, z0, 128, ctx);
    status |= gr_tower_lazy_get_acb(at, tau0, 128, ctx);
    for (i = 0; i < np && status == GR_SUCCESS; i++)
    {
        status = _mod_tp_z(GR_ENTRY(pts, i, sz), idx[i], ctx);
        status |= gr_tower_lazy_get_acb(ap + i, GR_ENTRY(pts, i, sz), 128, ctx);
    }

    /* rational multiples: dv z0 = s n u_i modulo the half lattice (the
       divisions dv > 1 and the multiples n > 2 not at the CM points of
       _mod_torsion_field_ok, where lambda is algebraic of large degree
       and the formulas of the other relations, with the inverses of
       their denominators in that field, would take minutes, or hours
       for n = 6; the multiples of all the points before any division) */
    maxdv = _mod_torsion_field_ok(tau0, ctx) ? THETA_DIV : 1;
    maxn = (maxdv > 1) ? THETA_MULT : 2;
    for (dv = 1; dv <= maxdv && !found; dv++)
    {
        for (i = 0; i < np && status == GR_SUCCESS && !found; i++)
        {
            for (n = 1; n <= maxn && !found; n++)
            {
                int s;
                if ((n == 1 && dv == 1) || n_gcd(n, dv) != 1)
                    continue;
                for (s = -1; s <= 1 && !found; s += 2)
                {
                    acb_mul_si(w, ap + i, s * n, 128);
                    acb_mul_si(dz, az, dv, 128);
                    if (!_mod_half_lattice_num(dz, w, at))
                        continue;
                    status = gr_mul_si(P, GR_ENTRY(pts, i, sz), s * n, ctx);
                    status |= gr_mul_si(Q, z0, dv, ctx);
                    if (status == GR_SUCCESS && _mod_half_lattice(Q, P, tau0, ctx))
                    {
                        found = 1;
                        e->kind = TP_MULT; e->i = idx[i]; e->n = n; e->dv = dv; e->s = s;
                    }
                }
            }
        }
    }

    /* (at the CM points of large degree, the other relations as well:
       independent generators at z0) */
    if (maxdv == 1)
        goto done;

    /* sums and differences: z0 = s (u_i + s2 u_j) */
    for (i = 0; i < np && status == GR_SUCCESS && !found; i++)
    {
        for (j = 0; j < np && status == GR_SUCCESS && !found; j++)
        {
            int s, s2;
            if (j == i)
                continue;
            for (s = -1; s <= 1 && !found; s += 2)
            {
                for (s2 = -1; s2 <= 1 && !found; s2 += 2)
                {
                    acb_mul_si(pp, ap + j, s2, 128);
                    acb_add(pp, ap + i, pp, 128);
                    acb_mul_si(pp, pp, s, 128);
                    if (!_mod_half_lattice_num(az, pp, at))
                        continue;
                    status = gr_mul_si(P, GR_ENTRY(pts, j, sz), s2, ctx);
                    status |= gr_add(P, GR_ENTRY(pts, i, sz), P, ctx);
                    status |= gr_mul_si(P, P, s, ctx);
                    if (status != GR_SUCCESS || !_mod_half_lattice(z0, P, tau0, ctx))
                        continue;
                    /* the partner u_i - s2 u_j */
                    acb_mul_si(qq, ap + j, s2, 128);
                    acb_sub(qq, ap + i, qq, 128);
                    for (l = 0; l < np && !found; l++)
                    {
                        acb_neg(m, ap + l);
                        if (_mod_half_lattice_num(qq, ap + l, at) || _mod_half_lattice_num(qq, m, at))
                        {
                            status = gr_mul_si(Q, GR_ENTRY(pts, j, sz), s2, ctx);
                            status |= gr_sub(Q, GR_ENTRY(pts, i, sz), Q, ctx);
                            status |= gr_neg(t, GR_ENTRY(pts, l, sz), ctx);
                            if (status == GR_SUCCESS &&
                                (_mod_half_lattice(Q, GR_ENTRY(pts, l, sz), tau0, ctx) ||
                                 _mod_half_lattice(Q, t, tau0, ctx)))
                            {
                                found = 1;
                                e->kind = TP_ADD; e->i = idx[i]; e->j = idx[j]; e->l = idx[l];
                                e->s = s; e->s2 = s2;
                            }
                        }
                    }
                    if (!found && (k5i < 0 || FLINT_MIN(i, j) > k5m))
                    {
                        k5i = idx[i]; k5j = idx[j]; k5s = s; k5s2 = s2;
                        k5m = FLINT_MIN(i, j); k5n = 1;
                    }
                }
            }
        }
    }

    /* z0 with u_i, u_i + z0 and u_i - z0 recorded */
    for (i = 0; i < np && status == GR_SUCCESS && !found; i++)
    {
        slong l1, l2;
        acb_add(pp, ap + i, az, 128);
        acb_sub(qq, ap + i, az, 128);
        for (l1 = 0; l1 < np && !found && status == GR_SUCCESS; l1++)
        {
            acb_neg(m, ap + l1);
            if (!(_mod_half_lattice_num(pp, ap + l1, at) || _mod_half_lattice_num(pp, m, at)))
                continue;
            for (l2 = 0; l2 < np && !found && status == GR_SUCCESS; l2++)
            {
                int ok;
                acb_neg(m, ap + l2);
                if (!(_mod_half_lattice_num(qq, ap + l2, at) || _mod_half_lattice_num(qq, m, at)))
                    continue;
                status = gr_add(P, GR_ENTRY(pts, i, sz), z0, ctx);
                status |= gr_sub(Q, GR_ENTRY(pts, i, sz), z0, ctx);
                if (status != GR_SUCCESS)
                    break;
                status = gr_neg(t, GR_ENTRY(pts, l1, sz), ctx);
                ok = _mod_half_lattice(P, GR_ENTRY(pts, l1, sz), tau0, ctx) || _mod_half_lattice(P, t, tau0, ctx);
                status |= gr_neg(t, GR_ENTRY(pts, l2, sz), ctx);
                ok = ok && (_mod_half_lattice(Q, GR_ENTRY(pts, l2, sz), tau0, ctx) || _mod_half_lattice(Q, t, tau0, ctx));
                if (ok)
                {
                    found = 1;
                    e->kind = TP_ADD_HALF; e->i = idx[i];
                }
            }
        }
    }

    /* z0 = s (n u_i + s2 u_j) with n = 2, 3 (Jacobi's formulas from the
       multiple: z0 = 2u - v with u and v = u + w recorded, say, whose
       values would otherwise be independent of those at z0); among the
       candidates (with n = 1 above), those with the latest points, from
       which the earlier ones may be derived (2u from u, say, rather than
       u by division from 2u, which would make the relations at z0 hold
       only modulo the division) */
    for (n = 2; n <= 3 && !found && status == GR_SUCCESS; n++)
    {
        for (i = 0; i < np && status == GR_SUCCESS; i++)
        {
            for (j = 0; j < np && status == GR_SUCCESS; j++)
            {
                int s, s2;
                if (j == i || FLINT_MIN(i, j) <= k5m)
                    continue;
                for (s = -1; s <= 1; s += 2)
                {
                    for (s2 = -1; s2 <= 1; s2 += 2)
                    {
                        if (FLINT_MIN(i, j) <= k5m)
                            continue;
                        acb_mul_si(pp, ap + i, n, 128);
                        acb_mul_si(qq, ap + j, s2, 128);
                        acb_add(pp, pp, qq, 128);
                        acb_mul_si(pp, pp, s, 128);
                        if (!_mod_half_lattice_num(az, pp, at))
                            continue;
                        status = gr_mul_si(P, GR_ENTRY(pts, i, sz), n, ctx);
                        status |= gr_mul_si(Q, GR_ENTRY(pts, j, sz), s2, ctx);
                        status |= gr_add(P, P, Q, ctx);
                        status |= gr_mul_si(P, P, s, ctx);
                        if (status == GR_SUCCESS && _mod_half_lattice(z0, P, tau0, ctx))
                        {
                            k5i = idx[i]; k5j = idx[j]; k5s = s; k5s2 = s2;
                            k5m = FLINT_MIN(i, j); k5n = n;
                        }
                    }
                }
            }
        }
    }

    if (!found && k5i >= 0 && status == GR_SUCCESS)
    {
        found = 1;
        e->kind = TP_JACOBI; e->i = k5i; e->j = k5j; e->s = k5s; e->s2 = k5s2;
        e->n = k5n; e->dv = 1;
    }

    /* 2 z0 = s (u_i + s2 u_j): the half of a sum or difference (z with
       z + w and z - w recorded) */
    for (i = 0; i < np && status == GR_SUCCESS && !found; i++)
    {
        for (j = 0; j < np && status == GR_SUCCESS && !found; j++)
        {
            int s, s2;
            if (j == i)
                continue;
            for (s = -1; s <= 1 && !found; s += 2)
            {
                for (s2 = -1; s2 <= 1 && !found; s2 += 2)
                {
                    acb_mul_si(qq, ap + j, s2, 128);
                    acb_add(pp, ap + i, qq, 128);
                    acb_mul_si(pp, pp, s, 128);
                    acb_mul_2exp_si(dz, az, 1);
                    if (!_mod_half_lattice_num(dz, pp, at))
                        continue;
                    status = gr_mul_si(Q, GR_ENTRY(pts, j, sz), s2, ctx);
                    status |= gr_add(P, GR_ENTRY(pts, i, sz), Q, ctx);
                    status |= gr_mul_si(P, P, s, ctx);
                    status |= gr_mul_ui(Q, z0, 2, ctx);
                    if (status == GR_SUCCESS && _mod_half_lattice(Q, P, tau0, ctx))
                    {
                        found = 1;
                        e->kind = TP_HALF; e->i = idx[i]; e->j = idx[j]; e->s = s; e->s2 = s2;
                    }
                }
            }
        }
    }

done:
    gr_heap_clear_vec(pts, np, ctx);
    GR_TMP_CLEAR3(P, Q, t, ctx);
    _acb_vec_clear(ap, np);
    acb_clear(az); acb_clear(at); acb_clear(w); acb_clear(dz);
    acb_clear(pp); acb_clear(qq); acb_clear(m);

    if (status != GR_SUCCESS)
        e->kind = TP_STANDARD;
    return GR_SUCCESS;
}

/* th[kk - 1] = val[k - 1] / (sg f), the values at the reduced point of P */
static int
_mod_tp_unmap(gr_ptr * th, gr_ptr * val, gr_srcptr P, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr f[4], zr;
    int kk[4], sg[4], k, status;

    GR_TMP_INIT4(f[0], f[1], f[2], f[3], ctx);
    GR_TMP_INIT(zr, ctx);
    status = _mod_reduce_map(zr, f, kk, sg, P, tau0, ctx);
    for (k = 0; k < 4 && status == GR_SUCCESS; k++)
    {
        status = gr_div(th[kk[k] - 1], val[k], f[k], ctx);
        if (sg[k] < 0)
            status |= gr_neg(th[kk[k] - 1], th[kk[k] - 1], ctx);
    }
    GR_TMP_CLEAR4(f[0], f[1], f[2], f[3], ctx);
    GR_TMP_CLEAR(zr, ctx);
    return status;
}

/* val[k - 1] = theta_k(P) from th, the values at the reduced point of P */
static int
_mod_tp_map(gr_ptr * val, gr_ptr * th, gr_srcptr P, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr f[4], zr;
    int kk[4], sg[4], k, status;

    GR_TMP_INIT4(f[0], f[1], f[2], f[3], ctx);
    GR_TMP_INIT(zr, ctx);
    status = _mod_reduce_map(zr, f, kk, sg, P, tau0, ctx);
    for (k = 0; k < 4 && status == GR_SUCCESS; k++)
    {
        status = gr_mul(val[k], f[k], th[kk[k] - 1], ctx);
        if (sg[k] < 0)
            status |= gr_neg(val[k], val[k], ctx);
    }
    GR_TMP_CLEAR4(f[0], f[1], f[2], f[3], ctx);
    GR_TMP_CLEAR(zr, ctx);
    return status;
}

/*
    Jacobi's addition formulas: val_k = theta_k(u + v) / theta_4(u + v)
    from the values tu at u and tv at v (val may alias tv), with sn =
    theta_1 / (c2 theta_4), cn = c4 theta_2 / (c2 theta_4), dn = c4
    theta_3 / theta_4, k^2 = c2^4:
        sn(u + v) = (sn u cn v dn v + sn v cn u dn u) / D,
        cn(u + v) = (cn u cn v - sn u sn v dn u dn v) / D,
        dn(u + v) = (dn u dn v - k^2 sn u sn v cn u cn v) / D,
    D = 1 - k^2 sn(u)^2 sn(v)^2
*/
static int
_mod_jacobi_add(gr_ptr * val, gr_ptr * tu, gr_ptr * tv, gr_srcptr c2, gr_srcptr c4, gr_ctx_t ctx)
{
    gr_ptr snu, cnu, dnu, snv, cnv, dnv, den, k2, x;
    int status;

    GR_TMP_INIT5(snu, cnu, dnu, snv, cnv, ctx);
    GR_TMP_INIT4(dnv, den, k2, x, ctx);

    status = gr_pow_ui(k2, c2, 4, ctx);
    status |= gr_div(snu, tu[0], tu[3], ctx);
    status |= gr_div(snu, snu, c2, ctx);
    status |= gr_div(cnu, tu[1], tu[3], ctx);
    status |= gr_mul(cnu, cnu, c4, ctx);
    status |= gr_div(cnu, cnu, c2, ctx);
    status |= gr_div(dnu, tu[2], tu[3], ctx);
    status |= gr_mul(dnu, dnu, c4, ctx);
    status |= gr_div(snv, tv[0], tv[3], ctx);
    status |= gr_div(snv, snv, c2, ctx);
    status |= gr_div(cnv, tv[1], tv[3], ctx);
    status |= gr_mul(cnv, cnv, c4, ctx);
    status |= gr_div(cnv, cnv, c2, ctx);
    status |= gr_div(dnv, tv[2], tv[3], ctx);
    status |= gr_mul(dnv, dnv, c4, ctx);
    status |= gr_mul(den, snu, snv, ctx);
    status |= gr_sqr(den, den, ctx);
    status |= gr_mul(den, den, k2, ctx);
    status |= gr_sub_ui(den, den, 1, ctx);
    status |= gr_neg(den, den, ctx);
    /* sn */
    status |= gr_mul(val[0], snu, cnv, ctx);
    status |= gr_mul(val[0], val[0], dnv, ctx);
    status |= gr_mul(x, snv, cnu, ctx);
    status |= gr_mul(x, x, dnu, ctx);
    status |= gr_add(val[0], val[0], x, ctx);
    status |= gr_div(val[0], val[0], den, ctx);
    /* cn */
    status |= gr_mul(val[1], cnu, cnv, ctx);
    status |= gr_mul(x, snu, snv, ctx);
    status |= gr_mul(x, x, dnu, ctx);
    status |= gr_mul(x, x, dnv, ctx);
    status |= gr_sub(val[1], val[1], x, ctx);
    status |= gr_div(val[1], val[1], den, ctx);
    /* dn */
    status |= gr_mul(val[2], dnu, dnv, ctx);
    status |= gr_mul(x, snu, snv, ctx);
    status |= gr_mul(x, x, cnu, ctx);
    status |= gr_mul(x, x, cnv, ctx);
    status |= gr_mul(x, x, k2, ctx);
    status |= gr_sub(val[2], val[2], x, ctx);
    status |= gr_div(val[2], val[2], den, ctx);
    /* the ratios theta_k / theta_4 */
    status |= gr_mul(val[0], val[0], c2, ctx);
    status |= gr_mul(val[1], val[1], c2, ctx);
    status |= gr_div(val[1], val[1], c4, ctx);
    status |= gr_div(val[2], val[2], c4, ctx);
    status |= gr_one(val[3], ctx);

    GR_TMP_CLEAR5(snu, cnu, dnu, snv, cnv, ctx);
    GR_TMP_CLEAR4(dnv, den, k2, x, ctx);
    return status;
}

/*
    The ratios rz_k = theta_k(z0) / theta_1(z0) from those at 2 z0 (v, in
    any normalization). In the units of _mod_theta_divide, with
    p = wp~(2 z0) and A_1 = c4 r_2, A_2 = c2 c4 r_3, A_3 = c2 r_4 the square
    roots of p - e_1, p - e_2, p - e_3 (r_k the ratios at 2 z0), the halving
    formula X = wp~(z0) = p + A_1 A_2 + A_1 A_3 + A_2 A_3 (with the signs
    of one of the four halves, chosen numerically) is rational in the
    values at 2 z0; then r_2(z0)^2 = (X - e_1)/c4^2, r_4(z0)^2 = (X - e_3)/c2^2
    (two square roots, their signs chosen numerically), and r_3(z0) from the duplication formulas theta_4(2z) theta_4^3
    = theta_4(z)^4 - theta_1(z)^4, theta_1(2z) theta_2 theta_3 theta_4 =
    2 theta_1(z) theta_2(z) theta_3(z) theta_4(z):
    r_3 = (r_4^4 - 1) c2 / (2 c4^2 r_2 r_4 q_4), q_4 = theta_4(2z0) / theta_1(2z0).
*/
static int
_mod_halve_ratios(gr_ptr * rz, gr_ptr * v, gr_srcptr z0, gr_srcptr tau0, gr_srcptr c2, gr_srcptr c4, gr_ctx_t ctx)
{
    gr_ptr e1, e2, e3, p, A1, A2, A3, X, cand, t, u, q4, R2, R4;
    acb_t az, at, w1, w2, w3, w4, z, ref, val;
    int status, e2s, e3s, nfound = 0;

    GR_TMP_INIT5(e1, e2, e3, p, A1, ctx);
    GR_TMP_INIT5(A2, A3, X, cand, t, ctx);
    GR_TMP_INIT4(u, q4, R2, R4, ctx);
    acb_init(az); acb_init(at); acb_init(w1); acb_init(w2); acb_init(w3); acb_init(w4);
    acb_init(z); acb_init(ref); acb_init(val);

    /* numerically: the ratios and wp~ at z0 */
    status = gr_tower_lazy_get_acb(az, z0, 128, ctx);
    status |= gr_tower_lazy_get_acb(at, tau0, 128, ctx);
    if (status == GR_SUCCESS)
    {
        acb_modular_theta(w1, w2, w3, w4, az, at, 128);
        acb_div(w2, w2, w1, 128);
        acb_div(w3, w3, w1, 128);
        acb_div(w4, w4, w1, 128);
        acb_elliptic_p(ref, az, at, 128);
        {
            /* / (pi^2 theta_3^4) */
            acb_t s1, s2, s4;
            acb_init(s1); acb_init(s2); acb_init(s4);
            acb_zero(val);
            acb_modular_theta(s1, s2, z, s4, val, at, 128);
            acb_clear(s1); acb_clear(s2); acb_clear(s4);
        }
        acb_pow_ui(z, z, 4, 128);
        acb_const_pi(val, 128);
        acb_mul(val, val, val, 128);
        acb_mul(z, z, val, 128);
        acb_div(ref, ref, z, 128);
    }

    /* e_k in the units of theta_3 */
    if (status == GR_SUCCESS)
    {
        status = gr_pow_ui(t, c2, 4, ctx);
        status |= gr_pow_ui(u, c4, 4, ctx);
        status |= gr_add_ui(e1, u, 1, ctx);
        status |= gr_div_ui(e1, e1, 3, ctx);
        status |= gr_sub(e2, t, u, ctx);
        status |= gr_div_ui(e2, e2, 3, ctx);
        status |= gr_add_ui(e3, t, 1, ctx);
        status |= gr_div_si(e3, e3, -3, ctx);
    }

    /* A_k, p from the ratios at 2 z0 */
    if (status == GR_SUCCESS)
    {
        status = gr_div(A1, v[1], v[0], ctx);
        status |= gr_mul(A1, A1, c4, ctx);
        status |= gr_div(A2, v[2], v[0], ctx);
        status |= gr_mul(A2, A2, c2, ctx);
        status |= gr_mul(A2, A2, c4, ctx);
        status |= gr_div(q4, v[3], v[0], ctx);
        status |= gr_mul(A3, q4, c2, ctx);
        status |= gr_sqr(p, A1, ctx);
        status |= gr_add(p, p, e1, ctx);
    }

    for (e2s = -1; e2s <= 1 && status == GR_SUCCESS; e2s += 2)
    {
        for (e3s = -1; e3s <= 1 && status == GR_SUCCESS; e3s += 2)
        {
            status = gr_mul(t, A1, A2, ctx);
            status |= gr_mul_si(cand, t, e2s, ctx);
            status |= gr_mul(t, A1, A3, ctx);
            status |= gr_mul_si(t, t, e3s, ctx);
            status |= gr_add(cand, cand, t, ctx);
            status |= gr_mul(t, A2, A3, ctx);
            status |= gr_mul_si(t, t, e2s * e3s, ctx);
            status |= gr_add(cand, cand, t, ctx);
            status |= gr_add(cand, cand, p, ctx);
            if (status == GR_SUCCESS && gr_tower_lazy_get_acb(val, cand, 128, ctx) == GR_SUCCESS &&
                acb_overlaps(val, ref))
            {
                nfound++;
                status = gr_set(X, cand, ctx);
            }
        }
    }
    if (status == GR_SUCCESS && nfound != 1)
        status = GR_UNABLE;

    /* r_2 = sqrt(R2), r_4 = sqrt(R4) */
    if (status == GR_SUCCESS)
    {
        status = gr_one(rz[0], ctx);
        status |= gr_sub(R2, X, e1, ctx);
        status |= gr_sqr(u, c4, ctx);
        status |= gr_div(R2, R2, u, ctx);
        if (status == GR_SUCCESS)
            status = _mod_sqrt_near(rz[1], R2, w2, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = gr_sub(R4, X, e3, ctx);
        status |= gr_sqr(u, c2, ctx);
        status |= gr_div(R4, R4, u, ctx);
        if (status == GR_SUCCESS)
            status = _mod_sqrt_near(rz[3], R4, w4, ctx);
    }

    /* r_3 = (r_4^4 - 1) c2 / (2 c4^2 r_2 r_4 q_4)
           = (R4^2 - 1) c2 r_2 r_4 / (2 c4^2 R2 R4 q_4) (no division by
       the new square roots) */
    if (status == GR_SUCCESS)
    {
        status = gr_sqr(t, R4, ctx);
        status |= gr_sub_ui(t, t, 1, ctx);
        status |= gr_mul(t, t, c2, ctx);
        status |= gr_sqr(u, c4, ctx);
        status |= gr_mul_ui(u, u, 2, ctx);
        status |= gr_mul(u, u, R2, ctx);
        status |= gr_mul(u, u, R4, ctx);
        status |= gr_mul(u, u, q4, ctx);
        status |= gr_div(t, t, u, ctx);
        status |= gr_mul(t, t, rz[1], ctx);
        status |= gr_mul(rz[2], t, rz[3], ctx);
        /* (a check of the choices) */
        if (status == GR_SUCCESS && (gr_tower_lazy_get_acb(val, rz[2], 128, ctx) != GR_SUCCESS || !acb_overlaps(val, w3)))
            status = GR_UNABLE;
    }

    GR_TMP_CLEAR5(e1, e2, e3, p, A1, ctx);
    GR_TMP_CLEAR5(A2, A3, X, cand, t, ctx);
    GR_TMP_CLEAR4(u, q4, R2, R4, ctx);
    acb_clear(az); acb_clear(at); acb_clear(w1); acb_clear(w2); acb_clear(w3); acb_clear(w4);
    acb_clear(z); acb_clear(ref); acb_clear(val);
    return status;
}

/*
    With P = z + w, Q = z - w and S_ab = theta_a(P) theta_b(Q) +
    theta_b(P) theta_a(Q), D_ab = theta_a(P) theta_b(Q) - theta_b(P)
    theta_a(Q), Jacobi's formulas for theta_a(z + w) theta_b(z - w) give

        S_14 = 2 theta_1 theta_4(z) theta_2 theta_3(w) / (c2 c3),
        D_14 = 2 theta_2 theta_3(z) theta_1 theta_4(w) / (c2 c3),
        S_23 = 2 theta_2 theta_3(z) theta_2 theta_3(w) / (c2 c3)

    and their images under the permutations of 2, 3, 4. So with r_k =
    theta_k / theta_1 at z = (P + Q) / 2, theta_1 theta_4(z) / theta_2
    theta_3(z) = r_4 / (r_2 r_3) = S_14 / S_23, and so on:

        r_2^2 = S_23 S_24 / (S_14 S_13),  r_3^2 = S_23 S_34 / (S_14 S_12),
        r_4^2 = S_34 S_24 / (S_12 S_13),  r_3 = (S_13 / S_24) r_2 r_4,

    of degree (2, 2) in the values at P, Q. With the generators h_ab =
    sqrt(S_ab), the ratios are monomials in them, r_2 = h_23 h_24 /
    (h_14 h_13), r_3 = h_23 h_34 / (h_14 h_12), r_4 = h_34 h_24 / (h_12
    h_13) (the signs chosen numerically: the square roots h_ab are only
    determined up to sign, and the products giving r_k are fixed by the
    enclosures of the r_k at z, then checked for consistency), so that the relations between the values at z, w, P, Q are
    rational in the S_ab, D_ab: the square roots of the ratios themselves
    would be roots of large rational functions. The divisions by S_ab are
    D_ab / E_ab, E_ab = S_ab D_ab a combination of squares, which is free
    of the algebraic generators (theta_2, theta_3 at P and Q, say). (The
    halving of the ratios at P + Q, _mod_halve_ratios, is the fallback when
    some S_ab vanishes.)
*/

/* res = x / S = x D / (S D) */
static int
_mod_div_sd(gr_ptr res, gr_srcptr x, gr_srcptr S, gr_srcptr D, gr_ctx_t ctx)
{
    gr_ptr E;
    int status;
    GR_TMP_INIT(E, ctx);
    status = gr_mul(E, S, D, ctx);
    status |= gr_mul(res, x, D, ctx);
    status |= gr_div(res, res, E, ctx);
    GR_TMP_CLEAR(E, ctx);
    return status;
}

static int
_mod_sd_forms(gr_ptr S[4][4], gr_ptr D[4][4], gr_ptr * tp, gr_ptr * tq, gr_ctx_t ctx)
{
    gr_ptr t, u;
    int a, b, status = GR_SUCCESS;
    GR_TMP_INIT2(t, u, ctx);
    for (a = 0; a < 4 && status == GR_SUCCESS; a++)
    {
        for (b = a + 1; b < 4 && status == GR_SUCCESS; b++)
        {
            status = gr_mul(t, tp[a], tq[b], ctx);
            status |= gr_mul(u, tp[b], tq[a], ctx);
            status |= gr_add(S[a][b], t, u, ctx);
            status |= gr_sub(D[a][b], t, u, ctx);
        }
    }
    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

#define SD_INIT(S, D) \
    do { int _a, _b; for (_a = 0; _a < 4; _a++) for (_b = _a + 1; _b < 4; _b++) \
        { GR_TMP_INIT(S[_a][_b], ctx); GR_TMP_INIT(D[_a][_b], ctx); } } while (0)
#define SD_CLEAR(S, D) \
    do { int _a, _b; for (_a = 0; _a < 4; _a++) for (_b = _a + 1; _b < 4; _b++) \
        { GR_TMP_CLEAR(S[_a][_b], ctx); GR_TMP_CLEAR(D[_a][_b], ctx); } } while (0)

static int
_mod_half_sum_ratios(gr_ptr * rz, gr_ptr * tp, gr_ptr * tq, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr S[4][4], D[4][4], h[4][4], t, u;
    acb_t az, at, w1, w[5], val;
    int a, b, k, status;
    /* r_k = h_(n1) h_(n2) / (h_(d1) h_(d2)), as index pairs */
    static const int mono[3][9] = {
        {2, 1, 2, 1, 3, 0, 3, 0, 2},    /* r_2 = h_23 h_24 / (h_14 h_13) */
        {3, 1, 2, 2, 3, 0, 3, 0, 1},    /* r_3 = h_23 h_34 / (h_14 h_12) */
        {4, 2, 3, 1, 3, 0, 1, 0, 2} };  /* r_4 = h_34 h_24 / (h_12 h_13) */

    SD_INIT(S, D);
    for (a = 0; a < 4; a++)
        for (b = a + 1; b < 4; b++)
            GR_TMP_INIT(h[a][b], ctx);
    GR_TMP_INIT2(t, u, ctx);
    acb_init(az); acb_init(at); acb_init(w1); acb_init(val);
    for (k = 0; k < 5; k++)
        acb_init(w[k]);

    status = _mod_sd_forms(S, D, tp, tq, ctx);

    /* the divisors S_12, S_13, S_14 (and S_24 below) */
    if (status == GR_SUCCESS && (gr_is_zero(S[0][1], ctx) != T_FALSE ||
        gr_is_zero(S[0][2], ctx) != T_FALSE || gr_is_zero(S[0][3], ctx) != T_FALSE ||
        gr_is_zero(S[1][3], ctx) != T_FALSE))
        status = GR_UNABLE;

    if (status == GR_SUCCESS)
    {
        status = gr_tower_lazy_get_acb(az, z0, 128, ctx);
        status |= gr_tower_lazy_get_acb(at, tau0, 128, ctx);
        if (status == GR_SUCCESS)
        {
            acb_modular_theta(w1, w[2], w[3], w[4], az, at, 128);
            for (k = 2; k <= 4; k++)
                acb_div(w[k], w[k], w1, 128);
        }
    }

    for (a = 0; a < 4 && status == GR_SUCCESS; a++)
        for (b = a + 1; b < 4 && status == GR_SUCCESS; b++)
            status = _mod_sqrt_any(h[a][b], S[a][b], ctx);

    if (status == GR_SUCCESS)
        status = gr_one(rz[0], ctx);
    for (k = 0; k < 3 && status == GR_SUCCESS; k++)
    {
        const int * m = mono[k];
        status = gr_mul(t, h[m[1]][m[2]], h[m[3]][m[4]], ctx);
        status |= gr_mul(t, t, h[m[5]][m[6]], ctx);
        status |= gr_mul(t, t, h[m[7]][m[8]], ctx);
        status |= _mod_div_sd(t, t, S[m[5]][m[6]], D[m[5]][m[6]], ctx);
        status |= _mod_div_sd(rz[m[0] - 1], t, S[m[7]][m[8]], D[m[7]][m[8]], ctx);
        if (status == GR_SUCCESS)
            status = _mod_fix_sign(rz[m[0] - 1], w[m[0]], ctx);
    }

    /* the consistency of the signs: r_3 S_24 = S_13 r_2 r_4 */
    if (status == GR_SUCCESS)
    {
        status = gr_mul(t, rz[1], rz[3], ctx);
        status |= gr_mul(t, t, S[0][2], ctx);
        status |= gr_mul(u, rz[2], S[1][3], ctx);
        if (status == GR_SUCCESS && gr_equal(t, u, ctx) != T_TRUE)
            status = GR_UNABLE;
    }

    SD_CLEAR(S, D);
    for (a = 0; a < 4; a++)
        for (b = a + 1; b < 4; b++)
            GR_TMP_CLEAR(h[a][b], ctx);
    GR_TMP_CLEAR2(t, u, ctx);
    acb_clear(az); acb_clear(at); acb_clear(w1); acb_clear(val);
    for (k = 0; k < 5; k++)
        acb_clear(w[k]);
    return status;
}

/*
    The ratios rw_k = theta_k(w) / theta_1(w) from those at z (values ty)
    with the values tp at P = z + w and tq at Q = z - w: by the formulas
    above (with D_12 = 2 theta_3 theta_4(z) theta_1 theta_2(w) / (c3 c4),
    and so on),

        r_2(w) / r_2(z) = (c4 / c3) S_13 / D_14,
        r_3(w) / r_3(z) = (c2 / c4) S_14 / D_12,
        r_4(w) / r_4(z) = (c3 / c2) S_12 / D_13,

    rational: the values at z, z + w and z - w determine the ratios at w
    (the halving of 2 w = P - Q from P, Q alone leaves four choices).
    c2, c4 normalized by theta_3.
*/
static int
_mod_half_diff_ratios(gr_ptr * rw, gr_ptr * ty, gr_ptr * tp, gr_ptr * tq, gr_srcptr c2, gr_srcptr c4, gr_ctx_t ctx)
{
    gr_ptr S[4][4], D[4][4], t;
    int k, status = GR_SUCCESS;
    /* (k, a, b of S_ab, c, d of D_cd) */
    static const int tab[3][5] = { {2, 1, 3, 1, 4}, {3, 1, 4, 1, 2}, {4, 1, 2, 1, 3} };

    SD_INIT(S, D);
    GR_TMP_INIT(t, ctx);

    if (gr_is_zero(ty[0], ctx) != T_FALSE)
        status = GR_UNABLE;
    if (status == GR_SUCCESS)
        status = _mod_sd_forms(S, D, tp, tq, ctx);
    if (status == GR_SUCCESS && (gr_is_zero(D[0][1], ctx) != T_FALSE ||
        gr_is_zero(D[0][2], ctx) != T_FALSE || gr_is_zero(D[0][3], ctx) != T_FALSE))
        status = GR_UNABLE;

    for (k = 0; k < 3 && status == GR_SUCCESS; k++)
    {
        int kk = tab[k][0], a = tab[k][1] - 1, b = tab[k][2] - 1, c = tab[k][3] - 1, d = tab[k][4] - 1;

        /* S_ab / D_cd = S_ab S_cd / (D_cd S_cd) */
        status = _mod_div_sd(t, S[a][b], D[c][d], S[c][d], ctx);
        if (kk == 2)
            status |= gr_mul(t, t, c4, ctx);
        else if (kk == 3)
        {
            status |= gr_mul(t, t, c2, ctx);
            status |= gr_div(t, t, c4, ctx);
        }
        else
            status |= gr_div(t, t, c2, ctx);
        status |= gr_mul(t, t, ty[kk - 1], ctx);
        status |= gr_div(rw[kk - 1], t, ty[0], ctx);
    }
    if (status == GR_SUCCESS)
        status = gr_one(rw[0], ctx);

    SD_CLEAR(S, D);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/*
    TP_ADD_HALF with y of type TP_HALF, y = (P + Q) / 2 with the values at
    P, Q and the generators h_ab = sqrt(S_ab) (_mod_half_sum_ratios), and
    z0 = +-(P - Q) / 2 modulo the half lattice: the values at w = (P - Q) / 2
    with r_k(w) = kappa_k r_k(y) (_mod_half_diff_ratios) and
    theta_1(y) theta_1(w) = sqrt(V),

        V = theta_4^2 S_12 S_13 D_12 D_13 / (2 S_23 S_24 S_34)

    (by S_12 S_13 - D_12 D_13 = 2 theta_1(P) theta_1(Q) S_23 in the
    addition formula for theta_1), a monomial in the h_ab and one new
    generator sqrt(D_12 D_13): the values at y and w then satisfy the
    addition formulas as rational identities in the S_ab, D_ab. The
    values at z0 by the reduction of w.
*/
static int
_mod_tp_add_half_sym(gr_ptr * th, slong iy, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_theta_point_struct * ey = LAZY(ctx)->theta_pts + iy;
    gr_ptr y, v, P, Q, w, t, m, c2, c4, t4, S[4][4], D[4][4], ty[4], tp[4], tq[4], rz[4];
    acb_t ay, aw, at, r1, r2, r3, r4, ref, val;
    slong j;
    int s, s2, k, status = GR_SUCCESS;

    if (ey->kind != TP_HALF)
        return GR_UNABLE;
    j = ey->j; s = ey->s; s2 = ey->s2;

    GR_TMP_INIT5(y, v, P, Q, w, ctx);
    GR_TMP_INIT5(t, m, c2, c4, t4, ctx);
    GR_TMP_INIT4(ty[0], ty[1], ty[2], ty[3], ctx);
    GR_TMP_INIT4(tp[0], tp[1], tp[2], tp[3], ctx);
    GR_TMP_INIT4(tq[0], tq[1], tq[2], tq[3], ctx);
    GR_TMP_INIT4(rz[0], rz[1], rz[2], rz[3], ctx);
    SD_INIT(S, D);
    acb_init(ay); acb_init(aw); acb_init(at); acb_init(ref); acb_init(val);
    acb_init(r1); acb_init(r2); acb_init(r3); acb_init(r4);

    /* P = 2y - v, Q = v (as in the evaluation of y), w = (P - Q) / 2 */
    status = _mod_tp_z(y, iy, ctx);
    status |= _mod_tp_z(v, j, ctx);
    status |= gr_mul_si(v, v, s * s2, ctx);
    status |= gr_mul_ui(P, y, 2, ctx);
    status |= gr_sub(P, P, v, ctx);
    status |= gr_set(Q, v, ctx);
    status |= gr_sub(w, P, Q, ctx);
    status |= gr_div_ui(w, w, 2, ctx);
    if (status == GR_SUCCESS && !_mod_half_lattice(z0, w, tau0, ctx))
    {
        status = gr_neg(t, w, ctx);
        if (status == GR_SUCCESS && !_mod_half_lattice(z0, t, tau0, ctx))
            status = GR_UNABLE;
    }

    if (status == GR_SUCCESS)
        status = _mod_theta_z_reduced(ty, y, tau0, ctx);
    if (status == GR_SUCCESS)
        status = _mod_theta_z_reduced(tp, P, tau0, ctx);
    if (status == GR_SUCCESS)
        status = _mod_theta_z_reduced(tq, Q, tau0, ctx);
    if (status == GR_SUCCESS)
        status = _mod_point(c2, NULL, c4, NULL, NULL, tau0, 1, ctx);
    if (status == GR_SUCCESS)
        status = _mod_theta0(NULL, NULL, t4, NULL, tau0, ctx);
    if (status == GR_SUCCESS)
        status = _mod_half_diff_ratios(rz, ty, tp, tq, c2, c4, ctx);
    if (status == GR_SUCCESS)
        status = _mod_sd_forms(S, D, tp, tq, ctx);
    if (status == GR_SUCCESS && (gr_is_zero(S[1][2], ctx) != T_FALSE ||
        gr_is_zero(S[1][3], ctx) != T_FALSE || gr_is_zero(S[2][3], ctx) != T_FALSE))
        status = GR_UNABLE;

    /* m = theta_4 h_12 h_13 sqrt(D_12 D_13) / (sqrt(2) h_23 h_24 h_34)
         = theta_4 h_12 h_13 sqrt(D_12 D_13) h_23 h_24 h_34 sqrt(2) / (2 S_23 S_24 S_34) */
    if (status == GR_SUCCESS)
    {
        status = gr_mul(t, D[0][1], D[0][2], ctx);
        if (status == GR_SUCCESS)
            status = _mod_sqrt_any(m, t, ctx);
        status |= gr_mul(m, m, t4, ctx);
        for (k = 0; k < 5 && status == GR_SUCCESS; k++)
        {
            static const int ab[5][2] = { {0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3} };
            status = _mod_sqrt_any(t, S[ab[k][0]][ab[k][1]], ctx);
            status |= gr_mul(m, m, t, ctx);
            if (k >= 2)
                status |= _mod_div_sd(m, m, S[ab[k][0]][ab[k][1]], D[ab[k][0]][ab[k][1]], ctx);
        }
        status |= gr_set_ui(t, 2, ctx);
        status |= gr_sqrt(t, t, ctx);
        status |= gr_mul(m, m, t, ctx);
        status |= gr_div_ui(m, m, 2, ctx);
    }

    /* the sign: m = theta_1(y) theta_1(w) */
    if (status == GR_SUCCESS)
    {
        status = gr_tower_lazy_get_acb(ay, y, 128, ctx);
        status |= gr_tower_lazy_get_acb(aw, w, 128, ctx);
        status |= gr_tower_lazy_get_acb(at, tau0, 128, ctx);
        if (status == GR_SUCCESS)
        {
            acb_modular_theta(r1, r2, r3, r4, ay, at, 128);
            acb_set(ref, r1);
            acb_modular_theta(r1, r2, r3, r4, aw, at, 128);
            acb_mul(ref, ref, r1, 128);
            status = _mod_fix_sign(m, ref, ctx);
        }
    }

    /* the values at w, then at its reduced point z0 */
    if (status == GR_SUCCESS)
    {
        status = gr_div(m, m, ty[0], ctx);
        for (k = 0; k < 4; k++)
            status |= gr_mul(tp[k], rz[k], m, ctx);
        if (status == GR_SUCCESS)
            status = _mod_tp_unmap(th, tp, w, tau0, ctx);
    }

    GR_TMP_CLEAR5(y, v, P, Q, w, ctx);
    GR_TMP_CLEAR5(t, m, c2, c4, t4, ctx);
    GR_TMP_CLEAR4(ty[0], ty[1], ty[2], ty[3], ctx);
    GR_TMP_CLEAR4(tp[0], tp[1], tp[2], tp[3], ctx);
    GR_TMP_CLEAR4(tq[0], tq[1], tq[2], tq[3], ctx);
    GR_TMP_CLEAR4(rz[0], rz[1], rz[2], rz[3], ctx);
    SD_CLEAR(S, D);
    acb_clear(ay); acb_clear(aw); acb_clear(at); acb_clear(ref); acb_clear(val);
    acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(r4);
    return status;
}

/* the values at z0 by the method e */
static int
_mod_tp_eval(gr_ptr * th, const gr_tower_lazy_theta_point_struct * e, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr u, v, P, Q, t, tu[4], tv[4], tq[4], val[4];
    int k, status = GR_SUCCESS;

    if (e->kind == TP_STANDARD)
        return _mod_theta_z0_gens(th, z0, tau0, ctx);
    if (e->kind == TP_TORSION)
        return _mod_theta_divide(th, e->n, z0, NULL, tau0, ctx);

    GR_TMP_INIT5(u, v, P, Q, t, ctx);
    GR_TMP_INIT4(tu[0], tu[1], tu[2], tu[3], ctx);
    GR_TMP_INIT4(tv[0], tv[1], tv[2], tv[3], ctx);
    GR_TMP_INIT4(tq[0], tq[1], tq[2], tq[3], ctx);
    GR_TMP_INIT4(val[0], val[1], val[2], val[3], ctx);

    if (e->kind == TP_MULT)
    {
        /* the values at P = s n u, then (dv > 1) at Q = dv z0 (the same
           reduced point) and at z0 by the division */
        status = _mod_tp_z(u, e->i, ctx);
        status |= gr_mul_si(u, u, e->s, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tu, u, tau0, ctx);
        if (status == GR_SUCCESS && e->n > 1)
            status = _mod_theta_mult_vals(val, tu, e->n, tau0, ctx);
        else
            for (k = 0; k < 4; k++)
                gr_swap(val[k], tu[k], ctx);
        status |= gr_mul_si(P, u, e->n, ctx);
        if (status == GR_SUCCESS && e->dv == 1)
            status = _mod_tp_unmap(th, val, P, tau0, ctx);
        else if (status == GR_SUCCESS)
        {
            status = _mod_tp_unmap(tu, val, P, tau0, ctx);
            status |= gr_mul_si(Q, z0, e->dv, ctx);
            if (status == GR_SUCCESS)
                status = _mod_tp_map(tq, tu, Q, tau0, ctx);
            if (status == GR_SUCCESS)
                status = _mod_theta_divide(th, e->dv, z0, tq, tau0, ctx);
        }
    }
    else if (e->kind == TP_ADD)
    {
        /* P = s (u + v), Q = u - v, v = s2 u_j:
           theta_k(u + v) = A_k / (theta_k(u - v) theta_4^2) */
        gr_ptr t4;
        GR_TMP_INIT(t4, ctx);
        status = _mod_tp_z(u, e->i, ctx);
        status |= _mod_tp_z(v, e->j, ctx);
        status |= gr_mul_si(v, v, e->s2, ctx);
        status |= gr_add(P, u, v, ctx);
        status |= gr_sub(Q, u, v, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tu, u, tau0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tv, v, tau0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tq, Q, tau0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta0(NULL, NULL, t4, NULL, tau0, ctx);
        status |= gr_sqr(t4, t4, ctx);
        for (k = 1; k <= 4 && status == GR_SUCCESS; k++)
        {
            /* A_k = theta_k(u)^2 theta_4(v)^2 - theta_{5-k}(u)^2 theta_1(v)^2 */
            status = gr_mul(t, tu[k - 1], tv[3], ctx);
            status |= gr_sqr(t, t, ctx);
            status |= gr_mul(val[k - 1], tu[4 - k], tv[0], ctx);
            status |= gr_sqr(val[k - 1], val[k - 1], ctx);
            status |= gr_sub(val[k - 1], t, val[k - 1], ctx);
            status |= gr_mul(t, tq[k - 1], t4, ctx);
            if (status == GR_SUCCESS && gr_is_zero(t, ctx) != T_FALSE)
                status = GR_UNABLE;
            if (status == GR_SUCCESS)
                status = gr_div(val[k - 1], val[k - 1], t, ctx);
        }
        /* theta at s P: theta_1 odd */
        if (status == GR_SUCCESS && e->s < 0)
        {
            status = gr_neg(val[0], val[0], ctx);
            status |= gr_neg(P, P, ctx);
        }
        if (status == GR_SUCCESS)
            status = _mod_tp_unmap(th, val, P, tau0, ctx);
        GR_TMP_CLEAR(t4, ctx);
    }
    else if (e->kind == TP_ADD_HALF && _mod_tp_add_half_sym(th, e->i, z0, tau0, ctx) == GR_SUCCESS)
    {
        status = GR_SUCCESS;
    }
    else if (e->kind == TP_ADD_HALF)
    {
        /* y = u_i, P = y + z0, Q = y - z0, a = theta_1(y)^2, b = theta_4(y)^2,
           alpha_k = theta_k(P) theta_k(Q) theta_4^2: alpha_1 = a Y - b X,
           alpha_4 = b Y - a X for X = theta_1(z0)^2, Y = theta_4(z0)^2;
           then theta_2(z0)^2 c4^2 = Y c2^2 - X, theta_3(z0)^2 c4^2 = Y - X c2^2
           (four square roots, the addition formulas then rational identities
           in the values at y, P, Q). If a^2 = b^2: the ratios at z0
           rational in the values at y, P, Q (_mod_half_diff_ratios; otherwise
           by the halving of P + (-Q)), and theta_1(z0) theta_1(y) =
           sqrt(alpha_1 / (r_4(z0)^2 - r_4(y)^2)) */
        gr_ptr c2, c4, t42, a, b, al1, al4, X, Y, rz[4];
        acb_t az, ay, at, r1, r2, r3, r4, ref;
        GR_TMP_INIT5(c2, c4, t42, a, b, ctx);
        GR_TMP_INIT4(al1, al4, X, Y, ctx);
        GR_TMP_INIT4(rz[0], rz[1], rz[2], rz[3], ctx);
        acb_init(az); acb_init(ay); acb_init(at); acb_init(ref);
        acb_init(r1); acb_init(r2); acb_init(r3); acb_init(r4);

        status = _mod_tp_z(u, e->i, ctx);
        status |= gr_add(P, u, z0, ctx);
        status |= gr_sub(Q, u, z0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tu, u, tau0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tv, P, tau0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tq, Q, tau0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_point(c2, NULL, c4, NULL, NULL, tau0, 1, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta0(NULL, NULL, t42, NULL, tau0, ctx);
        if (status == GR_SUCCESS)
        {
            status = gr_sqr(t42, t42, ctx);
            status |= gr_sqr(a, tu[0], ctx);
            status |= gr_sqr(b, tu[3], ctx);
            status |= gr_mul(al1, tv[0], tq[0], ctx);
            status |= gr_mul(al1, al1, t42, ctx);
            status |= gr_mul(al4, tv[3], tq[3], ctx);
            status |= gr_mul(al4, al4, t42, ctx);
            status |= gr_sqr(v, a, ctx);
            status |= gr_sqr(t, b, ctx);
            status |= gr_sub(v, v, t, ctx);
        }
        if (status == GR_SUCCESS)
        {
            status = gr_tower_lazy_get_acb(az, z0, 128, ctx);
            status |= gr_tower_lazy_get_acb(ay, u, 128, ctx);
            status |= gr_tower_lazy_get_acb(at, tau0, 128, ctx);
            if (status == GR_SUCCESS)
            {
                acb_modular_theta(r1, r2, r3, r4, ay, at, 128);
                acb_set(ref, r1);
                acb_modular_theta(r1, r2, r3, r4, az, at, 128);
                acb_mul(ref, ref, r1, 128);
            }
        }

        if (status == GR_SUCCESS && gr_is_zero(v, ctx) == T_FALSE)
        {
            /* X = (b alpha_1 - a alpha_4) / (a^2 - b^2), Y = (a alpha_1 - b alpha_4) / (a^2 - b^2) */
            status = gr_mul(X, b, al1, ctx);
            status |= gr_submul(X, a, al4, ctx);
            status |= gr_div(X, X, v, ctx);
            status |= gr_mul(Y, a, al1, ctx);
            status |= gr_submul(Y, b, al4, ctx);
            status |= gr_div(Y, Y, v, ctx);
            if (status == GR_SUCCESS)
                status = _mod_sqrt_near(th[0], X, r1, ctx);
            if (status == GR_SUCCESS)
                status = _mod_sqrt_near(th[3], Y, r4, ctx);
            if (status == GR_SUCCESS)
            {
                status = gr_sqr(v, c2, ctx);
                status |= gr_sqr(Q, c4, ctx);
                status |= gr_mul(t, Y, v, ctx);
                status |= gr_sub(t, t, X, ctx);
                status |= gr_div(t, t, Q, ctx);
                if (status == GR_SUCCESS)
                    status = _mod_sqrt_near(th[1], t, r2, ctx);
                status |= gr_mul(t, X, v, ctx);
                status |= gr_sub(t, Y, t, ctx);
                status |= gr_div(t, t, Q, ctx);
                if (status == GR_SUCCESS)
                    status = _mod_sqrt_near(th[2], t, r3, ctx);
            }
        }
        else if (status == GR_SUCCESS)
        {
            if (_mod_half_diff_ratios(rz, tu, tv, tq, c2, c4, ctx) != GR_SUCCESS)
            {
                for (k = 0; k < 4; k++)
                    status |= gr_set(val[k], tq[k], ctx);
                status |= gr_neg(val[0], val[0], ctx);
                if (status == GR_SUCCESS)
                    status = _mod_jacobi_add(val, tv, val, c2, c4, ctx);
                if (status == GR_SUCCESS)
                    status = _mod_halve_ratios(rz, val, z0, tau0, c2, c4, ctx);
            }
            /* X = alpha_1 / (r_4(z0)^2 - r_4(y)^2) = theta_1(y)^2 theta_1(z0)^2 */
            if (status == GR_SUCCESS)
            {
                status = gr_div(t, tu[3], tu[0], ctx);
                status |= gr_sqr(t, t, ctx);
                status |= gr_sqr(v, rz[3], ctx);
                status |= gr_sub(v, v, t, ctx);
                if (status == GR_SUCCESS && gr_is_zero(v, ctx) != T_FALSE)
                    status = GR_UNABLE;
                if (status == GR_SUCCESS)
                    status = gr_div(X, al1, v, ctx);
            }
            if (status == GR_SUCCESS)
                status = _mod_sqrt_near(t, X, ref, ctx);
            status |= gr_div(t, t, tu[0], ctx);
            for (k = 0; k < 4 && status == GR_SUCCESS; k++)
                status = gr_mul(th[k], rz[k], t, ctx);
        }

        GR_TMP_CLEAR5(c2, c4, t42, a, b, ctx);
        GR_TMP_CLEAR4(al1, al4, X, Y, ctx);
        GR_TMP_CLEAR4(rz[0], rz[1], rz[2], rz[3], ctx);
        acb_clear(az); acb_clear(ay); acb_clear(at); acb_clear(ref);
        acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(r4);
    }
    else if (e->kind == TP_JACOBI || e->kind == TP_HALF)
    {
        /* TP_JACOBI: T = s (n u + s2 u_j), the ratios at T (the map from
           T to z0 is linear); TP_HALF: T = s (u + s2 u_j) = 2 z0 modulo
           the half lattice, the ratios at T mapped to 2 z0 and halved;
           theta_1(z0) the generator */
        gr_ptr c2, c4, rho[4], g1, args;
        slong sz = ctx->sizeof_elem;
        GR_TMP_INIT3(c2, c4, g1, ctx);
        GR_TMP_INIT4(rho[0], rho[1], rho[2], rho[3], ctx);
        args = gr_heap_init_vec(2, ctx);

        status = _mod_tp_z(u, e->i, ctx);
        status |= _mod_tp_z(v, e->j, ctx);
        status |= gr_mul_si(u, u, e->s, ctx);
        status |= gr_mul_si(v, v, e->s * e->s2, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tu, u, tau0, ctx);
        /* (the multiple n u) */
        if (status == GR_SUCCESS && e->kind == TP_JACOBI && e->n > 1)
        {
            status = _mod_theta_mult_vals(tq, tu, e->n, tau0, ctx);
            for (k = 0; k < 4; k++)
                gr_swap(tu[k], tq[k], ctx);
            status |= gr_mul_si(u, u, e->n, ctx);
        }
        status |= gr_add(P, u, v, ctx);
        if (status == GR_SUCCESS)
            status = _mod_theta_z_reduced(tv, v, tau0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_point(c2, NULL, c4, NULL, NULL, tau0, 1, ctx);

        if (status == GR_SUCCESS && e->kind == TP_JACOBI)
        {
            status = _mod_jacobi_add(val, tu, tv, c2, c4, ctx);
            if (status == GR_SUCCESS)
                status = _mod_tp_unmap(rho, val, P, tau0, ctx);
        }
        else if (status == GR_SUCCESS)
        {
            /* 2 z0 = (u + h) + v, h = 2 z0 - T in the half lattice: the
               ratios at z0 from the values at u + h and v, else by the
               halving */
            status = gr_mul_ui(Q, z0, 2, ctx);
            status |= gr_sub(t, Q, P, ctx);
            status |= gr_add(t, t, u, ctx);
            if (status == GR_SUCCESS)
                status = _mod_theta_z_reduced(tq, t, tau0, ctx);
            if (status == GR_SUCCESS && _mod_half_sum_ratios(rho, tq, tv, z0, tau0, ctx) != GR_SUCCESS)
            {
                status = _mod_jacobi_add(val, tu, tv, c2, c4, ctx);
                if (status == GR_SUCCESS)
                    status = _mod_tp_unmap(tq, val, P, tau0, ctx);
                if (status == GR_SUCCESS)
                    status = _mod_tp_map(val, tq, Q, tau0, ctx);
                if (status == GR_SUCCESS)
                    status = _mod_halve_ratios(rho, val, z0, tau0, c2, c4, ctx);
            }
        }

        if (status == GR_SUCCESS)
        {
            status = gr_set(args, z0, ctx);
            status |= gr_set(GR_ENTRY(args, 1, sz), tau0, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_multi_locked(g1, args, 2, GR_TOWER_JACOBI_THETA, 1, ctx);
        }
        for (k = 0; k < 4 && status == GR_SUCCESS; k++)
        {
            status = gr_div(th[k], rho[k], rho[0], ctx);
            status |= gr_mul(th[k], th[k], g1, ctx);
        }
        GR_TMP_CLEAR3(c2, c4, g1, ctx);
        GR_TMP_CLEAR4(rho[0], rho[1], rho[2], rho[3], ctx);
        gr_heap_clear_vec(args, 2, ctx);
    }
    else
        status = GR_UNABLE;

    GR_TMP_CLEAR5(u, v, P, Q, t, ctx);
    GR_TMP_CLEAR4(tu[0], tu[1], tu[2], tu[3], ctx);
    GR_TMP_CLEAR4(tv[0], tv[1], tv[2], tv[3], ctx);
    GR_TMP_CLEAR4(tq[0], tq[1], tq[2], tq[3], ctx);
    GR_TMP_CLEAR4(val[0], val[1], val[2], val[3], ctx);
    return status;
}

/*
    Whether the division polynomials are solved at tau0: not at the CM
    points where lambda is algebraic of degree more than 2 (with large
    coefficients, from a class polynomial: lambda(1/2 + i) a root of a
    quartic, over which the roots of psi_3 take 20 seconds); the torsion
    points are then generators.
*/
#define TORSION_LAMBDA_DEGREE 2

static int
_mod_torsion_field_ok(gr_srcptr tau0, gr_ctx_t ctx)
{
    fmpz_t D, A, B, C;
    int r, ok;
    fmpz_init(D);
    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(C);
    r = _gr_tower_lazy_modular_cm_form(D, A, B, C, tau0, ctx);
    ok = (r == 0) || (r == 1 && fmpz_cmp_si(D, -4) >= 0);
    if (r == 1 && !ok)
    {
        /* lambda of small degree (through an anchor: lambda(4i) from
           lambda(i) = 1/2, say) */
        gr_ptr L;
        GR_TMP_INIT(L, ctx);
        if (_mod_lambda0(L, tau0, ctx) == GR_SUCCESS)
        {
            slong level;
            gr_tower_struct * T = gr_tower_lazy_get_tower(&level, L, ctx);
            ok = (T == NULL) || (gr_tower_degree_at(T, level) <= TORSION_LAMBDA_DEGREE);
        }
        GR_TMP_CLEAR(L, ctx);
    }
    fmpz_clear(D);
    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(C);
    return ok;
}

/*
    Rebasing a division: z0 with dv z0 = s n u modulo the half lattice
    (dv > 1, gcd(n, dv) = 1), u a recorded point of type TP_STANDARD (its
    values the generators theta_k(u)), would be evaluated by the division
    by dv (roots of polynomials of degree dv^2, of large coefficients).
    Instead the point b = a z0 + c s u (a n + c dv = 1, so that n b = z0
    and dv b = s u modulo the half lattice) gets generators of its own,
    the record of u becomes the multiplication by dv from b and z0 that
    by n from b (type TP_MULT: values rational in those at b), and the
    context records the values of the generators theta_k(u) at b
    (_gr_tower_lazy_rebase_add): the towers containing both make them
    rational functions of the generators at b. So with u evaluated first,
    the later values at u, z0 and b are those of the order b first.
    Sets th to the values at z0 and returns the index of its record, or
    -1 (no change made; th undefined).
*/
static slong
_mod_tp_rebase(gr_ptr * th, const gr_tower_lazy_theta_point_struct * e, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_lazy_theta_point_struct es, eu, ez;
    gr_ptr u, b, b1, t, nv[4];
    ulong def[4], trigger;
    slong ku = e->i, kb = -1, kz = -1, a, c, p, q;
    int k, neg, status;

    if (e->kind != TP_MULT || e->dv <= 1 ||
        L->theta_pts[ku].kind != TP_STANDARD || !L->theta_pts[ku].have_vals)
        return -1;

    /* the generators at u */
    for (k = 0; k < 4; k++)
    {
        def[k] = _gr_tower_lazy_gen_def_of(&L->theta_pts[ku].vals[k], ctx);
        if (def[k] == 0)
            return -1;
    }

    for (a = 0; a < e->dv; a++)
        if ((a * e->n) % e->dv == 1 % e->dv)
            break;
    c = (1 - a * e->n) / e->dv;

    GR_TMP_INIT4(u, b, b1, t, ctx);
    GR_TMP_INIT4(nv[0], nv[1], nv[2], nv[3], ctx);

    /* the reduced point of b, not recorded yet */
    status = _mod_tp_z(u, ku, ctx);
    status |= gr_mul_si(b, z0, a, ctx);
    status |= gr_mul_si(t, u, c * e->s, ctx);
    status |= gr_add(b, b, t, ctx);
    if (status == GR_SUCCESS)
        status = _mod_z_reduce(b, b1, &p, &q, &neg, b, tau0, ctx);
    if (status == GR_SUCCESS && (gr_is_zero(b, ctx) != T_FALSE || _mod_tp_find(b, tau0, ctx) >= 0))
        status = GR_UNABLE;

    /* (the reduced points of dv b and n b are u and z0) */
    if (status == GR_SUCCESS)
    {
        status = gr_mul_si(t, b, e->dv, ctx);
        status |= _mod_z_reduce(t, b1, &p, &q, &neg, t, tau0, ctx);
        if (status == GR_SUCCESS && gr_equal(t, u, ctx) != T_TRUE)
            status = GR_UNABLE;
    }
    if (status == GR_SUCCESS)
    {
        status = gr_mul_si(t, b, e->n, ctx);
        status |= _mod_z_reduce(t, b1, &p, &q, &neg, t, tau0, ctx);
        if (status == GR_SUCCESS && gr_equal(t, z0, ctx) != T_TRUE)
            status = GR_UNABLE;
    }

    /* the generators at b, recorded */
    es.kind = TP_STANDARD;
    es.i = es.j = es.l = -1;
    es.n = es.dv = 0;
    es.s = es.s2 = 1;
    es.have_vals = 0;
    es.stamp = 0;
    if (status == GR_SUCCESS)
        status = _mod_tp_eval(th, &es, b, tau0, ctx);
    trigger = (status == GR_SUCCESS) ? _gr_tower_lazy_gen_def_of(th[0], ctx) : 0;
    if (trigger != 0)
    {
        kb = _mod_tp_record(b, tau0, &es, ctx);
        _mod_tp_store(kb, th, ctx);

        /* the values at u from b */
        eu = es;
        eu.kind = TP_MULT;
        eu.i = kb;
        eu.n = e->dv;
        eu.dv = 1;
        status = _mod_tp_eval(nv, &eu, u, tau0, ctx);
    }

    if (trigger != 0 && status == GR_SUCCESS)
    {
        /* (cached values of u keep their towers: copies) */
        L->theta_pts[ku].kind = TP_MULT;
        L->theta_pts[ku].i = kb;
        L->theta_pts[ku].j = L->theta_pts[ku].l = -1;
        L->theta_pts[ku].n = e->dv;
        L->theta_pts[ku].dv = 1;
        L->theta_pts[ku].s = 1;
        if (L->theta_pts[ku].have_vals)
            for (k = 0; k < 4; k++)
                GR_MUST_SUCCEED(gr_set(&L->theta_pts[ku].vals[k], nv[k], ctx));

        for (k = 0; k < 4; k++)
            _gr_tower_lazy_rebase_add(def[k], trigger, 0, nv[k], ctx);

        /* z0 = n b */
        if (e->n == 1)
            kz = kb;
        else
        {
            ez = eu;
            ez.n = e->n;
            for (k = 0; k < 4; k++)
                GR_MUST_SUCCEED(gr_set(nv[k], th[k], ctx));
            status = _mod_tp_eval(th, &ez, z0, tau0, ctx);
            if (status != GR_SUCCESS)
            {
                /* (not expected: the division then) */
                ez = *e;
                status = _mod_tp_eval(th, &ez, z0, tau0, ctx);
            }
            if (status == GR_SUCCESS)
            {
                kz = _mod_tp_record(z0, tau0, &ez, ctx);
                _mod_tp_store(kz, th, ctx);
            }
        }
    }
    else if (trigger != 0)
    {
        /* (not expected) b is then the division of u, as without the
           rebase (its generators unused), and so is z0 */
        gr_tower_lazy_theta_point_struct * f = L->theta_pts + kb;
        f->kind = TP_MULT;
        f->i = ku;
        f->n = 1;
        f->dv = e->dv;
        f->s = 1;
        for (k = 0; k < 4; k++)
            _gr_tower_lazy_clear(&f->vals[k], ctx);
        f->have_vals = 0;
        kz = -1;
    }

    GR_TMP_CLEAR4(u, b, b1, t, ctx);
    GR_TMP_CLEAR4(nv[0], nv[1], nv[2], nv[3], ctx);
    return kz;
}

/* whether the method of the record k refers to the record t (through
   the recorded points it uses; conservatively 1 when unsure) */
static int
_mod_tp_refers(slong k, slong t, int depth, gr_ctx_t ctx)
{
    const gr_tower_lazy_theta_point_struct * e;

    if (k == t || depth > 8)
        return 1;
    e = LAZY(ctx)->theta_pts + k;
    switch (e->kind)
    {
        case TP_STANDARD:
        case TP_TORSION:
            return 0;
        case TP_MULT:
            return _mod_tp_refers(e->i, t, depth + 1, ctx);
        case TP_ADD:
            return _mod_tp_refers(e->i, t, depth + 1, ctx) ||
                   _mod_tp_refers(e->j, t, depth + 1, ctx) ||
                   (e->l >= 0 && _mod_tp_refers(e->l, t, depth + 1, ctx));
        case TP_JACOBI:
        case TP_HALF:
            return _mod_tp_refers(e->i, t, depth + 1, ctx) ||
                   _mod_tp_refers(e->j, t, depth + 1, ctx);
        default:
            return 1;
    }
}

/* the sign sigma with u + sigma x reducing to the point p, or 0 */
static int
_mod_tp_sum_sign_pt(gr_srcptr u, gr_srcptr x, gr_srcptr p, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr t, t1;
    slong pp, qq;
    int sigma, neg, res = 0;

    GR_TMP_INIT2(t, t1, ctx);
    for (sigma = 1; sigma >= -1 && res == 0; sigma -= 2)
    {
        if (gr_mul_si(t, x, sigma, ctx) == GR_SUCCESS &&
            gr_add(t, u, t, ctx) == GR_SUCCESS &&
            _mod_z_reduce(t, t1, &pp, &qq, &neg, t, tau0, ctx) == GR_SUCCESS &&
            gr_equal(t, p, ctx) == T_TRUE)
            res = sigma;
    }
    GR_TMP_CLEAR2(t, t1, ctx);
    return res;
}

/* the same for the point of the record k */
static int
_mod_tp_sum_sign(gr_srcptr u, gr_srcptr x, slong k, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr p;
    int res = 0;
    GR_TMP_INIT(p, ctx);
    if (_mod_tp_z(p, k, ctx) == GR_SUCCESS)
        res = _mod_tp_sum_sign_pt(u, x, p, tau0, ctx);
    GR_TMP_CLEAR(p, ctx);
    return res;
}

/* removes the record k if it is the last one (returns 1); otherwise
   keeps it (returns 0: a record of type TP_STANDARD with its values) */
static int
_mod_tp_pop(slong k, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_lazy_theta_point_struct * e = L->theta_pts + k;
    slong q;

    /* (points recorded after it, by nested evaluations, may refer to it:
       kept) */
    if (k != L->num_theta_pts - 1)
        return 0;
    if (e->have_vals)
        for (q = 0; q < 4; q++)
            _gr_tower_lazy_clear(&e->vals[q], ctx);
    _gr_tower_lazy_clear(&e->z0, ctx);
    _gr_tower_lazy_clear(&e->tau0, ctx);
    L->num_theta_pts--;
    return 1;
}

/*
    Whether the addition of the point x to the point y can be rebased
    (see _mod_tp_rebase_add), with the records kP, kQ of the reduced
    points of y + x and y - x (in some order; kQ = -1: the point zq, not
    recorded): sets *kJ (the one to become of type TP_JACOBI), *kA (of
    type TP_ADD; -1 for zq) and the sign sJ with y + sJ x reducing to
    the point of kJ. The record ky of y (-1: not recorded) must not refer
    to them.
*/
static int
_mod_tp_add_check(slong * kJ, slong * kA, int * sJ, gr_srcptr y, slong ky, gr_srcptr x, slong kP, slong kQ, gr_srcptr zq, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    int sA;

    if (kP < 0 || kP == kQ || kP == ky || (kQ >= 0 && kQ == ky) || (kQ < 0 && zq == NULL))
        return 0;
    /* (of type TP_STANDARD: a point of type TP_JACOBI, Q = 2y - P say,
       has values rational in those at y and P, which would become large
       rational functions of the values at y and x, in the elements
       already built on them) */
    if (!L->theta_pts[kP].have_vals || L->theta_pts[kP].kind != TP_STANDARD)
        return 0;
    if (kQ >= 0 && (!L->theta_pts[kQ].have_vals || L->theta_pts[kQ].kind != TP_STANDARD))
        return 0;
    if (ky >= 0 && (_mod_tp_refers(ky, kP, 0, ctx) || (kQ >= 0 && _mod_tp_refers(ky, kQ, 0, ctx))))
        return 0;

    *kJ = kP;
    *kA = kQ;
    *sJ = _mod_tp_sum_sign(y, x, kP, tau0, ctx);
    sA = (kQ >= 0) ? _mod_tp_sum_sign(y, x, kQ, tau0, ctx) : _mod_tp_sum_sign_pt(y, x, zq, tau0, ctx);
    return *sJ != 0 && sA == -*sJ;
}

/*
    Rebasing an addition: the points y, P = y + x and Q = y - x recorded
    (modulo the half lattice, up to sign; P and Q of type TP_STANDARD),
    x not recorded; or kQ = -1: Q = zq a new point. Without the rebase,
    the values at x would be algebraic over those at y, P, Q (square
    roots: TP_ADD_HALF, TP_HALF), or a new Q = 2y - P of type TP_JACOBI
    from 2y and P (large formulas, and the values at x later algebraic).
    Instead x gets generators of its own (th, the values at x), P becomes
    of type TP_JACOBI from y and x (keeping its generator theta_1), Q of
    type TP_ADD (rational; a new Q is recorded, its values set in thq and
    its record in *kq), and the context records the values of the
    generators of P and Q (_gr_tower_lazy_rebase_add, triggered by
    theta_1(x) or trigger2): the values at the points are then those of
    the order y, x first, whatever the order of evaluation. (Elements
    built on the generators of P and Q before become rational functions
    of the values at y and x: for the theta functions, the values
    themselves.) Returns the record of x, or -1 (no
    change made; th, thq undefined).
*/
static slong
_mod_tp_rebase_add(gr_ptr * th, slong ky, gr_srcptr x0, slong kP, slong kQ, gr_srcptr zq, gr_ptr * thq, slong * kq, ulong trigger2, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_lazy_theta_point_struct es, ej, ea, oldj;
    gr_ptr y, zj, za, nvj[4], nva[4], oldv[4];
    ulong defj[4], defa[4], trigger;
    slong kJ, kA, kx;
    int sJ, k, status, pin_y = 0, pin_j = 0, pin_a = 0;

    GR_TMP_INIT3(y, zj, za, ctx);
    GR_TMP_INIT4(nvj[0], nvj[1], nvj[2], nvj[3], ctx);
    GR_TMP_INIT4(nva[0], nva[1], nva[2], nva[3], ctx);
    GR_TMP_INIT4(oldv[0], oldv[1], oldv[2], oldv[3], ctx);

    kx = -1;
    kJ = kA = -1;
    status = _mod_tp_z(y, ky, ctx);
    if (status != GR_SUCCESS || !_mod_tp_add_check(&kJ, &kA, &sJ, y, ky, x0, kP, kQ, zq, tau0, ctx))
    {
        kJ = kA = -1;
        goto cleanup;
    }

    /* (the values at y, P and Q are read below, after the stores of new
       points, which must not evict them) */
    pin_y = L->theta_pts[ky].pinned;
    pin_j = L->theta_pts[kJ].pinned;
    pin_a = (kA >= 0) ? L->theta_pts[kA].pinned : 0;
    L->theta_pts[ky].pinned = 1;
    L->theta_pts[kJ].pinned = 1;
    if (kA >= 0)
        L->theta_pts[kA].pinned = 1;

    /* the generators of P and Q (theta_1 at P stays) */
    for (k = 0; k < 4; k++)
    {
        defj[k] = (k == 0) ? 0 : _gr_tower_lazy_gen_def_of(&L->theta_pts[kJ].vals[k], ctx);
        defa[k] = (kA < 0) ? 0 : _gr_tower_lazy_gen_def_of(&L->theta_pts[kA].vals[k], ctx);
    }

    /* the generators at x, recorded */
    es.kind = TP_STANDARD;
    es.i = es.j = es.l = -1;
    es.n = es.dv = 0;
    es.s = es.s2 = 1;
    es.have_vals = 0;
    es.stamp = 0;
    es.pinned = 0;
    status = _mod_tp_eval(th, &es, x0, tau0, ctx);
    trigger = (status == GR_SUCCESS) ? _gr_tower_lazy_gen_def_of(th[0], ctx) : 0;
    if (trigger == 0)
        goto cleanup;
    kx = _mod_tp_record(x0, tau0, &es, ctx);
    _mod_tp_store(kx, th, ctx);

    /* P of type TP_JACOBI, y + sJ x */
    ej = L->theta_pts[kJ];
    ej.kind = TP_JACOBI;
    ej.i = ky; ej.j = kx; ej.l = -1;
    ej.n = 1; ej.dv = 1;
    ej.s = 1; ej.s2 = sJ;
    status = _mod_tp_z(zj, kJ, ctx);
    if (status == GR_SUCCESS)
        status = _mod_tp_eval(nvj, &ej, zj, tau0, ctx);
    if (status != GR_SUCCESS)
    {
        /* (kept if not the last record: x of type TP_STANDARD, th its
           values, no rebase) */
        if (_mod_tp_pop(kx, ctx))
            kx = -1;
        goto cleanup;
    }
    oldj = L->theta_pts[kJ];
    for (k = 0; k < 4; k++)
    {
        GR_MUST_SUCCEED(gr_set(oldv[k], &L->theta_pts[kJ].vals[k], ctx));
        GR_MUST_SUCCEED(gr_set(&L->theta_pts[kJ].vals[k], nvj[k], ctx));
    }
    L->theta_pts[kJ].kind = ej.kind;
    L->theta_pts[kJ].i = ej.i; L->theta_pts[kJ].j = ej.j; L->theta_pts[kJ].l = ej.l;
    L->theta_pts[kJ].n = ej.n; L->theta_pts[kJ].dv = ej.dv;
    L->theta_pts[kJ].s = ej.s; L->theta_pts[kJ].s2 = ej.s2;

    /* Q of type TP_ADD, y - sJ x, with the partner P */
    ea = (kA >= 0) ? L->theta_pts[kA] : es;
    ea.kind = TP_ADD;
    ea.i = ky; ea.j = kx; ea.l = kJ;
    ea.n = 0; ea.dv = 0;
    ea.s = 1; ea.s2 = -sJ;
    status = (kA >= 0) ? _mod_tp_z(za, kA, ctx) : gr_set(za, zq, ctx);
    if (status == GR_SUCCESS)
        status = _mod_tp_eval(nva, &ea, za, tau0, ctx);
    if (status != GR_SUCCESS)
    {
        /* (not expected) */
        for (k = 0; k < 4; k++)
            GR_MUST_SUCCEED(gr_set(&L->theta_pts[kJ].vals[k], oldv[k], ctx));
        L->theta_pts[kJ].kind = oldj.kind;
        L->theta_pts[kJ].i = oldj.i; L->theta_pts[kJ].j = oldj.j; L->theta_pts[kJ].l = oldj.l;
        L->theta_pts[kJ].n = oldj.n; L->theta_pts[kJ].dv = oldj.dv;
        L->theta_pts[kJ].s = oldj.s; L->theta_pts[kJ].s2 = oldj.s2;
        if (_mod_tp_pop(kx, ctx))
            kx = -1;
        goto cleanup;
    }
    if (kA >= 0)
    {
        for (k = 0; k < 4; k++)
            GR_MUST_SUCCEED(gr_set(&L->theta_pts[kA].vals[k], nva[k], ctx));
        L->theta_pts[kA].kind = ea.kind;
        L->theta_pts[kA].i = ea.i; L->theta_pts[kA].j = ea.j; L->theta_pts[kA].l = ea.l;
        L->theta_pts[kA].n = ea.n; L->theta_pts[kA].dv = ea.dv;
        L->theta_pts[kA].s = ea.s; L->theta_pts[kA].s2 = ea.s2;
    }
    else
    {
        ea.have_vals = 0;
        *kq = _mod_tp_record(zq, tau0, &ea, ctx);
        _mod_tp_store(*kq, nva, ctx);
        for (k = 0; k < 4; k++)
            gr_swap(thq[k], nva[k], ctx);
    }

    for (k = 0; k < 4; k++)
    {
        if (defj[k] != 0)
            _gr_tower_lazy_rebase_add(defj[k], trigger, trigger2, nvj[k], ctx);
        if (defa[k] != 0)
            _gr_tower_lazy_rebase_add(defa[k], trigger, trigger2, nva[k], ctx);
    }

cleanup:
    if (kJ >= 0)
    {
        L->theta_pts[ky].pinned = pin_y;
        L->theta_pts[kJ].pinned = pin_j;
        if (kA >= 0)
            L->theta_pts[kA].pinned = pin_a;
    }
    GR_TMP_CLEAR3(y, zj, za, ctx);
    GR_TMP_CLEAR4(nvj[0], nvj[1], nvj[2], nvj[3], ctx);
    GR_TMP_CLEAR4(nva[0], nva[1], nva[2], nva[3], ctx);
    GR_TMP_CLEAR4(oldv[0], oldv[1], oldv[2], oldv[3], ctx);
    return kx;
}

/* the record of the reduced point of y + sigma x, or -1 */
static slong
_mod_tp_find_sum(gr_srcptr y, gr_srcptr x, int sigma, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_ptr t, t1;
    slong pp, qq, k = -1;
    int neg;

    GR_TMP_INIT2(t, t1, ctx);
    if (gr_mul_si(t, x, sigma, ctx) == GR_SUCCESS &&
        gr_add(t, y, t, ctx) == GR_SUCCESS &&
        _mod_z_reduce(t, t1, &pp, &qq, &neg, t, tau0, ctx) == GR_SUCCESS &&
        gr_is_zero(t, ctx) == T_FALSE)
        k = _mod_tp_find(t, tau0, ctx);
    GR_TMP_CLEAR2(t, t1, ctx);
    return k;
}

/*
    The rebases of the additions for the new point z0 with the method e
    (_mod_tp_rebase_add):

    TP_ADD_HALF   y, y + z0, y - z0 recorded (the latter two of type
                  TP_STANDARD): z0 the new point x;
    TP_HALF       2 z0 = s (u_i + s2 u_j): z0 gets generators, and
                  w = s u_i - z0 is the new point x for y = z0 (the
                  rebase triggered by theta_1 at either);
    TP_JACOBI     z0 = 2 U + V (n = 2) with V of type TP_STANDARD: the
                  new point x = -(U + V), with y = U, P = -V and Q = z0.

    Sets th to the values at z0 and returns its record, or -1 (no change
    made).
*/
static slong
_mod_tp_rebase_sum(gr_ptr * th, const gr_tower_lazy_theta_point_struct * e, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    slong kz = -1, kP, kQ, kJ, kA, N;
    int sJ;


    if (e->kind == TP_ADD_HALF)
    {
        gr_ptr y;
        GR_TMP_INIT(y, ctx);
        if (_mod_tp_z(y, e->i, ctx) == GR_SUCCESS)
        {
            kP = _mod_tp_find_sum(y, z0, 1, tau0, ctx);
            kQ = _mod_tp_find_sum(y, z0, -1, tau0, ctx);
            kz = _mod_tp_rebase_add(th, e->i, z0, kP, kQ, NULL, NULL, NULL, 0, tau0, ctx);
        }
        GR_TMP_CLEAR(y, ctx);
    }
    else if (e->kind == TP_HALF)
    {
        gr_tower_lazy_theta_point_struct es;
        gr_ptr w, w1, tw[4];
        slong pp, qq;
        int neg, status;
        ulong trig;

        GR_TMP_INIT2(w, w1, ctx);
        GR_TMP_INIT4(tw[0], tw[1], tw[2], tw[3], ctx);

        /* w = s u_i - z0, reduced (z0 + w = s u_i, z0 - w = s s2 u_j
           modulo the half lattice) */
        status = _mod_tp_z(w, e->i, ctx);
        status |= gr_mul_si(w, w, e->s, ctx);
        status |= gr_sub(w, w, z0, ctx);
        if (status == GR_SUCCESS)
            status = _mod_z_reduce(w, w1, &pp, &qq, &neg, w, tau0, ctx);
        if (status == GR_SUCCESS && gr_is_zero(w, ctx) == T_FALSE &&
            _mod_tp_find(w, tau0, ctx) < 0 && !_mod_torsion_point(&N, w, tau0, ctx) &&
            _mod_tp_add_check(&kJ, &kA, &sJ, z0, -1, w, e->i, e->j, NULL, tau0, ctx))
        {
            es.kind = TP_STANDARD;
            es.i = es.j = es.l = -1;
            es.n = es.dv = 0;
            es.s = es.s2 = 1;
            es.have_vals = 0;
            es.stamp = 0;
            es.pinned = 0;
            status = _mod_tp_eval(th, &es, z0, tau0, ctx);
            trig = (status == GR_SUCCESS) ? _gr_tower_lazy_gen_def_of(th[0], ctx) : 0;
            if (trig != 0)
            {
                kz = _mod_tp_record(z0, tau0, &es, ctx);
                _mod_tp_store(kz, th, ctx);
                if (_mod_tp_rebase_add(tw, kz, w, e->i, e->j, NULL, NULL, NULL, trig, tau0, ctx) < 0)
                {
                    /* (not expected; kept, of type TP_STANDARD with the
                       values th, if not the last record) */
                    if (_mod_tp_pop(kz, ctx))
                        kz = -1;
                }
            }
        }

        GR_TMP_CLEAR2(w, w1, ctx);
        GR_TMP_CLEAR4(tw[0], tw[1], tw[2], tw[3], ctx);
    }
    else if (e->kind == TP_JACOBI && e->n == 2)
    {
        /* z0 = 2 U + V (U = s u_i, V = s s2 u_j): z0 = U - x, V = -(U + x)
           for the new point x = -(U + V) */
        gr_ptr x, x1, t, tx[4];
        slong pp, qq, kq = -1;
        int neg, status;

        GR_TMP_INIT3(x, x1, t, ctx);
        GR_TMP_INIT4(tx[0], tx[1], tx[2], tx[3], ctx);
        status = _mod_tp_z(x, e->i, ctx);
        status |= _mod_tp_z(t, e->j, ctx);
        status |= gr_mul_si(t, t, e->s2, ctx);
        status |= gr_add(x, x, t, ctx);
        status |= gr_mul_si(x, x, -e->s, ctx);
        if (status == GR_SUCCESS)
            status = _mod_z_reduce(x, x1, &pp, &qq, &neg, x, tau0, ctx);
        if (status == GR_SUCCESS && gr_is_zero(x, ctx) == T_FALSE &&
            _mod_tp_find(x, tau0, ctx) < 0 && !_mod_torsion_point(&N, x, tau0, ctx) &&
            _mod_tp_rebase_add(tx, e->i, x, e->j, -1, z0, th, &kq, 0, tau0, ctx) >= 0)
            kz = kq;
        GR_TMP_CLEAR3(x, x1, t, ctx);
        GR_TMP_CLEAR4(tx[0], tx[1], tx[2], tx[3], ctx);
    }

    return kz;
}

/*
    The values theta_1, ..., theta_4 at the reduced point z0 != 0, by the
    method recorded for z0, or chosen now (torsion point, a relation to
    recorded points, otherwise generators) and recorded.
*/
static int
_mod_theta_z0(gr_ptr * th, gr_srcptr z0, gr_srcptr tau0, gr_ctx_t ctx)
{
    gr_tower_lazy_theta_point_struct e;
    slong k, N;
    int status;

    k = _mod_tp_find(z0, tau0, ctx);
    if (k >= 0)
    {
        e = LAZY(ctx)->theta_pts[k];
        if (e.have_vals)
        {
            slong q;
            status = GR_SUCCESS;
            for (q = 0; q < 4; q++)
            {
                _gr_tower_lazy_cache_refresh(&LAZY(ctx)->theta_pts[k].vals[q], ctx);
                status |= gr_set(th[q], &LAZY(ctx)->theta_pts[k].vals[q], ctx);
            }
            LAZY(ctx)->theta_pts[k].stamp = ++LAZY(ctx)->cache_clock;
            return status;
        }
        status = _mod_tp_eval(th, &e, z0, tau0, ctx);
        if (status != GR_SUCCESS && e.kind != TP_STANDARD)
        {
            /* (the method recorded failed this time, e.g. a sign test of
               differently represented values: the generators at z0, the
               values given before then independent of these) */
            e.kind = TP_STANDARD;
            status = _mod_tp_eval(th, &e, z0, tau0, ctx);
        }
        if (status == GR_SUCCESS)
            _mod_tp_store(k, th, ctx);
        return status;
    }

    if (_mod_torsion_point(&N, z0, tau0, ctx) && _mod_torsion_field_ok(tau0, ctx))
    {
        e.kind = TP_TORSION;
        e.i = e.j = e.l = -1;
        e.n = N; e.dv = 0; e.s = e.s2 = 1;
    }
    else
    {
        _mod_tp_search(&e, z0, tau0, ctx);
        if (_mod_tp_rebase(th, &e, z0, tau0, ctx) >= 0)
            return GR_SUCCESS;
        if (_mod_tp_rebase_sum(th, &e, z0, tau0, ctx) >= 0)
            return GR_SUCCESS;
    }

    status = _mod_tp_eval(th, &e, z0, tau0, ctx);
    if (status != GR_SUCCESS && e.kind != TP_STANDARD)
    {
        e.kind = TP_STANDARD;
        status = _mod_tp_eval(th, &e, z0, tau0, ctx);
    }
    if (status == GR_SUCCESS)
    {
        e.have_vals = 0;
        k = _mod_tp_record(z0, tau0, &e, ctx);
        _mod_tp_store(k, th, ctx);
    }
    return status;
}


/*
    z' = -z / (c tau + d) for c != 0, otherwise z (the convention of
    acb_modular_theta_transform)
*/
static int
_mod_zprime(gr_ptr zp, const psl2z_t g, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
{
    int status;
    if (fmpz_is_zero(&g->c))
        return gr_set(zp, z, ctx);
    status = _mod_cd(zp, g, tau, ctx);
    status |= gr_div(zp, z, zp, ctx);
    status |= gr_neg(zp, zp, ctx);
    return status;
}

/* the order of the reduced points: Im, then Re */
static int
_mod_z_cmp(int * sgn, gr_srcptr a, gr_srcptr b, gr_ctx_t ctx)
{
    gr_ptr d, t;
    int status;
    GR_TMP_INIT2(d, t, ctx);
    status = gr_sub(d, a, b, ctx);
    status |= gr_im(t, d, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_real_sign_fast(sgn, t, ctx);
    if (status == GR_SUCCESS && *sgn == 0)
    {
        status = gr_re(t, d, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_real_sign_fast(sgn, t, ctx);
    }
    GR_TMP_CLEAR2(d, t, ctx);
    return status;
}

/*
    For tau0 = g tau in {i, rho}: g := h g with h in the stabilizer of
    tau0 giving the least reduced z0 (so that theta(z, tau) is expressed
    through the same generators for every z in the orbit).
*/
static int
_mod_stabilizer_choice(psl2z_t g, gr_srcptr z, gr_srcptr tau, gr_srcptr tau0, gr_ctx_t ctx)
{
    psl2z_t h, gh, best;
    gr_ptr t, zp, z0, z1, zbest;
    slong p, q, n = 1, s;
    int neg, status = GR_SUCCESS, sgn;
    truth_t eq;

    GR_TMP_INIT5(t, zp, z0, z1, zbest, ctx);
    psl2z_init(h);
    psl2z_init(gh);
    psl2z_init(best);

    /* tau0 = i: h = S; tau0 = rho: h = [[1, -1], [1, 0]] (order 3) */
    status = gr_i(t, ctx);
    eq = gr_equal(tau0, t, ctx);
    if (eq == T_TRUE)
    {
        n = 2;
        fmpz_zero(&h->a); fmpz_set_si(&h->b, -1);
        fmpz_one(&h->c); fmpz_zero(&h->d);
    }
    else if (eq == T_FALSE)
    {
        status = gr_set_si(t, -3, ctx);
        status |= gr_sqrt(t, t, ctx);
        status |= gr_add_ui(t, t, 1, ctx);
        status |= gr_div_ui(t, t, 2, ctx);
        eq = (status == GR_SUCCESS) ? gr_equal(tau0, t, ctx) : T_UNKNOWN;
        if (eq == T_TRUE)
        {
            n = 3;
            fmpz_one(&h->a); fmpz_set_si(&h->b, -1);
            fmpz_one(&h->c); fmpz_zero(&h->d);
        }
    }
    if (eq == T_UNKNOWN)
        status = GR_UNABLE;

    if (status == GR_SUCCESS && n > 1)
    {
        psl2z_set(gh, g);
        psl2z_set(best, g);
        for (s = 0; s < n && status == GR_SUCCESS; s++)
        {
            if (s > 0)
                _mod_left_mul(gh, h);
            status = _mod_zprime(zp, gh, z, tau, ctx);
            if (status == GR_SUCCESS)
                status = _mod_z_reduce(z0, z1, &p, &q, &neg, zp, tau0, ctx);
            if (status == GR_SUCCESS)
            {
                sgn = -1;
                if (s > 0)
                    status = _mod_z_cmp(&sgn, z0, zbest, ctx);
                if (status == GR_SUCCESS && sgn < 0)
                {
                    psl2z_set(best, gh);
                    status = gr_set(zbest, z0, ctx);
                }
            }
        }
        if (status == GR_SUCCESS)
            psl2z_set(g, best);
    }

    GR_TMP_CLEAR5(t, zp, z0, z1, zbest, ctx);
    psl2z_clear(h);
    psl2z_clear(gh);
    psl2z_clear(best);
    return status;
}

/*
    theta_k(z, tau): with tau0 = g tau, theta_{1+i}(z, tau) = exp(pi i R_i / 4)
    A B theta_{1+S_i}(z', tau0), z' = -z / (c tau + d), A = sqrt(i / (c tau + d)),
    B = exp(-pi i c z^2 / (c tau + d)) (z' = z, A = B = 1 if C = 0).
*/
static int
_mod_theta_z(gr_ptr * out, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
{
    psl2z_t g;
    gr_ptr tau0, zp, cd, A, B, u, th[4];
    int R[4], S[4], C, i;
    int status;

    psl2z_init(g);
    GR_TMP_INIT5(tau0, zp, cd, A, B, ctx);
    GR_TMP_INIT(u, ctx);
    GR_TMP_INIT4(th[0], th[1], th[2], th[3], ctx);

    status = _gr_tower_lazy_modular_reduce(tau0, g, tau, ctx);

    /* at tau0 = i and rho, the stabilizer acts on z0 (by z0 -> -i z0, and
       z0 -> rho^-1 z0, rho^-2 z0): the canonical z0 of the orbit */
    if (status == GR_SUCCESS)
        status = _mod_stabilizer_choice(g, z, tau, tau0, ctx);

    if (status == GR_SUCCESS)
    {
        acb_modular_theta_transform(R, S, &C, g);
        if (C)
        {
            status = _mod_cd(cd, g, tau, ctx);
            status |= gr_div(zp, z, cd, ctx);
            status |= gr_neg(zp, zp, ctx);
            status |= gr_inv(A, cd, ctx);
            status |= gr_i(u, ctx);
            status |= gr_mul(A, A, u, ctx);
            status |= gr_sqrt(A, A, ctx);
            /* B = exp(-pi i c z^2 / (c tau + d)) */
            status |= gr_sqr(B, z, ctx);
            status |= gr_div(B, B, cd, ctx);
            status |= gr_mul_fmpz(B, B, &g->c, ctx);
            status |= gr_neg(B, B, ctx);
            status |= _mod_exp_pi_i(B, B, ctx);
            status |= gr_mul(A, A, B, ctx);
        }
        else
        {
            status = gr_set(zp, z, ctx);
            status |= gr_one(A, ctx);
        }
    }

    if (status == GR_SUCCESS)
        status = _mod_theta_z_reduced(th, zp, tau0, ctx);

    for (i = 0; i < 4 && status == GR_SUCCESS; i++)
    {
        if (out[i] == NULL)
            continue;
        status = _mod_root_of_unity(u, R[i], 4, ctx);
        status |= gr_mul(u, u, A, ctx);
        status |= gr_mul(out[i], u, th[S[i]], ctx);
    }

    psl2z_clear(g);
    GR_TMP_CLEAR5(tau0, zp, cd, A, B, ctx);
    GR_TMP_CLEAR(u, ctx);
    GR_TMP_CLEAR4(th[0], th[1], th[2], th[3], ctx);
    return status;
}

/* theta_j(z, tau), for the definition of a generator (conjugation and
   transfer between contexts) */
int
_gr_tower_lazy_jacobi_theta_j(gr_ptr res, slong j, gr_srcptr z_in, gr_srcptr tau_in, gr_ctx_t ctx)
{
    gr_ptr out[4], z, tau;
    int status;
    _gr_tower_lazy_lock(ctx);
    GR_TMP_INIT2(z, tau, ctx);
    out[0] = out[1] = out[2] = out[3] = NULL;
    status = gr_set(z, z_in, ctx);
    status |= gr_set(tau, tau_in, ctx);
    if (j < 1 || j > 4)
        status = GR_DOMAIN;
    else if (status == GR_SUCCESS)
    {
        out[j - 1] = res;
        status = _mod_theta_z(out, z, tau, ctx);
    }
    GR_TMP_CLEAR2(z, tau, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* the theta functions */
int
gr_tower_lazy_jacobi_theta(gr_tower_lazy_elem_t res1, gr_tower_lazy_elem_t res2, gr_tower_lazy_elem_t res3,
    gr_tower_lazy_elem_t res4, const gr_tower_lazy_elem_t z_in, const gr_tower_lazy_elem_t tau_in, gr_ctx_t ctx)
{
    int status, real, alg;
    truth_t zero;
    gr_ptr z, tau;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    GR_TMP_INIT2(z, tau, ctx);
    if (gr_set(z, z_in, ctx) != GR_SUCCESS || gr_set(tau, tau_in, ctx) != GR_SUCCESS)
        zero = T_UNKNOWN;
    else
        zero = gr_is_zero(z, ctx);
    if (zero == T_TRUE)
    {
        status = _mod_theta_consts(res2, res3, res4, tau, ctx);
        status |= gr_zero(res1, ctx);
    }
    else if (zero == T_FALSE)
    {
        gr_ptr out[4];
        out[0] = res1; out[1] = res2; out[2] = res3; out[3] = res4;
        status = _mod_theta_z(out, z, tau, ctx);
    }
    else
        status = GR_UNABLE;
    if (status == GR_SUCCESS)
    {
        status = _gr_tower_lazy_view_finish(status, res1, real, alg, ctx);
        status |= _gr_tower_lazy_view_finish(GR_SUCCESS, res2, real, alg, ctx);
        status |= _gr_tower_lazy_view_finish(GR_SUCCESS, res3, real, alg, ctx);
        status |= _gr_tower_lazy_view_finish(GR_SUCCESS, res4, real, alg, ctx);
    }
    GR_TMP_CLEAR2(z, tau, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

#define THETA_WRAP(name, which) \
int name(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t tau, gr_ctx_t ctx) \
{ \
    gr_ptr t[4]; \
    int status, i; \
    _gr_tower_lazy_lock(ctx); \
    for (i = 0; i < 4; i++) \
        t[i] = gr_heap_init(ctx); \
    status = gr_tower_lazy_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, ctx); \
    if (status == GR_SUCCESS) \
        status = gr_set(res, t[which], ctx); \
    for (i = 0; i < 4; i++) \
        gr_heap_clear(t[i], ctx); \
    _gr_tower_lazy_unlock(ctx); \
    return status; \
}

THETA_WRAP(gr_tower_lazy_jacobi_theta_1, 0)
THETA_WRAP(gr_tower_lazy_jacobi_theta_2, 1)
THETA_WRAP(gr_tower_lazy_jacobi_theta_3, 2)
THETA_WRAP(gr_tower_lazy_jacobi_theta_4, 3)

/* -------------------------------------------------------------------- */
/* Weierstrass functions                                                 */
/* -------------------------------------------------------------------- */

/*
    For the lattice Z + tau Z, with theta_k = theta_k(0, tau) and
    theta_k(z) = theta_k(z, tau):

        e_1 = pi^2 (theta_3^4 + theta_4^4) / 3,  e_2 = pi^2 (theta_2^4 - theta_4^4) / 3,
        e_3 = -pi^2 (theta_2^4 + theta_3^4) / 3,
        g_2 = 4 pi^4 E_4 / 3,  g_3 = 8 pi^6 E_6 / 27,
        wp(z) = e_1 + (pi theta_3 theta_4 theta_2(z) / theta_1(z))^2,
        wp'(z) = -2 pi^3 (theta_2 theta_3 theta_4)^2 theta_2(z) theta_3(z) theta_4(z) / theta_1(z)^3,
        sigma(z) = theta_1(z) exp(pi^2 E_2 z^2 / 6) / (pi theta_2 theta_3 theta_4)

    (the theta functions of z and E_2 from the theta layer, so that the
    algebraic relations between wp, wp', the e_k and g_k hold by
    construction). The Weierstrass zeta function needs theta_1'(z) /
    theta_1(z), which is not a theta value: it is not implemented, nor is
    the inverse of wp.
*/

/* e_1, e_2, e_3 (any may be NULL) */
static int
_mod_roots(gr_ptr e1, gr_ptr e2, gr_ptr e3, gr_srcptr tau, gr_ctx_t ctx)
{
    gr_ptr a2, a3, a4, p, u;
    int status;

    GR_TMP_INIT5(a2, a3, a4, p, u, ctx);
    status = _mod_theta_consts(a2, a3, a4, tau, ctx);
    status |= gr_pow_ui(a2, a2, 4, ctx);
    status |= gr_pow_ui(a3, a3, 4, ctx);
    status |= gr_pow_ui(a4, a4, 4, ctx);
    status |= gr_pi(p, ctx);
    status |= gr_sqr(p, p, ctx);
    status |= gr_div_ui(p, p, 3, ctx);
    if (status == GR_SUCCESS && e1 != NULL)
    {
        status = gr_add(u, a3, a4, ctx);
        status |= gr_mul(e1, u, p, ctx);
    }
    if (status == GR_SUCCESS && e2 != NULL)
    {
        status = gr_sub(u, a2, a4, ctx);
        status |= gr_mul(e2, u, p, ctx);
    }
    if (status == GR_SUCCESS && e3 != NULL)
    {
        status = gr_add(u, a2, a3, ctx);
        status |= gr_mul(e3, u, p, ctx);
        status |= gr_neg(e3, e3, ctx);
    }
    GR_TMP_CLEAR5(a2, a3, a4, p, u, ctx);
    return status;
}

/* which: 0 (wp), 1 (wp'), 2 (sigma) */
static int
_mod_weierstrass(gr_ptr res, int which, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx)
{
    gr_ptr th[4], a2, a3, a4, p, u;
    truth_t zero;
    int status;

    zero = gr_is_zero(z, ctx);
    if (zero == T_UNKNOWN)
        return GR_UNABLE;
    if (zero == T_TRUE)
        return (which == 2) ? gr_zero(res, ctx) : GR_DOMAIN;

    GR_TMP_INIT4(th[0], th[1], th[2], th[3], ctx);
    GR_TMP_INIT5(a2, a3, a4, p, u, ctx);

    status = _mod_theta_z(th, z, tau, ctx);
    status |= _mod_theta_consts(a2, a3, a4, tau, ctx);
    status |= gr_pi(p, ctx);

    /* (theta_1(z) = 0 at the lattice points: a pole) */
    if (status == GR_SUCCESS && which != 2)
    {
        truth_t t = gr_is_zero(th[0], ctx);
        if (t == T_TRUE)
            status = GR_DOMAIN;
        else if (t != T_FALSE)
            status = GR_UNABLE;
    }

    if (status == GR_SUCCESS && which == 0)
    {
        status = _mod_roots(res, NULL, NULL, tau, ctx);
        status |= gr_mul(u, p, a3, ctx);
        status |= gr_mul(u, u, a4, ctx);
        status |= gr_mul(u, u, th[1], ctx);
        status |= gr_div(u, u, th[0], ctx);
        status |= gr_sqr(u, u, ctx);
        status |= gr_add(res, res, u, ctx);
    }
    else if (status == GR_SUCCESS && which == 1)
    {
        status = gr_mul(u, a2, a3, ctx);
        status |= gr_mul(u, u, a4, ctx);
        status |= gr_sqr(u, u, ctx);
        status |= gr_pow_ui(a2, p, 3, ctx);
        status |= gr_mul(u, u, a2, ctx);
        status |= gr_mul_si(u, u, -2, ctx);
        status |= gr_mul(u, u, th[1], ctx);
        status |= gr_mul(u, u, th[2], ctx);
        status |= gr_mul(u, u, th[3], ctx);
        status |= gr_pow_ui(a2, th[0], 3, ctx);
        status |= gr_div(res, u, a2, ctx);
    }
    else if (status == GR_SUCCESS)
    {
        /* theta_1(z) exp(pi^2 E_2 z^2 / 6) / (pi theta_2 theta_3 theta_4) */
        status = _mod_eisenstein(u, 2, tau, ctx);
        status |= gr_mul(u, u, p, ctx);
        status |= gr_mul(u, u, p, ctx);
        status |= gr_mul(u, u, z, ctx);
        status |= gr_mul(u, u, z, ctx);
        status |= gr_div_ui(u, u, 6, ctx);
        status |= gr_exp(u, u, ctx);
        status |= gr_mul(u, u, th[0], ctx);
        status |= gr_mul(p, p, a2, ctx);
        status |= gr_mul(p, p, a3, ctx);
        status |= gr_mul(p, p, a4, ctx);
        status |= gr_div(res, u, p, ctx);
    }

    GR_TMP_CLEAR4(th[0], th[1], th[2], th[3], ctx);
    GR_TMP_CLEAR5(a2, a3, a4, p, u, ctx);
    return status;
}

#define WEIERSTRASS_WRAP(name, which) \
int name(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z_in, const gr_tower_lazy_elem_t tau_in, gr_ctx_t ctx) \
{ \
    int status, real, alg; \
    gr_ptr z, tau; \
    _gr_tower_lazy_lock(ctx); \
    real = REAL(ctx); \
    alg = ALG(ctx); \
    GR_TMP_INIT2(z, tau, ctx); \
    status = gr_set(z, z_in, ctx); \
    status |= gr_set(tau, tau_in, ctx); \
    if (status == GR_SUCCESS) \
        status = _gr_tower_lazy_view_finish(_mod_weierstrass(res, which, z, tau, ctx), res, real, alg, ctx); \
    GR_TMP_CLEAR2(z, tau, ctx); \
    _gr_tower_lazy_unlock(ctx); \
    return status; \
}

WEIERSTRASS_WRAP(gr_tower_lazy_weierstrass_p, 0)
WEIERSTRASS_WRAP(gr_tower_lazy_weierstrass_p_prime, 1)
WEIERSTRASS_WRAP(gr_tower_lazy_weierstrass_sigma, 2)

int
gr_tower_lazy_elliptic_roots(gr_tower_lazy_elem_t e1, gr_tower_lazy_elem_t e2, gr_tower_lazy_elem_t e3,
    const gr_tower_lazy_elem_t tau_in, gr_ctx_t ctx)
{
    int status, real, alg;
    gr_ptr tau;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    GR_TMP_INIT(tau, ctx);
    status = gr_set(tau, tau_in, ctx);
    if (status == GR_SUCCESS)
        status = _mod_roots(e1, e2, e3, tau, ctx);
    if (status == GR_SUCCESS)
    {
        status = _gr_tower_lazy_view_finish(GR_SUCCESS, e1, real, alg, ctx);
        status |= _gr_tower_lazy_view_finish(GR_SUCCESS, e2, real, alg, ctx);
        status |= _gr_tower_lazy_view_finish(GR_SUCCESS, e3, real, alg, ctx);
    }
    GR_TMP_CLEAR(tau, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* g_2 = 4 pi^4 E_4 / 3, g_3 = 8 pi^6 E_6 / 27 */
int
gr_tower_lazy_elliptic_invariants(gr_tower_lazy_elem_t g2, gr_tower_lazy_elem_t g3,
    const gr_tower_lazy_elem_t tau_in, gr_ctx_t ctx)
{
    int status, real, alg;
    gr_ptr tau, p;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    GR_TMP_INIT2(tau, p, ctx);
    status = gr_set(tau, tau_in, ctx);
    status |= gr_pi(p, ctx);
    status |= gr_sqr(p, p, ctx);
    if (status == GR_SUCCESS)
        status = _mod_eisenstein(g2, 4, tau, ctx);
    if (status == GR_SUCCESS)
        status = _mod_eisenstein(g3, 6, tau, ctx);
    if (status == GR_SUCCESS)
    {
        gr_ptr q;
        GR_TMP_INIT(q, ctx);
        status = gr_sqr(q, p, ctx);
        status |= gr_mul(g2, g2, q, ctx);
        status |= gr_mul_ui(g2, g2, 4, ctx);
        status |= gr_div_ui(g2, g2, 3, ctx);
        status |= gr_mul(q, q, p, ctx);
        status |= gr_mul(g3, g3, q, ctx);
        status |= gr_mul_ui(g3, g3, 8, ctx);
        status |= gr_div_ui(g3, g3, 27, ctx);
        GR_TMP_CLEAR(q, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = _gr_tower_lazy_view_finish(GR_SUCCESS, g2, real, alg, ctx);
        status |= _gr_tower_lazy_view_finish(GR_SUCCESS, g3, real, alg, ctx);
    }
    GR_TMP_CLEAR2(tau, p, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* the generator lambda(tau) for tau in the fundamental domain, with
   canonicalization (conjugation and transfer of generators) */
int
_gr_tower_lazy_modular_lambda(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t tau_in, gr_ctx_t ctx)
{
    int status;
    gr_ptr tau;
    _gr_tower_lazy_lock(ctx);
    GR_TMP_INIT(tau, ctx);
    status = gr_set(tau, tau_in, ctx);
    if (status == GR_SUCCESS)
        status = _mod_lambda(res, tau, ctx);
    GR_TMP_CLEAR(tau, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* periods at CM points (Chowla-Selberg)                                 */
/* -------------------------------------------------------------------- */

/* the fundamental discriminants of class number one */
static int
_mod_cm_fundamental_h1(slong d)
{
    return d == -3 || d == -4 || d == -7 || d == -8 || d == -11 ||
           d == -19 || d == -43 || d == -67 || d == -163;
}

/*
    theta_3(tau0)^2 at the CM point tau0 in F of fundamental discriminant
    d with h(d) = 1 (tau0 = (b + sqrt(d))/2, b = 0, 1), with lam0 =
    lambda(tau0): by the Chowla-Selberg formula
    Delta(tau0) Im(tau0)^6 = (4 pi sqrt|d|)^(-6) prod_j Gamma(j/|d|)^(3 w chi(j))
    (up to the sign of Delta), so that
    |eta(tau0)|^2 = (4 pi sqrt|d|)^(-1/2) prod Gamma(j/|d|)^(w chi(j)/4) (Im tau0)^(-1/2),
    eta(tau0)^2 = exp(pi i Re(tau0) / 6) |eta(tau0)|^2 (q is real), and
    theta_3^2 = 2^(2/3) eta^2 (lam0 (1 - lam0))^(-1/6) up to a cube root of
    unity, chosen numerically.
*/
static int
_mod_cm_theta3_sq(gr_ptr res, gr_srcptr tau0, gr_srcptr lam0, slong d, gr_ctx_t ctx)
{
    slong w = (d == -3) ? 6 : (d == -4) ? 4 : 2, j, n = -d;
    gr_ptr t, u, G;
    fmpq_t e;
    fmpz_t fd, fj;
    int status = GR_SUCCESS;

    GR_TMP_INIT3(t, u, G, ctx);
    fmpq_init(e);
    fmpz_init_set_si(fd, d);
    fmpz_init(fj);

    /* G = prod Gamma(j/n)^(w chi(j)/4) */
    status |= gr_one(G, ctx);
    for (j = 1; j < n && status == GR_SUCCESS; j++)
    {
        int c;
        fmpz_set_si(fj, j);
        c = fmpz_kronecker(fd, fj);
        if (c == 0)
            continue;
        fmpq_set_si(e, j, n);
        status = gr_set_fmpq(t, e, ctx);
        status |= gr_gamma(t, t, ctx);
        fmpq_set_si(e, w * c, 4);
        status |= gr_pow_fmpq(t, t, e, ctx);
        status |= gr_mul(G, G, t, ctx);
    }

    /* |eta|^2 = (4 pi sqrt(n))^(-1/2) G (sqrt(n)/2)^(-1/2) = G / sqrt(2 pi n) */
    if (status == GR_SUCCESS)
    {
        status = gr_pi(t, ctx);
        status |= gr_mul_si(t, t, 2 * n, ctx);
        status |= gr_rsqrt(t, t, ctx);
        status |= gr_mul(res, G, t, ctx);
    }

    /* exp(pi i Re(tau0) / 6): Re(tau0) = 0 or 1/2 */
    if (status == GR_SUCCESS && (n % 2 == 1))
    {
        status = _mod_root_of_unity(t, 1, 12, ctx);
        status |= gr_mul(res, res, t, ctx);
    }

    /* theta_3^2 = 2^(2/3) eta^2 (lam0 (1 - lam0))^(-1/6) */
    if (status == GR_SUCCESS)
    {
        fmpq_set_si(e, 2, 3);
        status = gr_set_ui(t, 2, ctx);
        status |= gr_pow_fmpq(t, t, e, ctx);
        status |= gr_mul(res, res, t, ctx);
        status |= gr_sub_ui(t, lam0, 1, ctx);
        status |= gr_neg(t, t, ctx);
        status |= gr_mul(t, t, lam0, ctx);
        fmpq_set_si(e, 1, 6);
        status |= gr_pow_fmpq(t, t, e, ctx);
        status |= gr_div(res, res, t, ctx);
    }

    /* the sixth root of unity (theta_3^12 = 16 eta^12 / (lam0 (1 - lam0))
       by 2 eta^3 = theta_2 theta_3 theta_4): theta_3(tau0)^2 numerically,
       the one root overlapping */
    if (status == GR_SUCCESS)
    {
        acb_t x, y, th2, th3, th4, th1, tz;
        slong k, found = -1;
        acb_init(x); acb_init(y); acb_init(th1); acb_init(th2); acb_init(th3); acb_init(th4); acb_init(tz);
        if (gr_tower_lazy_get_acb(x, res, 128, ctx) == GR_SUCCESS &&
            gr_tower_lazy_get_acb(tz, tau0, 128, ctx) == GR_SUCCESS)
        {
            acb_t zero;
            acb_init(zero);
            acb_modular_theta(th1, th2, th3, th4, zero, tz, 128);
            acb_sqr(th3, th3, 128);
            for (k = 0; k < 6; k++)
            {
                acb_set_si(y, k);
                acb_div_ui(y, y, 3, 128);
                acb_exp_pi_i(y, y, 128);
                acb_mul(y, y, x, 128);
                if (acb_overlaps(y, th3))
                {
                    if (found >= 0)
                        found = -2;
                    else if (found == -1)
                        found = k;
                }
            }
            acb_clear(zero);
        }
        if (found >= 1)
        {
            status = _mod_root_of_unity(t, found, 3, ctx);
            status |= gr_mul(res, res, t, ctx);
        }
        else if (found < 0)
            status = GR_UNABLE;
        acb_clear(x); acb_clear(y); acb_clear(th1); acb_clear(th2); acb_clear(th3); acb_clear(th4); acb_clear(tz);
    }

    GR_TMP_CLEAR3(t, u, G, ctx);
    fmpq_clear(e);
    fmpz_clear(fd);
    fmpz_clear(fj);
    return status;
}

/*
    K(m) at a CM modulus: m = lambda(tau_m), tau_m = i K(1 - m) / K(m), at a
    CM point of a discriminant listed above (found numerically, then
    verified exactly: tau_m = g^(-1) tau0 with tau0 the exact quadratic
    point in F and lambda(tau_m) = m exactly), K(m) = (pi/2)
    theta_3(tau_m)^2 by the transformation of theta_3^2 from tau0. A
    final numerical check guards the choice of tau_m (the other points of
    its orbit give other values). Sets *done.
*/
static int
_mod_elliptic_cm_K(gr_ptr res, int * done, gr_srcptr m, gr_ptr taum_out, slong * dout, gr_ctx_t ctx)
{
    acb_t mz, k1, k2, tz, w;
    arf_t eps;
    psl2z_t g;
    slong prec = 128;
    int status = GR_SUCCESS;

    *done = 0;
    if (_gr_tower_lazy_is_algebraic_repr_locked(m, ctx) != T_TRUE)
        return GR_SUCCESS;

    acb_init(mz); acb_init(k1); acb_init(k2); acb_init(tz); acb_init(w);
    arf_init(eps);
    psl2z_init(g);

    if (gr_tower_lazy_get_acb(mz, m, prec, ctx) == GR_SUCCESS)
    {
        /* tau_m = i K(1 - m) / K(m), reduced */
        acb_sub_ui(k1, mz, 1, prec);
        acb_neg(k1, k1);
        acb_elliptic_k(k1, k1, prec);
        acb_elliptic_k(k2, mz, prec);
        acb_div(tz, k1, k2, prec);
        acb_mul_onei(tz, tz);
        arf_set_d(eps, 1.0 - 1.0 / 64);

        if (acb_is_finite(tz) && arb_is_positive(acb_imagref(tz)))
        {
            acb_modular_fundamental_domain_approx(w, g, tz, eps, prec);

            /* tau0 approximately (-B + sqrt(D)) / (2A): A = 1 for h = 1 */
            {
                double re = arf_get_d(arb_midref(acb_realref(w)), ARF_RND_NEAR);
                double im = arf_get_d(arb_midref(acb_imagref(w)), ARF_RND_NEAR);
                slong d = (slong) floor(-4.0 * im * im + 0.5) - 0;
                slong b = (fabs(fabs(re) - 0.5) < 1e-8) ? 1 : (fabs(re) < 1e-8) ? 0 : -1;
                if (b == 1 && re < 0)
                {
                    /* (Re(tau0) = 1/2 in F: tau0 = w + 1) */
                    psl2z_t T1;
                    psl2z_init(T1);
                    fmpz_one(&T1->b);
                    psl2z_mul(g, T1, g);
                    psl2z_clear(T1);
                    re += 1.0;
                }
                /* D = b^2 - 4 c, |tau0|^2 = c */
                if (b >= 0)
                    d = b * b - 4 * (slong) floor(re * re + im * im + 0.5);
                /* (the Gamma values at level |d| in the normal form of
                   the gamma lattice: otherwise K(m) is a generator, which
                   is cheaper than as many Gamma generators) */
                if (b >= 0 && _mod_cm_fundamental_h1(d) && fabs(4.0 * im * im - (double) (-d)) < 1e-6 &&
                    -d <= gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_GAMMA_LATTICE_LIMIT))
                {
                    gr_ptr tau0, taum, lam0, lamm, th[4], t, A;
                    psl2z_t ginv;
                    int R[4], S[4], C;

                    GR_TMP_INIT4(tau0, taum, lam0, lamm, ctx);
                    GR_TMP_INIT4(th[0], th[1], th[2], th[3], ctx);
                    GR_TMP_INIT2(t, A, ctx);
                    psl2z_init(ginv);

                    /* tau0 = (b + i sqrt(-d)) / 2 */
                    status = gr_set_si(t, -d, ctx);
                    status |= gr_sqrt(t, t, ctx);
                    status |= gr_i(A, ctx);
                    status |= gr_mul(t, t, A, ctx);
                    status |= gr_add_si(t, t, b, ctx);
                    status |= gr_div_ui(tau0, t, 2, ctx);

                    /* tau_m = g^(-1) tau0 */
                    psl2z_inv(ginv, g);
                    status |= _mod_apply(taum, ginv, tau0, ctx);

                    /* lambda(tau_m) = m exactly */
                    if (status == GR_SUCCESS)
                        status = _mod_lambda(lamm, taum, ctx);
                    if (status == GR_SUCCESS && gr_equal(lamm, m, ctx) == T_TRUE)
                    {
                        psl2z_t h;
                        psl2z_init(h);

                        /* tau0 = h tau_m with the exact reduction, theta^2 at tau0 */
                        status = _gr_tower_lazy_modular_reduce(t, h, taum, ctx);
                        if (status == GR_SUCCESS)
                            status = _mod_lambda0(lam0, t, ctx);
                        if (status == GR_SUCCESS)
                            status = _mod_cm_theta3_sq(th[2], t, lam0, d, ctx);
                        if (status == GR_SUCCESS)
                        {
                            fmpq_t e;
                            fmpq_init(e);
                            fmpq_set_si(e, 1, 2);
                            status = gr_zero(th[0], ctx);
                            status |= gr_pow_fmpq(th[1], lam0, e, ctx);
                            status |= gr_mul(th[1], th[1], th[2], ctx);
                            status |= gr_sub_ui(th[3], lam0, 1, ctx);
                            status |= gr_neg(th[3], th[3], ctx);
                            status |= gr_pow_fmpq(th[3], th[3], e, ctx);
                            status |= gr_mul(th[3], th[3], th[2], ctx);
                            fmpq_clear(e);
                        }

                        /* theta_3(tau_m)^2 = exp(pi i R_2 / 2) A^2 theta_{S_2}(tau0)^2,
                           A^2 = i / (c tau_m + d) */
                        if (status == GR_SUCCESS)
                        {
                            acb_modular_theta_transform(R, S, &C, h);
                            status = _mod_root_of_unity(res, R[2], 2, ctx);
                            status |= gr_mul(res, res, th[S[2]], ctx);
                            if (C)
                            {
                                status |= _mod_cd(A, h, taum, ctx);
                                status |= gr_div(res, res, A, ctx);
                                status |= gr_i(A, ctx);
                                status |= gr_mul(res, res, A, ctx);
                            }
                            status |= gr_pi(A, ctx);
                            status |= gr_mul(res, res, A, ctx);
                            status |= gr_div_ui(res, res, 2, ctx);
                        }

                        /* the numerical check */
                        if (status == GR_SUCCESS)
                        {
                            acb_t x;
                            acb_init(x);
                            if (gr_tower_lazy_get_acb(x, res, prec, ctx) == GR_SUCCESS && acb_overlaps(x, k2))
                            {
                                *done = 1;
                                *dout = d;
                                if (taum_out != NULL)
                                    status = gr_set(taum_out, taum, ctx);
                            }
                            acb_clear(x);
                        }
                        psl2z_clear(h);
                    }

                    /* (a failure here is not an error: K(m) is a generator) */
                    if (!*done)
                        status = GR_SUCCESS;

                    psl2z_clear(ginv);
                    GR_TMP_CLEAR4(tau0, taum, lam0, lamm, ctx);
                    GR_TMP_CLEAR4(th[0], th[1], th[2], th[3], ctx);
                    GR_TMP_CLEAR2(t, A, ctx);
                }
            }
        }
    }

    acb_clear(mz); acb_clear(k1); acb_clear(k2); acb_clear(tz); acb_clear(w);
    arf_clear(eps);
    psl2z_clear(g);
    return status;
}

/*
    E(m) at the CM moduli of discriminant -3 and -4 (m = 1/2, -1, 2,
    exp(+-pi i/3)): with E2 = theta_3^4 (3 E / K - 2 + m) at tau_m and
    E2(tau_m) = 3 / (pi Im(tau_m)) (the nonholomorphic E2* vanishes on the
    orbits of i and rho, being of weight 2 under their stabilizers),

        E(m) = K(m) (2 - m) / 3 + pi / (4 Im(tau_m) K(m)),

    so that E2 there is the rational multiple of 1/pi which the
    transformations of E2 under the stabilizers of i, rho require.
*/
int
_gr_tower_lazy_elliptic_cm(gr_ptr res, int * done, gr_srcptr m, int kind, gr_ctx_t ctx)
{
    gr_ptr K, taum, t, u;
    slong d = 0;
    int status;

    if (kind == GR_TOWER_ELLIPTIC_K)
        return _mod_elliptic_cm_K(res, done, m, NULL, &d, ctx);

    *done = 0;
    if (kind != GR_TOWER_ELLIPTIC_E)
        return GR_SUCCESS;

    GR_TMP_INIT4(K, taum, t, u, ctx);
    status = _mod_elliptic_cm_K(K, done, m, taum, &d, ctx);
    if (status == GR_SUCCESS && *done && (d == -3 || d == -4))
    {
        status = gr_sub_ui(t, m, 2, ctx);
        status |= gr_neg(t, t, ctx);
        status |= gr_mul(t, t, K, ctx);
        status |= gr_div_ui(t, t, 3, ctx);
        status |= gr_im(u, taum, ctx);
        status |= gr_mul(u, u, K, ctx);
        status |= gr_mul_ui(u, u, 4, ctx);
        status |= gr_pi(res, ctx);
        status |= gr_div(res, res, u, ctx);
        status |= gr_add(res, res, t, ctx);
        if (status != GR_SUCCESS)
        {
            /* (a failure is not an error: E(m) is a generator) */
            *done = 0;
            status = GR_SUCCESS;
        }
    }
    else
        *done = 0;
    GR_TMP_CLEAR4(K, taum, t, u, ctx);
    return status;
}

POP_OPTIONS
