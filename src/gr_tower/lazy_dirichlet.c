/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Dirichlet L-functions, the Hurwitz zeta function at rational
    parameters and the Lerch transcendent in the lazy field.

    Normal form: L(s, chi) for a Dirichlet character chi mod q is
    L(s, chi*) prod_{p | q} (1 - chi*(p) p^(-s)) for the primitive
    character chi* of conductor f inducing chi (zeta(s) for f = 1). At an
    integer s, L(s, chi*) = f^(-s) sum_{a <= f} chi*(a) zeta(s, a/f) in
    terms of the Hurwitz (polygamma) normal form, with -(1/f) sum chi*(a)
    psi(a/f) at s = 1 and Bernoulli polynomials at s <= 0. At other s,
    the functional equation

        L(s, chi) = eps(chi) (f/pi)^(1/2 - s) Gamma((1 - s + a)/2) / Gamma((s + a)/2) L(1 - s, conj(chi)),

    eps(chi) = tau(chi) / (i^a sqrt(f)) (tau the Gauss sum, a the parity),
    gives the canonical half Re(s) > 1/2, or Re(s) = 1/2 with Im(s) > 0,
    or s = 1/2 with the smaller Conrey label of chi, conj(chi); there
    L(s, chi*) is a generator of kind GR_TOWER_DIRICHLET_L (zeta(s) one of
    kind GR_TOWER_ZETA). So zeta(1/3) is expressed through zeta(2/3).

    The Hurwitz zeta function at a rational parameter p/q (q <= the
    expansion limit) is the character sum

        zeta(s, p/q) = q^s / phi(q) sum_{chi mod q} conj(chi)(p) L(s, chi),

    after a shift into (0, 1], so that its distribution and the relations
    with L-functions hold by construction; Li_s and the Lerch
    transcendent Phi(z, s, a) at roots of unity z reduce to it:
    Phi(z, s, a) = v^(-s) sum_{k < v} z^k zeta(s, (a + k)/v) for z^v = 1.
*/

#include "fmpq.h"
#include "fmpq_poly.h"
#include "arith.h"
#include "ulong_extras.h"
#include "acb.h"
#include "acb_dirichlet.h"
#include "dirichlet.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

#define FLAGS(ctx) gr_tower_lazy_ctx_field_flags(ctx)
#define ALG(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_ALGEBRAIC) && _gr_tower_lazy_outermost(ctx))
#define REAL(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_REAL) && _gr_tower_lazy_outermost(ctx))

#define CHECK_PREC 64

/* characters mod q, roots of unity of order v, in the expansions */
#define DIRICHLET_EXPAND_LIMIT 240

/* shifts of the Hurwitz parameter */
#define DIRICHLET_SHIFT_LIMIT 1000

/* whether x is represented by a rational number (then c) */
static int
_rational(fmpq_t c, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_tower_lazy_rational_recognize(c, x, ctx) == 1;
}

/* res = exp(2 pi i r) */
static int
_unit(gr_ptr res, const fmpq_t r, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    GR_TMP_INIT(t, ctx);
    status = gr_set_fmpq(res, r, ctx);
    status |= gr_pi(t, ctx);
    status |= gr_mul(res, res, t, ctx);
    status |= gr_i(t, ctx);
    status |= gr_mul(res, res, t, ctx);
    status |= gr_mul_2exp_si(res, res, 1, ctx);
    status |= gr_exp(res, res, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* res = chi(n) (conj(chi)(n) if conj) */
static int
_chi(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, ulong n, int conj, gr_ctx_t ctx)
{
    ulong v = dirichlet_chi(G, chi, n % G->q);
    fmpq_t r;
    int status;

    if (v == DIRICHLET_CHI_NULL)
        return gr_zero(res, ctx);
    if (v == 0)
        return gr_one(res, ctx);
    fmpq_init(r);
    fmpq_set_si(r, conj ? -(slong) v : (slong) v, G->expo);
    status = _unit(res, r, ctx);
    fmpq_clear(r);
    return status;
}

/* zeta(s, a) for an integer s <= 0: -B_{1-s}(a) / (1 - s) */
static int
_hurwitz_nonpositive(gr_ptr res, slong s, gr_srcptr a, gr_ctx_t ctx)
{
    fmpq_poly_t B;
    fmpq_t c;
    gr_ptr t;
    slong k, n = 1 - s;
    int status = GR_SUCCESS;

    fmpq_poly_init(B);
    fmpq_init(c);
    GR_TMP_INIT(t, ctx);
    arith_bernoulli_polynomial(B, n);
    status = gr_zero(t, ctx);
    for (k = fmpq_poly_degree(B); k >= 0 && status == GR_SUCCESS; k--)
    {
        status = gr_mul(t, t, a, ctx);
        fmpq_poly_get_coeff_fmpq(c, B, k);
        status |= gr_add_fmpq(t, t, c, ctx);
    }
    status |= gr_div_si(res, t, -n, ctx);
    GR_TMP_CLEAR(t, ctx);
    fmpq_poly_clear(B);
    fmpq_clear(c);
    return status;
}

/* L(n, chi) for the primitive chi = chi_q(k), q >= 3, at an integer n */
static int
_dirichlet_l_prim_integer(gr_ptr res, slong n, ulong q, ulong k, gr_ctx_t ctx)
{
    dirichlet_group_t G;
    dirichlet_char_t chi;
    gr_ptr t, u, a, ss;
    ulong j;
    int status = GR_SUCCESS;

    if (q > DIRICHLET_EXPAND_LIMIT)
        return GR_UNABLE;

    dirichlet_group_init(G, q);
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, k);
    GR_TMP_INIT4(t, u, a, ss, ctx);

    status = gr_zero(res, ctx);
    status |= gr_set_si(ss, n, ctx);
    for (j = 1; j <= q && status == GR_SUCCESS; j++)
    {
        if (n_gcd(j, q) != 1)
            continue;
        status = _chi(t, G, chi, j, 0, ctx);
        status |= gr_set_ui(a, j, ctx);
        status |= gr_div_ui(a, a, q, ctx);
        if (status != GR_SUCCESS)
            break;
        if (n == 1)
            status = gr_digamma(u, a, ctx);
        else if (n >= 2)
            status = gr_hurwitz_zeta(u, ss, a, ctx);
        else
            status = _hurwitz_nonpositive(u, n, a, ctx);
        status |= gr_mul(t, t, u, ctx);
        status |= gr_add(res, res, t, ctx);
    }

    /* L(1, chi) = -(1/q) sum chi(a) psi(a/q); otherwise q^(-n) sum chi(a) zeta(n, a/q) */
    if (status == GR_SUCCESS)
    {
        if (n == 1)
            status = gr_div_si(res, res, -(slong) q, ctx);
        else
        {
            status = gr_set_ui(t, q, ctx);
            status |= gr_pow_si(t, t, -n, ctx);
            status |= gr_mul(res, res, t, ctx);
        }
    }

    GR_TMP_CLEAR4(t, u, a, ss, ctx);
    dirichlet_char_clear(chi);
    dirichlet_group_clear(G);
    return status;
}

/* eps(chi) = tau(chi) / (i^a sqrt(q)) for the primitive chi = chi_q(k) */
static int
_root_number(gr_ptr res, ulong q, ulong k, int * parity, gr_ctx_t ctx)
{
    dirichlet_group_t G;
    dirichlet_char_t chi;
    gr_ptr t, u;
    fmpq_t r;
    ulong n;
    int status = GR_SUCCESS;

    dirichlet_group_init(G, q);
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, k);
    *parity = dirichlet_parity_char(G, chi);
    GR_TMP_INIT2(t, u, ctx);
    fmpq_init(r);

    /* tau(chi) = sum_n chi(n) e(n/q) = sum_n e(v_n/expo + n/q) */
    status = gr_zero(res, ctx);
    for (n = 1; n < q && status == GR_SUCCESS; n++)
    {
        ulong v = dirichlet_chi(G, chi, n);
        fmpq_t r2;
        if (v == DIRICHLET_CHI_NULL)
            continue;
        fmpq_init(r2);
        fmpq_set_si(r, v, G->expo);
        fmpq_set_si(r2, n, q);
        fmpq_add(r, r, r2);
        fmpq_clear(r2);
        status = _unit(t, r, ctx);
        status |= gr_add(res, res, t, ctx);
    }

    /* / (i^a sqrt(q)) */
    status |= gr_set_ui(t, q, ctx);
    status |= gr_sqrt(t, t, ctx);
    if (*parity)
    {
        status |= gr_i(u, ctx);
        status |= gr_mul(t, t, u, ctx);
    }
    status |= gr_div(res, res, t, ctx);

    GR_TMP_CLEAR2(t, u, ctx);
    fmpq_clear(r);
    dirichlet_char_clear(chi);
    dirichlet_group_clear(G);
    return status;
}

/*
    L(s, chi) for the primitive character chi = chi_q(k) (q >= 3, or
    q = 1: zeta(s)), in the normal form.
*/
int
_gr_tower_lazy_dirichlet_l_prim(gr_ptr res, gr_srcptr s_in, ulong q, ulong k, gr_ctx_t ctx)
{
    gr_ptr s, t, u, v, w;
    fmpq_t c;
    ulong kbar;
    int status, sg, canonical, parity = 0;

    if (q == 2 || (q >= 3 && (k >= q || n_gcd(k, q) != 1)))
        return GR_DOMAIN;
    /* (q, k packed in the definition parameter) */
    if (q >= (UWORD(1) << (FLINT_BITS / 2 - 1)))
        return GR_UNABLE;

    fmpq_init(c);
    GR_TMP_INIT5(s, t, u, v, w, ctx);
    status = gr_set(s, s_in, ctx);      /* (res may alias s) */

    if (_rational(c, s, ctx) && fmpz_is_one(fmpq_denref(c)))
    {
        if (fmpz_cmp_si(fmpq_numref(c), -100000) < 0 || fmpz_cmp_si(fmpq_numref(c), 100000) > 0)
            status = GR_UNABLE;
        else if (q == 1)
        {
            /* (the integer itself: gr_zeta of a representation it does not
               recognize as an integer would come back here) */
            status = gr_set_fmpz(t, fmpq_numref(c), ctx);
            if (status == GR_SUCCESS)
                status = gr_zeta(res, t, ctx);
        }
        else
            status = _dirichlet_l_prim_integer(res, fmpz_get_si(fmpq_numref(c)), q, k, ctx);
        goto cleanup;
    }

    /* the canonical half */
    kbar = (q == 1) ? 1 : n_invmod(k, q);
    fmpq_set_si(c, 1, 2);
    status = _gr_tower_lazy_re_cmp(&sg, s, c, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;
    if (sg > 0)
        canonical = 1;
    else if (sg < 0)
        canonical = 0;
    else
    {
        status = _gr_tower_lazy_im_sign(&sg, s, ctx);
        if (status != GR_SUCCESS)
            goto cleanup;
        canonical = (sg > 0) || (sg == 0 && k <= kbar);
    }

    if (canonical)
    {
        if (q == 1)
            status = _gr_tower_lazy_special_gen_locked(res, s, GR_TOWER_ZETA, 0, ctx);
        else
            status = _gr_tower_lazy_special_gen_locked(res, s, GR_TOWER_DIRICHLET_L, GR_TOWER_DIRICHLET_PARAM(q, k), ctx);
        goto cleanup;
    }

    /* eps (q/pi)^(1/2 - s) Gamma((1 - s + a)/2) / Gamma((s + a)/2) L(1 - s, conj(chi));
       (the root number from the Gauss sum: q - 1 roots of unity of order
       lcm(q, order of chi), bounded like the other expansions) */
    if (q > DIRICHLET_EXPAND_LIMIT)
    {
        status = GR_UNABLE;
        goto cleanup;
    }
    if (q == 1)
        status = gr_one(w, ctx);
    else
        status = _root_number(w, q, k, &parity, ctx);

    /* (q/pi)^(1/2 - s) */
    status |= gr_set_ui(t, q, ctx);
    status |= gr_pi(u, ctx);
    status |= gr_div(t, t, u, ctx);
    fmpq_set_si(c, 1, 2);
    status |= gr_set_fmpq(u, c, ctx);
    status |= gr_sub(u, u, s, ctx);
    status |= gr_pow(t, t, u, ctx);
    status |= gr_mul(w, w, t, ctx);

    /* Gamma((1 - s + a)/2) / Gamma((s + a)/2) */
    status |= gr_neg(t, s, ctx);
    status |= gr_add_ui(t, t, 1 + parity, ctx);
    status |= gr_div_ui(t, t, 2, ctx);
    if (status == GR_SUCCESS)
        status = gr_gamma(t, t, ctx);
    status |= gr_add_ui(u, s, parity, ctx);
    status |= gr_div_ui(u, u, 2, ctx);
    if (status == GR_SUCCESS)
        status = gr_gamma(u, u, ctx);
    status |= gr_div(t, t, u, ctx);
    status |= gr_mul(w, w, t, ctx);

    /* L(1 - s, conj(chi)) */
    status |= gr_neg(t, s, ctx);
    status |= gr_add_ui(t, t, 1, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_dirichlet_l_prim(v, t, q, kbar, ctx);
    if (status == GR_SUCCESS)
        status = gr_mul(res, w, v, ctx);

cleanup:
    GR_TMP_CLEAR5(s, t, u, v, w, ctx);
    fmpq_clear(c);
    return status;
}

/* L(s, chi) for any character chi mod q (unlocked) */
static int
_dirichlet_l(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr s_in, gr_ctx_t ctx)
{
    ulong f = dirichlet_conductor_char(G, chi), k = 1, q = G->q;
    gr_ptr s, t, u;
    n_factor_t fac;
    slong i;
    int status = GR_SUCCESS;
    dirichlet_group_t H;
    dirichlet_char_t y;

    GR_TMP_INIT3(s, t, u, ctx);
    status = gr_set(s, s_in, ctx);

    if (f > 1)
    {
        dirichlet_group_init(H, f);
        dirichlet_char_init(y, H);
        dirichlet_char_lower(y, H, chi, G);
        k = dirichlet_char_exp(H, y);
    }

    /* L(s, chi*) */
    if (status == GR_SUCCESS)
    {
        if (f == 1)
        {
            /* (the pole at s = 1) */
            truth_t one = gr_is_one(s, ctx);
            if (one == T_TRUE)
                status = GR_DOMAIN;
            else if (one == T_UNKNOWN)
                status = GR_UNABLE;
            else
                status = _gr_tower_lazy_dirichlet_l_prim(res, s, 1, 1, ctx);
        }
        else
            status = _gr_tower_lazy_dirichlet_l_prim(res, s, f, k, ctx);
    }

    /* the Euler factors at p | q, p not dividing f: (1 - chi*(p) p^(-s)) */
    n_factor_init(&fac);
    n_factor(&fac, q, 1);
    for (i = 0; i < fac.num && status == GR_SUCCESS; i++)
    {
        ulong p = fac.p[i];
        if (f % p == 0)
            continue;
        status = gr_set_ui(t, p, ctx);
        status |= gr_neg(u, s, ctx);
        status |= gr_pow(t, t, u, ctx);
        if (f > 1)
        {
            status |= _chi(u, H, y, p, 0, ctx);
            status |= gr_mul(t, t, u, ctx);
        }
        status |= gr_sub_ui(t, t, 1, ctx);
        status |= gr_neg(t, t, ctx);
        status |= gr_mul(res, res, t, ctx);
    }

    if (f > 1)
    {
        dirichlet_char_clear(y);
        dirichlet_group_clear(H);
    }
    GR_TMP_CLEAR3(s, t, u, ctx);
    return status;
}

int
gr_tower_lazy_dirichlet_l(gr_tower_lazy_elem_t res, const dirichlet_group_t G, const dirichlet_char_t chi,
    const gr_tower_lazy_elem_t s, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_dirichlet_l(res, G, chi, s, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/*
    zeta(s, p/q), 0 < p <= q, gcd(p, q) = 1, s not an integer:
    q^s / phi(q) sum_{chi mod q} conj(chi)(p) L(s, chi)
*/
int
_gr_tower_lazy_hurwitz_rational(gr_ptr res, gr_srcptr s, slong p, slong q, gr_ctx_t ctx)
{
    dirichlet_group_t G;
    dirichlet_char_t chi;
    gr_ptr t, u;
    int status = GR_SUCCESS;

    if (q == 1)
        return _gr_tower_lazy_dirichlet_l_prim(res, s, 1, 1, ctx);
    if (q > DIRICHLET_EXPAND_LIMIT)
        return GR_UNABLE;

    dirichlet_group_init(G, q);
    dirichlet_char_init(chi, G);
    GR_TMP_INIT2(t, u, ctx);

    status = gr_zero(res, ctx);
    dirichlet_char_one(chi, G);
    do
    {
        status = _chi(t, G, chi, p, 1, ctx);
        if (status == GR_SUCCESS)
            status = _dirichlet_l(u, G, chi, s, ctx);
        status |= gr_mul(t, t, u, ctx);
        status |= gr_add(res, res, t, ctx);
    }
    while (status == GR_SUCCESS && dirichlet_char_next(chi, G) >= 0);

    if (status == GR_SUCCESS)
    {
        status = gr_set_ui(t, q, ctx);
        status |= gr_pow(t, t, s, ctx);
        status |= gr_mul(res, res, t, ctx);
        status |= gr_div_ui(res, res, G->phi_q, ctx);
    }

    GR_TMP_CLEAR2(t, u, ctx);
    dirichlet_char_clear(chi);
    dirichlet_group_clear(G);
    return status;
}

/*
    zeta(s, a) for irrational a, s not an integer: shifted to
    0 < Re(a0) <= 1 (zeta(s, a) = a^(-s) + zeta(s, a + 1), principal
    powers), where it is a generator of kind GR_TOWER_HURWITZ_ZETA
*/
static int
_hurwitz_irrational(gr_ptr res, gr_srcptr s, gr_srcptr a, gr_ctx_t ctx)
{
    acb_t z;
    fmpz_t n;
    fmpq_t c;
    gr_ptr a0, t, e, args;
    slong nn, i;
    int status, sg;

    acb_init(z);
    fmpz_init(n);
    fmpq_init(c);
    GR_TMP_INIT3(a0, t, e, ctx);

    /* n = ceil(Re(a)) - 1 numerically, then 0 < Re(a - n) <= 1 exactly */
    status = gr_tower_lazy_get_acb(z, a, CHECK_PREC, ctx);
    if (status == GR_SUCCESS && arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 20) < 0)
    {
        arf_get_fmpz(n, arb_midref(acb_realref(z)), ARF_RND_CEIL);
        fmpz_sub_ui(n, n, 1);
        nn = fmpz_get_si(n);
        status = gr_sub_si(a0, a, nn, ctx);
        fmpq_zero(c);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_re_cmp(&sg, a0, c, ctx);
        if (status == GR_SUCCESS && sg <= 0)
        {
            nn--;
            status = gr_add_si(a0, a0, 1, ctx);
        }
        fmpq_one(c);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_re_cmp(&sg, a0, c, ctx);
        if (status == GR_SUCCESS && sg > 0)
        {
            nn++;
            status = gr_sub_si(a0, a0, 1, ctx);
        }
        if (status == GR_SUCCESS && FLINT_ABS(nn) > DIRICHLET_SHIFT_LIMIT)
            status = GR_UNABLE;
    }
    else if (status == GR_SUCCESS)
        status = GR_UNABLE;

    if (status == GR_SUCCESS)
    {
        args = gr_heap_init_vec(2, ctx);
        status = gr_set(args, s, ctx);
        status |= gr_set(GR_ENTRY(args, 1, ctx->sizeof_elem), a0, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_special_gen_multi_locked(t, args, 2, GR_TOWER_HURWITZ_ZETA, 0, ctx);
        gr_heap_clear_vec(args, 2, ctx);

        /* a = a0 + nn: zeta(s, a) = zeta(s, a0) - sum_{k<nn} (a0 + k)^(-s),
           or + sum_{k=1}^{-nn} (a0 - k)^(-s) */
        status |= gr_neg(e, s, ctx);
        for (i = 0; i < FLINT_ABS(nn) && status == GR_SUCCESS; i++)
        {
            gr_ptr u;
            GR_TMP_INIT(u, ctx);
            status = (nn > 0) ? gr_add_si(u, a0, i, ctx) : gr_sub_si(u, a0, i + 1, ctx);
            status |= gr_pow(u, u, e, ctx);
            status |= (nn > 0) ? gr_sub(t, t, u, ctx) : gr_add(t, t, u, ctx);
            GR_TMP_CLEAR(u, ctx);
        }
        if (status == GR_SUCCESS)
            status = gr_set(res, t, ctx);
    }

    acb_clear(z);
    fmpz_clear(n);
    fmpq_clear(c);
    GR_TMP_CLEAR3(a0, t, e, ctx);
    return status;
}

/*
    zeta(s, a) for s not an integer >= 2 (those are polygamma values):
    -B_{1-s}(a)/(1-s) at integers s <= 0; the character sums at
    rational a (shifted into (0, 1]: zeta(s, a) = zeta(s, a + n) +
    sum_{k<n} (a + k)^(-s)); a generator at irrational a.
*/
int
_gr_tower_lazy_hurwitz_general(gr_ptr res, gr_srcptr s, gr_srcptr a, gr_ctx_t ctx)
{
    fmpq_t c, b, b0;
    fmpz_t n;
    gr_ptr t, u, e;
    slong i, nn;
    int status = GR_SUCCESS;

    fmpq_init(c);
    fmpq_init(b);
    fmpq_init(b0);
    fmpz_init(n);

    if (_rational(c, s, ctx) && fmpz_is_one(fmpq_denref(c)))
    {
        if (fmpz_is_one(fmpq_numref(c)))
            status = GR_DOMAIN;
        else if (fmpz_cmp_ui(fmpq_numref(c), 2) >= 0)
        {
            /* (polygamma; with the integer itself: gr_hurwitz_zeta would
               come back here when it cannot) */
            if (fmpz_fits_si(fmpq_numref(c)))
                status = _gr_tower_lazy_hurwitz_zeta_int(res, fmpz_get_si(fmpq_numref(c)), a, ctx);
            else
                status = GR_UNABLE;
        }
        else if (fmpz_sgn(fmpq_numref(c)) <= 0 && fmpz_cmp_si(fmpq_numref(c), -1000) >= 0)
            status = _hurwitz_nonpositive(res, fmpz_get_si(fmpq_numref(c)), a, ctx);
        else
            status = GR_UNABLE;
        goto cleanup;
    }

    if (!_rational(b, a, ctx))
    {
        status = _hurwitz_irrational(res, s, a, ctx);
        goto cleanup;
    }

    /* a = b0 + n with 0 < b0 <= 1 */
    fmpz_cdiv_q(n, fmpq_numref(b), fmpq_denref(b));
    fmpz_sub_ui(n, n, 1);
    fmpq_sub_fmpz(b0, b, n);
    if (fmpz_cmp_si(n, -DIRICHLET_SHIFT_LIMIT) < 0 || fmpz_cmp_si(n, DIRICHLET_SHIFT_LIMIT) > 0 || !fmpz_fits_si(fmpq_denref(b0)) ||
        fmpz_cmp_ui(fmpq_denref(b0), DIRICHLET_EXPAND_LIMIT) > 0)
    {
        status = GR_UNABLE;
        goto cleanup;
    }
    /* (a nonpositive integer: the term n + a = 0) */
    if (fmpz_sgn(n) < 0 && fmpz_is_one(fmpq_denref(b)))
    {
        status = GR_DOMAIN;
        goto cleanup;
    }

    GR_TMP_INIT3(t, u, e, ctx);
    status = _gr_tower_lazy_hurwitz_rational(res, s, fmpz_get_si(fmpq_numref(b0)), fmpz_get_si(fmpq_denref(b0)), ctx);

    /* zeta(s, b0 + n) = zeta(s, b0) - sum_{k<n} (b0 + k)^(-s);
       zeta(s, b0 - n) = zeta(s, b0) + sum_{k=1}^{n} (b0 - k)^(-s) */
    nn = fmpz_get_si(n);
    status |= gr_neg(e, s, ctx);
    for (i = 0; i < FLINT_ABS(nn) && status == GR_SUCCESS; i++)
    {
        fmpq_set(c, b0);
        if (nn > 0)
            fmpq_add_si(c, c, i);
        else
            fmpq_sub_si(c, c, i + 1);
        status = gr_set_fmpq(t, c, ctx);
        status |= gr_pow(t, t, e, ctx);
        if (nn > 0)
            status |= gr_sub(res, res, t, ctx);
        else
            status |= gr_add(res, res, t, ctx);
    }
    GR_TMP_CLEAR3(t, u, e, ctx);

cleanup:
    fmpq_clear(c);
    fmpq_clear(b);
    fmpq_clear(b0);
    fmpz_clear(n);
    return status;
}

/*
    Phi(z, s, a) = sum_n z^n (n + a)^(-s): zeta(s, a) at z = 1, the
    expansion at the roots of unity z = e(u/v) of order v <= the limit,
    Li_s(z)/z at a = 1 (s an integer); GR_UNABLE otherwise.
*/
static int
_lerch_phi(gr_ptr res, gr_srcptr z, gr_srcptr s, gr_srcptr a, gr_ctx_t ctx)
{
    fmpq_t r, c;
    truth_t one;
    int status = GR_SUCCESS;

    fmpq_init(r);
    fmpq_init(c);

    one = gr_is_one(z, ctx);
    if (one == T_TRUE)
    {
        status = gr_hurwitz_zeta(res, s, a, ctx);
    }
    else if (_gr_tower_lazy_root_of_unity_angle_locked(r, z, ctx) &&
             fmpz_cmp_ui(fmpq_denref(r), DIRICHLET_EXPAND_LIMIT) <= 0)
    {
        /* v^(-s) sum_{k<v} z^k zeta(s, (a + k)/v), z = e(r) */
        slong v = fmpz_get_si(fmpq_denref(r)), k;
        gr_ptr t, u, w;
        /* (s = 1: the poles cancel, sum z^k = 0, and -psi stands for zeta(1, .)) */
        int s1 = (gr_is_one(s, ctx) == T_TRUE);
        GR_TMP_INIT3(t, u, w, ctx);
        status = gr_zero(res, ctx);
        for (k = 0; k < v && status == GR_SUCCESS; k++)
        {
            fmpq_mul_si(c, r, k);
            status = _unit(t, c, ctx);
            status |= gr_add_si(u, a, k, ctx);
            status |= gr_div_si(u, u, v, ctx);
            if (status == GR_SUCCESS && s1)
            {
                status = gr_digamma(w, u, ctx);
                status |= gr_neg(w, w, ctx);
            }
            else if (status == GR_SUCCESS)
                status = gr_hurwitz_zeta(w, s, u, ctx);
            status |= gr_mul(t, t, w, ctx);
            status |= gr_add(res, res, t, ctx);
        }
        status |= gr_set_si(t, v, ctx);
        status |= gr_neg(u, s, ctx);
        status |= gr_pow(t, t, u, ctx);
        status |= gr_mul(res, res, t, ctx);
        GR_TMP_CLEAR3(t, u, w, ctx);
    }
    else if (gr_is_one(a, ctx) == T_TRUE && _rational(c, s, ctx) && fmpz_is_one(fmpq_denref(c)))
    {
        /* Li_s(z) / z */
        status = gr_polylog(res, s, z, ctx);
        status |= gr_div(res, res, z, ctx);
    }
    else
        status = GR_UNABLE;

    fmpq_clear(r);
    fmpq_clear(c);
    return status;
}

int
gr_tower_lazy_lerch_phi(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t z, const gr_tower_lazy_elem_t s,
    const gr_tower_lazy_elem_t a, gr_ctx_t ctx)
{
    gr_ptr zz, ss, aa;
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    /* (res may alias the arguments) */
    GR_TMP_INIT3(zz, ss, aa, ctx);
    status = gr_set(zz, z, ctx);
    status |= gr_set(ss, s, ctx);
    status |= gr_set(aa, a, ctx);
    if (status == GR_SUCCESS)
        status = _lerch_phi(res, zz, ss, aa, ctx);
    GR_TMP_CLEAR3(zz, ss, aa, ctx);
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

POP_OPTIONS
