/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Rigorous evaluation and bounds of the incomplete gamma functions by
    integration along steepest-descent contours (sd_plan.c,
    sd_integrate.c).  In the logarithmic variable w = log t,

        Gamma(s, z) = int_{Log z}^{V+(0)} exp(phi(w)) dw,
        gamma(s, z) = int_{V}^{Log z} exp(phi(w)) dw,

    with phi(w) = -e^w + s w, that is, alpha = -1, beta = 0, nu = -s in
    the notation of sd.h, where V+(0) is the valley Re w -> +inf, Im w
    near 0, and V the valley to the left, w -> inf in a direction along
    which Re(s w) -> -inf: for Re s > 0 this is the defining integral,
    and otherwise (in a direction turned away from the real axis, when
    Im s != 0) its analytic continuation in s.  The principal branch
    Log z gives the principal branches of the functions; a ball z
    crossing the cut (-inf, 0) is split along it.

    The path of steepest descent from Log z runs into one of the two
    valleys (or another one, or a saddle).  The function whose contour
    is found by the planner (directly, or through the saddles) is
    integrated; the other one is obtained from Gamma(s, z) + gamma(s, z)
    = Gamma(s), which is then free from cancellation.  This matters
    when |Im s| is large: the valley V is then separated from V+(0) by
    infinitely many saddles (with contributions in geometric
    progression, summing to Gamma(s)).  If the descent flow leads
    nowhere (s near 0), the upper function is integrated along a fixed
    polygon.
*/

#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"
#include "acb_hypgeom/sd.h"

#define BOUND_DROP 24.0
#define BOUND_TAILTOL 16.0
#define QUAD_MAXDROP 16.0

static int
crosses_cut(const acb_t z)
{
    return arb_contains_zero(acb_imagref(z)) && !arb_is_exact(acb_imagref(z))
        && !arb_is_positive(acb_realref(z));
}

/* The problem with a contour for the upper function or for minus the
   lower function (as flagged by upper) */
typedef struct
{
    sd_problem P;
    sd_work W;
    sd_contour C;
    acb_t alpha, nu, start;
    int upper;
}
gamma_ctx;

static int
gamma_ctx_init(gamma_ctx * G, const acb_t s, const acb_t z, int side, int upper,
    double drop, double tailtol, double maxdrop, slong wp)
{
    sd_end a = sd_end_make(SD_START, 0), vp = sd_end_make(SD_PLUS, 0), vl = sd_end_make(SD_LINEAR, 0);
    int ok;

    acb_init(G->alpha);
    acb_init(G->nu);
    acb_init(G->start);

    acb_set_si(G->alpha, -1);
    acb_neg(G->nu, s);
    /* Log z, continuous on a ball touching the cut from the given side
       (or on either side, if side = 0): Log z = Log(-z) +- pi i for
       +-Im z >= 0 */
    if (arb_is_negative(acb_realref(z)))
    {
        arb_t pi;
        arb_init(pi);
        arb_const_pi(pi, wp);
        acb_neg(G->start, z);
        acb_log(G->start, G->start, wp);
        if (side > 0 || (side == 0 && arb_is_nonnegative(acb_imagref(z))))
            arb_add(acb_imagref(G->start), acb_imagref(G->start), pi, wp);
        else
            arb_sub(acb_imagref(G->start), acb_imagref(G->start), pi, wp);
        arb_clear(pi);
    }
    else
        acb_log(G->start, z, wp);

    if (_sd_problem_init(&G->P, G->alpha, 0, G->nu, G->start, drop, tailtol, maxdrop))
    {
        _sd_work_init(&G->W, &G->P);

        G->upper = upper;
        ok = _sd_plan(&G->C, &G->W, &G->P, a, upper ? vp : vl, 0);
        if (!ok)
        {
            G->upper = !upper;
            ok = _sd_plan(&G->C, &G->W, &G->P, a, upper ? vl : vp, 0);
        }
        if (!ok)
        {
            G->upper = 1;
            ok = _sd_plan(&G->C, &G->W, &G->P, a, vp, 1);
        }
        if (ok)
            return 1;

        _sd_work_clear(&G->W);
    }

    acb_clear(G->alpha);
    acb_clear(G->nu);
    acb_clear(G->start);
    return 0;
}

static void
gamma_ctx_clear(gamma_ctx * G)
{
    _sd_contour_clear(&G->C);
    _sd_work_clear(&G->W);
    acb_clear(G->alpha);
    acb_clear(G->nu);
    acb_clear(G->start);
}

static int
gamma_nosplit(acb_t res, acb_t FLINT_UNUSED(res2), const acb_t s, const acb_t z, int side, slong prec, void * param)
{
    gamma_ctx G;
    double tol = prec * SD_LN2 + 12.0;
    int ok, upper = *(int *) param;

    if (!gamma_ctx_init(&G, s, z, side, upper, tol, tol, QUAD_MAXDROP, prec + 30))
        return 0;

    ok = _sd_integrate(res, &G.C, &G.P, G.alpha, 0, G.nu, G.start, prec);
    if (!G.upper)
        acb_neg(res, res);

    if (ok && G.upper != upper)
    {
        acb_t g;
        acb_init(g);
        acb_gamma(g, s, prec + 20);
        acb_sub(res, g, res, prec + 20);
        acb_clear(g);
    }

    gamma_ctx_clear(&G);
    return ok;
}

int
_acb_hypgeom_gamma_incomplete_integration(acb_t res, const acb_t s, const acb_t z, int upper, slong prec)
{
    acb_t t;
    int ok;

    if (!acb_is_finite(s) || !acb_is_finite(z) || acb_contains_zero(z))
        return 0;

    acb_init(t);

    if (crosses_cut(z))
        ok = _acb_hypgeom_split_at_cut(t, NULL, s, z, prec, gamma_nosplit, &upper);
    else
        ok = gamma_nosplit(t, NULL, s, z, 0, prec, &upper);

    ok = ok && acb_is_finite(t);
    if (ok)
        acb_set_round(res, t, prec);

    acb_clear(t);
    return ok;
}

/* Bounds for sup |F| and sup |dF/ds| over the input balls from the L1
   norms along the contour (dF/ds has the integrand w exp(phi)), with
   bounds for Gamma(s) and Gamma'(s) = Gamma(s) psi(s) added for the
   complementary function. */
int
_acb_hypgeom_gamma_incomplete_integration_bounds(mag_t M, mag_t Ds,
    const acb_t s, const acb_t z, int upper)
{
    gamma_ctx G;
    int ok = 1;

    mag_inf(M);
    mag_inf(Ds);

    if (!acb_is_finite(s) || !acb_is_finite(z) || acb_contains_zero(z) || crosses_cut(z)
        || !gamma_ctx_init(&G, s, z, 0, upper, BOUND_DROP, BOUND_TAILTOL, 1.0, 64))
        return 0;

    _sd_mag_set_exp_d(M, G.C.logL[0]);
    _sd_mag_set_exp_d(Ds, G.C.logL[2]);

    if (G.upper != upper)
    {
        acb_t g, d;
        mag_t t;

        acb_init(g);
        acb_init(d);
        mag_init(t);
        acb_gamma(g, s, 32);
        acb_digamma(d, s, 32);
        acb_mul(d, d, g, 32);
        acb_get_mag(t, g);
        mag_add(M, M, t);
        acb_get_mag(t, d);
        mag_add(Ds, Ds, t);
        ok = mag_is_finite(M) && mag_is_finite(Ds);
        acb_clear(g);
        acb_clear(d);
        mag_clear(t);
    }

    gamma_ctx_clear(&G);
    return ok;
}
