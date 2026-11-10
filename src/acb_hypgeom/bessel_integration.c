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
    Rigorous evaluation and bounds of the Bessel functions from the
    Schlafli integrals along steepest-descent contours (sd_plan.c,
    sd_integrate.c).  With phi(w) = zeta sinh w - nu w, that is,
    alpha = -beta = zeta / 2 in the notation of sd.h, and the valleys
    V+(k): Re w -> +inf, Im w near (2k+1) pi - arg zeta and
    V-(k): Re w -> -inf, Im w near arg zeta + 2 pi k
    (DLMF 10.9.17-18, 10.27.6, 10.32.8),

        J_nu(z)  = (1/2 pi i) int_{V+(-1)}^{V+(0)},         zeta = z
        H1_nu(z) = (1/pi i)   int_{V-(0)}^{V+(0)},          zeta = z
        H2_nu(z) = -(1/pi i)  int_{V-(0)}^{V+(-1)},         zeta = z
        Y_nu(z)  = (H1_nu(z) - H2_nu(z)) / (2 i)
        I_nu(z)  = e^{-+ i pi nu/2} J_nu(+- i z)
        K_nu(z)  = (1/2) e^{- i pi nu/2} int_{V-(0)}^{V+(-1)},  zeta = -i z
                 = (1/2) e^{+ i pi nu/2} int_{V-(0)}^{V+(0)},   zeta = i z

    Since the valleys follow arg zeta in (-pi, pi], this gives the
    principal branches, continuous from above on the cut (-inf, 0).  The
    integrals are analytic in zeta, so for a ball z the enclosures and
    bounds are valid for the continuation from the midpoint: a ball
    crossing the cut is split along it.
*/

#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"
#include "acb_hypgeom/sd.h"

#define KIND_J ACB_HYPGEOM_KIND_BESSEL_J
#define KIND_Y ACB_HYPGEOM_KIND_BESSEL_Y
#define KIND_JY ACB_HYPGEOM_KIND_BESSEL_JY
#define KIND_H1 ACB_HYPGEOM_KIND_HANKEL_1
#define KIND_H2 ACB_HYPGEOM_KIND_HANKEL_2
#define KIND_I ACB_HYPGEOM_KIND_BESSEL_I
#define KIND_K ACB_HYPGEOM_KIND_BESSEL_K

/* the bounds: trace until exp(Re phi - top) < e^-24 */
#define BOUND_DROP 24.0
#define BOUND_TAILTOL 16.0
/* the largest drop of Re phi per planner step for the quadrature (the
   chords are chosen from the planned polygon) */
#define QUAD_MAXDROP 16.0

/* the rotation zeta = z, i z, -i z (rot = 0, 1, -1) for each function */
static int
bessel_rot(int kind, const acb_t z)
{
    double arg = atan2(arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR),
                       arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR));

    /* I_nu(z) = e^{-i pi nu/2} J_nu(i z) (arg z <= pi/2),
               = e^{+i pi nu/2} J_nu(-i z) (arg z >= -pi/2) */
    if (kind == KIND_I)
        return (arg <= SD_PI / 4) ? 1 : -1;
    /* the first form of K for arg z >= 0 */
    if (kind == KIND_K)
        return (arg >= 0) ? -1 : 1;
    return 0;
}

/* alpha = zeta / 2 */
static void
bessel_alpha(acb_t alpha, const acb_t z, int rot)
{
    if (rot == 0) acb_set(alpha, z);
    else if (rot == 1) acb_mul_onei(alpha, z);
    else acb_div_onei(alpha, z);
    acb_mul_2exp_si(alpha, alpha, -1);
}

/* whether the ball z crosses the branch cut */
static int
crosses_cut(const acb_t z)
{
    return arb_contains_zero(acb_imagref(z)) && !arb_is_exact(acb_imagref(z))
        && !arb_is_positive(acb_realref(z));
}

typedef struct
{
    sd_problem P;
    sd_work W;
    acb_t alpha;
    acb_srcptr nu;
    slong prec;
}
bessel_ctx;

static int
bessel_ctx_init(bessel_ctx * B, const acb_t nu, const acb_t z, int rot, slong prec)
{
    double tol = prec * SD_LN2 + 12.0;

    acb_init(B->alpha);
    bessel_alpha(B->alpha, z, rot);
    B->nu = nu;
    B->prec = prec;
    if (!_sd_problem_init(&B->P, B->alpha, 1, nu, NULL, tol, tol, QUAD_MAXDROP))
    {
        acb_clear(B->alpha);
        return 0;
    }
    _sd_work_init(&B->W, &B->P);
    return 1;
}

static void
bessel_ctx_clear(bessel_ctx * B)
{
    _sd_work_clear(&B->W);
    acb_clear(B->alpha);
}

/* the integral from the valley a to the valley b */
static int
bessel_integral(acb_t res, bessel_ctx * B, sd_end a, sd_end b)
{
    sd_contour C;
    int ok;

    if (!_sd_plan(&C, &B->W, &B->P, a, b, 1))
        return 0;
    /* by the saddle-point expansions if they apply and are expected to
       cost less than the quadrature (with a margin for the accuracy of
       the estimates), and otherwise by quadrature */
    ok = _sd_expand(res, &C, &B->P, B->alpha, 1, B->nu, B->prec, 0.9)
        || _sd_integrate(res, &C, &B->P, B->alpha, 1, B->nu, NULL, B->prec);
    _sd_contour_clear(&C);
    return ok;
}

/* J, Y, H1 or H2, or J and Y together */
static int
bessel_jyh(acb_t res, acb_t res2, int kind, const acb_t nu, const acb_t z, slong prec)
{
    bessel_ctx B;
    acb_t I1, I2, t;
    sd_end vm = sd_end_make(SD_MINUS, 0), vp0 = sd_end_make(SD_PLUS, 0), vp1 = sd_end_make(SD_PLUS, -1);
    slong wp = prec + 20;
    int ok = 1, need1, need2, real, direct_j, want_j, want_y;

    if (!bessel_ctx_init(&B, nu, z, 0, prec))
        return 0;

    acb_init(I1);
    acb_init(I2);
    acb_init(t);

    want_j = (kind == KIND_J || kind == KIND_JY);
    want_y = (kind == KIND_Y || kind == KIND_JY);

    /* For real nu and z > 0, the Hankel integrals int_1 (to V+(0)) and
       int_2 (to V+(-1)) are complex conjugates.  J = (int_1 - int_2) /
       (2 pi i) from these cancels where J is exponentially smaller than
       Y (|z| < |nu|); it is then integrated directly. */
    real = acb_is_real(nu) && acb_is_real(z) && arb_is_positive(acb_realref(z));
    direct_j = want_j && !((real || want_y) && 2.0 * B.P.aal >= B.P.anu);
    need2 = !real && (kind == KIND_H2 || want_y || (want_j && !direct_j));
    need1 = (kind == KIND_H1 || want_y || (want_j && !direct_j) || (real && kind == KIND_H2));

    if (need1)
        ok = bessel_integral(I1, &B, vm, vp0);
    if (ok && need2)
        ok = bessel_integral(I2, &B, vm, vp1);
    if (ok && real)
        acb_conj(I2, I1);
    if (ok && direct_j)
        ok = bessel_integral(t, &B, vp1, vp0);

    if (ok)
    {
        arb_t pi;
        arb_init(pi);
        arb_const_pi(pi, wp);

        if (kind == KIND_H1)
        {
            acb_div_onei(res, I1);
            acb_div_arb(res, res, pi, wp);
        }
        else if (kind == KIND_H2)
        {
            acb_mul_onei(res, I2);
            acb_div_arb(res, res, pi, wp);
        }

        if (want_j)
        {
            /* (int_1 - int_2) / (2 pi i) */
            if (!direct_j)
                acb_sub(t, I1, I2, wp);
            acb_div_onei(t, t);
            acb_div_arb(t, t, pi, wp);
            acb_mul_2exp_si(res, t, -1);
        }

        if (want_y)
        {
            /* -(int_1 + int_2) / (2 pi) */
            acb_ptr y = (kind == KIND_JY) ? res2 : res;
            acb_add(t, I1, I2, wp);
            acb_div_arb(t, t, pi, wp);
            acb_mul_2exp_si(t, t, -1);
            acb_neg(y, t);
        }

        arb_clear(pi);
    }

    bessel_ctx_clear(&B);
    acb_clear(I1);
    acb_clear(I2);
    acb_clear(t);
    return ok;
}

/* the number of bits lost in exp(pi i nu) through the imaginary part
   of nu */
static slong
bessel_order_bits(const acb_t nu)
{
    slong e;

    if (arf_is_zero(arb_midref(acb_imagref(nu))) || !arf_is_finite(arb_midref(acb_imagref(nu))))
        return 0;
    e = arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(nu)));
    return FLINT_MIN(FLINT_MAX(e, 0), 128);
}

/* I or K, optionally scaled by exp(-+z) */
static int
bessel_ik(acb_t res, int kind, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    bessel_ctx B;
    acb_t t, u;
    slong wp = prec + 20;
    int ok, rot = bessel_rot(kind, z);

    if (!bessel_ctx_init(&B, nu, z, rot, prec))
        return 0;

    acb_init(t);
    acb_init(u);

    if (kind == KIND_K)
        ok = bessel_integral(t, &B, sd_end_make(SD_MINUS, 0), sd_end_make(SD_PLUS, (rot == -1) ? -1 : 0));
    else
        ok = bessel_integral(t, &B, sd_end_make(SD_PLUS, -1), sd_end_make(SD_PLUS, 0));

    if (ok)
    {
        /* the factor 1/2 or 1/(2 pi i), times e^{-+ i pi nu / 2} (and
           e^{-+ z} if scaled) */
        acb_const_pi(u, wp);
        if (kind == KIND_I)
        {
            acb_div(t, t, u, wp);
            acb_div_onei(t, t);
        }
        acb_mul_2exp_si(t, t, -1);

        /* (with the real part of the argument reduced exactly and
           extra precision for its imaginary part: no loss of accuracy
           for large orders) */
        acb_mul_2exp_si(u, nu, -1);
        if ((kind == KIND_I) == (rot == 1))
            acb_neg(u, u);
        acb_exp_pi_i(u, u, wp + bessel_order_bits(nu));
        acb_mul(t, t, u, wp);
        if (scaled)
        {
            if (kind == KIND_I)
                acb_neg(u, z);
            else
                acb_set(u, z);
            acb_exp(u, u, wp);
            acb_mul(t, t, u, wp);
        }
        acb_set(res, t);
    }

    bessel_ctx_clear(&B);
    acb_clear(t);
    acb_clear(u);
    return ok;
}

static int
bessel_nosplit(acb_t res, acb_t res2, int kind, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    if (kind == KIND_I || kind == KIND_K)
        return bessel_ik(res, kind, nu, z, scaled, prec);
    return bessel_jyh(res, res2, kind, nu, z, prec);
}

/* res = the union of f(nu, z+) and f(nu, z-), where z+ and z- are the
   parts Im >= 0 and Im <= 0 of a ball z crossing the cut (-inf, 0); f
   is told which side it is called for (the balls z+ and z- may reach
   slightly across the real axis due to rounding) */
int
_acb_hypgeom_split_at_cut(acb_t res, acb_t res2, const acb_t nu, const acb_t z, slong prec,
    int (*f)(acb_t, acb_t, const acb_t, const acb_t, int, slong, void *), void * param)
{
    arf_t a, b, zero;
    acb_t zz, t, t2;
    int ok;

    arf_init(a);
    arf_init(b);
    arf_init(zero);
    acb_init(zz);
    acb_init(t);
    acb_init(t2);

    /* (rounded outwards) */
    arb_get_interval_arf(a, b, acb_imagref(z), 2 * MAG_BITS);
    arb_set(acb_realref(zz), acb_realref(z));
    arb_set_interval_arf(acb_imagref(zz), zero, b, 2 * MAG_BITS);
    ok = f(res, res2, nu, zz, 1, prec, param);

    /* (points with Im = 0 belong to the upper part) */
    if (ok && arf_sgn(a) < 0)
    {
        arb_set_interval_arf(acb_imagref(zz), a, zero, 2 * MAG_BITS);
        ok = f(t, t2, nu, zz, -1, prec, param);
        if (ok)
        {
            acb_union(res, res, t, prec);
            if (res2 != NULL)
                acb_union(res2, res2, t2, prec);
        }
    }

    arf_clear(a);
    arf_clear(b);
    arf_clear(zero);
    acb_clear(zz);
    acb_clear(t);
    acb_clear(t2);
    return ok;
}

typedef struct { int kind; int scaled; } bessel_param;

static int
bessel_func(acb_t res, acb_t res2, const acb_t nu, const acb_t z, int FLINT_UNUSED(side), slong prec, void * param)
{
    return bessel_nosplit(res, res2, ((bessel_param *) param)->kind, nu, z,
        ((bessel_param *) param)->scaled, prec);
}

int
_acb_hypgeom_bessel_integration(acb_t res, acb_t res2, int kind, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    acb_t t, t2;
    int ok;

    if (!acb_is_finite(nu) || !acb_is_finite(z) || acb_contains_zero(z))
        return 0;

    acb_init(t);
    acb_init(t2);

    /* K_{-nu} = K_nu, H1_{-nu} = e^{i pi nu} H1_nu, H2_{-nu} = e^{-i pi nu}
       H2_nu, and for integer n, J_{-n} = (-1)^n J_n, Y_{-n} = (-1)^n Y_n,
       I_{-n} = I_n (where the integrals for J and I would cancel) */
    if (arf_sgn(arb_midref(acb_realref(nu))) < 0
        && (kind == KIND_K || kind == KIND_H1 || kind == KIND_H2 || acb_is_int(nu)))
    {
        int odd = acb_is_int(nu) && !arf_is_int_2exp_si(arb_midref(acb_realref(nu)), 1);

        acb_neg(t, nu);
        ok = _acb_hypgeom_bessel_integration(res, res2, kind, t, z, scaled, prec);

        if (ok && (kind == KIND_H1 || kind == KIND_H2))
        {
            acb_exp_pi_i(t, t, prec + 10 + bessel_order_bits(nu));
            if (kind == KIND_H1)
                acb_mul(res, res, t, prec);
            else
                acb_div(res, res, t, prec);
        }
        else if (ok && odd && kind != KIND_K && kind != KIND_I)
        {
            acb_neg(res, res);
            if (kind == KIND_JY)
                acb_neg(res2, res2);
        }

        acb_clear(t);
        acb_clear(t2);
        return ok;
    }

    if (crosses_cut(z))
    {
        bessel_param param;
        param.kind = kind;
        param.scaled = scaled;
        ok = _acb_hypgeom_split_at_cut(t, t2, nu, z, prec, bessel_func, &param);
    }
    else
        ok = bessel_nosplit(t, t2, kind, nu, z, scaled, prec);

    ok = ok && acb_is_finite(t) && acb_is_finite(t2);

    if (ok)
    {
        acb_set_round(res, t, prec);
        if (kind == KIND_JY)
            acb_set_round(res2, t2, prec);
    }

    acb_clear(t);
    acb_clear(t2);
    return ok;
}

/* Bounds for sup |F|, sup |dF/dz|, sup |dF/dnu| over the input balls
   from the L1 norms along the contour: dF/dzeta has the integrand
   sinh(w) exp(phi), with |sinh w| <= (e^x + e^(-x)) / 2, and dF/dnu the
   integrand -w exp(phi). */
int
_acb_hypgeom_bessel_integration_bounds(mag_t M, mag_t Dz, mag_t Dnu,
    int kind, const acb_t nu, const acb_t z)
{
    sd_problem P;
    sd_work W;
    sd_contour C;
    acb_t alpha;
    sd_end vm = sd_end_make(SD_MINUS, 0), vp0 = sd_end_make(SD_PLUS, 0), vp1 = sd_end_make(SD_PLUS, -1);
    double L[3], pre;
    int ok, rot, k;

    mag_inf(M);
    mag_inf(Dz);
    mag_inf(Dnu);

    if (!acb_is_finite(nu) || !acb_is_finite(z) || acb_contains_zero(z) || crosses_cut(z))
        return 0;

    acb_init(alpha);
    rot = bessel_rot(kind, z);
    bessel_alpha(alpha, z, rot);
    ok = _sd_problem_init(&P, alpha, 1, nu, NULL, BOUND_DROP, BOUND_TAILTOL, 1.0);
    acb_clear(alpha);
    if (!ok)
        return 0;

    _sd_work_init(&W, &P);

    if (kind == KIND_K)
        ok = _sd_plan(&C, &W, &P, vm, (rot == -1) ? vp1 : vp0, 1);
    else if (kind == KIND_Y)
        ok = _sd_plan(&C, &W, &P, vm, vp0, 1);
    else
        ok = _sd_plan(&C, &W, &P, vp1, vp0, 1);

    if (ok)
    {
        for (k = 0; k < 3; k++)
            L[k] = C.logL[k];
        _sd_contour_clear(&C);
    }

    /* Y = (H1 - H2) / (2 i) */
    if (ok && kind == KIND_Y)
    {
        ok = _sd_plan(&C, &W, &P, vm, vp1, 1);
        if (ok)
        {
            for (k = 0; k < 3; k++)
                L[k] = _sd_logaddexp(L[k], C.logL[k]);
            _sd_contour_clear(&C);
        }
    }

    _sd_work_clear(&W);

    if (!ok)
        return 0;

    /* the prefactors; for I and K, |e^{-+ i pi nu/2}|, whose
       nu-derivative adds (pi/2) |F| */
    pre = (kind == KIND_K) ? -log(2.0) : -log(2.0 * SD_PI);
    if (kind == KIND_I || kind == KIND_K)
    {
        pre += ((kind == KIND_I) ? rot : -rot) * SD_PI * P.nu.im / 2 + SD_PI * P.rnu / 2;
        L[2] = _sd_logaddexp(L[2], L[0] + log(SD_PI / 2));
    }

    _sd_mag_set_exp_d(M, L[0] + pre);
    _sd_mag_set_exp_d(Dz, L[1] + pre - log(2.0));
    _sd_mag_set_exp_d(Dnu, L[2] + pre);
    return 1;
}

static void
wrap(acb_t res, int kind, const acb_t nu, const acb_t z, slong prec)
{
    if (!_acb_hypgeom_bessel_integration(res, NULL, kind, nu, z, 0, prec))
        acb_indeterminate(res);
}

void acb_hypgeom_bessel_j_integration(acb_t res, const acb_t nu, const acb_t z, slong prec)
{ wrap(res, KIND_J, nu, z, prec); }

void acb_hypgeom_bessel_y_integration(acb_t res, const acb_t nu, const acb_t z, slong prec)
{ wrap(res, KIND_Y, nu, z, prec); }

void acb_hypgeom_bessel_h1_integration(acb_t res, const acb_t nu, const acb_t z, slong prec)
{ wrap(res, KIND_H1, nu, z, prec); }

void acb_hypgeom_bessel_h2_integration(acb_t res, const acb_t nu, const acb_t z, slong prec)
{ wrap(res, KIND_H2, nu, z, prec); }

void acb_hypgeom_bessel_i_integration(acb_t res, const acb_t nu, const acb_t z, slong prec)
{ wrap(res, KIND_I, nu, z, prec); }

void acb_hypgeom_bessel_k_integration(acb_t res, const acb_t nu, const acb_t z, slong prec)
{ wrap(res, KIND_K, nu, z, prec); }

void
acb_hypgeom_bessel_jy_integration(acb_t res1, acb_t res2, const acb_t nu, const acb_t z, slong prec)
{
    if (res1 == NULL)
        wrap(res2, KIND_Y, nu, z, prec);
    else if (!_acb_hypgeom_bessel_integration(res1, res2, KIND_JY, nu, z, 0, prec))
    {
        acb_indeterminate(res1);
        acb_indeterminate(res2);
    }
}
