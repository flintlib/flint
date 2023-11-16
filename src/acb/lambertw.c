/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

#ifdef __GNUC__
# define atan2 __builtin_atan2
# define cos __builtin_cos
# define fabs __builtin_fabs
# define log __builtin_log
# define pow __builtin_pow
# define sin __builtin_sin
#else
# include <math.h>
#endif

/* Check if z crosses a branch cut. */
int
acb_lambertw_branch_crossing(const acb_t z, const acb_t ez1, const fmpz_t k)
{
    if (arb_contains_zero(acb_imagref(z)) && !arb_is_nonnegative(acb_imagref(z)))
    {
        if (fmpz_is_zero(k))
        {
            if (!arb_is_positive(acb_realref(ez1)))
            {
                return 1;
            }
        }
        else if (!arb_is_positive(acb_realref(z)))
        {
            return 1;
        }
    }

    return 0;
}

/* todo: remove radii */
void
acb_lambertw_halley_step(acb_t res, acb_t ew, const acb_t z, const acb_t w, slong prec)
{
    acb_t t, u, v;

    acb_init(t);
    acb_init(u);
    acb_init(v);

    acb_exp(ew, w, prec);
    acb_add_ui(u, w, 2, prec);
    acb_add_ui(v, w, 1, prec);
    acb_mul_2exp_si(v, v, 1);
    acb_div(v, u, v, prec);
    acb_mul(t, ew, w, prec);
    acb_sub(u, t, z, prec);
    acb_mul(v, v, u, prec);
    acb_neg(v, v);
    acb_add(v, v, t, prec);
    acb_add(v, v, ew, prec);
    acb_div(t, u, v, prec);

    acb_sub(t, w, t, prec);

    acb_swap(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

/* assumes no aliasing of w and p */
void
acb_lambertw_branchpoint_series(acb_t w, const acb_t t, int bound, slong prec)
{
    slong i;
    static const int coeffs[] = {-130636800,130636800,-43545600,19958400,
        -10402560,5813640,-3394560,2042589,-1256320};

    acb_zero(w);

    for (i = 8; i >= 0; i--)
    {
        acb_mul(w, w, t, prec);
        acb_add_si(w, w, coeffs[i], prec);
    }

    acb_div_si(w, w, -coeffs[0], prec);

    if (bound)
    {
        mag_t err;
        mag_init(err);
        acb_get_mag(err, t);
        mag_geom_series(err, err, 9);

        if (acb_is_real(t))
            arb_add_error_mag(acb_realref(w), err);
        else
            acb_add_error_mag(w, err);
        mag_clear(err);
    }
}

void
acb_lambertw_principal_d(acb_t res, const acb_t z)
{
    double za, zb, wa, wb, ewa, ewb, t, u, q, r;
    int k, maxk = 15;

    za = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    zb = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    /* make sure we end up on the right branch */
    if (za < -0.367 && zb > -1e-20 && zb <= 0.0
                  && arf_sgn(arb_midref(acb_imagref(z))) < 0)
        zb = -1e-20;

    wa = za;
    wb = zb;

    if (fabs(wa) > 2.0 || fabs(wb) > 2.0)
    {
        t = atan2(wb, wa);
        wa = 0.5 * log(wa * wa + wb * wb);
        wb = t;
    }
    else if (fabs(wa) > 0.25 || fabs(wb) > 0.25)
    {
        /* We have W(z) ~= -1 + (2(ez+1))^(1/2) near the branch point.
           Changing the exponent to 1/4 gives a much worse local guess
           which however does the job on a larger domain. */
        wa *= 5.43656365691809;
        wb *= 5.43656365691809;
        wa += 2.0;
        t = atan2(wb, wa);
        r = pow(wa * wa + wb * wb, 0.125);
        wa = r * cos(0.25 * t);
        wb = r * sin(0.25 * t);
        wa -= 1.0;
    }

    for (k = 0; k < maxk; k++)
    {
#if __GNUC__
        t = __builtin_exp(wa);
#else
        t = exp(wa);
#endif
        ewa = t * cos(wb);
        ewb = t * sin(wb);
        t = (ewa * wa - ewb * wb); q = t + ewa; t -= za;
        u = (ewb * wa + ewa * wb); r = u + ewb; u -= zb;
        ewa = q * t + r * u; ewb = q * u - r * t;
        r = 1.0 / (q * q + r * r);
        ewa *= r; ewb *= r;
        if ((ewa*ewa + ewb*ewb) < (wa*wa + wb*wb) * 1e-12)
            maxk = FLINT_MIN(maxk, k + 2);
        wa -= ewa; wb -= ewb;
    }

    acb_set_d_d(res, wa, wb);
}

void
acb_lambertw_initial_asymp(acb_t w, const acb_t z, const fmpz_t k, slong prec)
{
    acb_t L1, L2, t;

    acb_init(L1);
    acb_init(L2);
    acb_init(t);

    acb_const_pi(L2, prec);
    acb_mul_2exp_si(L2, L2, 1);
    acb_mul_fmpz(L2, L2, k, prec);
    acb_mul_onei(L2, L2);
    acb_log(L1, z, prec);
    acb_add(L1, L1, L2, prec);
    acb_log(L2, L1, prec);

    /* L1 - L2 + L2/L1 + L2(L2-2)/(2 L1^2) */
    acb_inv(t, L1, prec);
    acb_mul_2exp_si(w, L2, 1);
    acb_submul(w, L2, L2, prec);
    acb_neg(w, w);
    acb_mul(w, w, t, prec);
    acb_mul_2exp_si(w, w, -1);
    acb_add(w, w, L2, prec);
    acb_mul(w, w, t, prec);
    acb_sub(w, w, L2, prec);
    acb_add(w, w, L1, prec);

    acb_clear(L1);
    acb_clear(L2);
    acb_clear(t);
}

/* assumes no aliasing */
slong
acb_lambertw_initial(acb_t res, const acb_t z, const acb_t ez1, const fmpz_t k, slong prec)
{
    /* Handle z very close to 0 on the principal branch. */
    if (fmpz_is_zero(k) &&
            (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -20) <= 0 &&
             arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -20) <= 0))
    {
        acb_set(res, z);
        acb_submul(res, res, res, prec);
        return 40;  /* could be tightened... */
    }

    /* For moderate input not close to the branch point, compute a double
       approximation as the initial value. */
    if (fmpz_is_zero(k) &&
        arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 400) < 0 &&
        arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 400) < 0 &&
          (arf_cmp_d(arb_midref(acb_realref(z)), -0.37) < 0 ||
           arf_cmp_d(arb_midref(acb_realref(z)), -0.36) > 0 ||
           arf_cmpabs_d(arb_midref(acb_imagref(z)), 0.01) > 0))
    {
        acb_lambertw_principal_d(res, z);
        return 48;
    }

    /* Check if we are close to the branch point at -1/e. */
    if ((fmpz_is_zero(k) || (fmpz_is_one(k) && arb_is_negative(acb_imagref(z)))
                         || (fmpz_equal_si(k, -1) && arb_is_nonnegative(acb_imagref(z))))
        && ((arf_cmpabs_2exp_si(arb_midref(acb_realref(ez1)), -2) <= 0 &&
             arf_cmpabs_2exp_si(arb_midref(acb_imagref(ez1)), -2) <= 0)))
    {
        acb_t t;
        acb_init(t);
        acb_mul_2exp_si(t, ez1, 1);
        mag_zero(arb_radref(acb_realref(t)));
        mag_zero(arb_radref(acb_imagref(t)));
        acb_mul_ui(t, t, 3, prec);
        acb_sqrt(t, t, prec);
        if (!fmpz_is_zero(k))
            acb_neg(t, t);
        acb_lambertw_branchpoint_series(res, t, 0, prec);
        acb_clear(t);
        return 1;  /* todo: estimate */
    }

    acb_lambertw_initial_asymp(res, z, k, prec);
    return 1;  /* todo: estimate */
}

/* note: z should be exact here */
void acb_lambertw_main(acb_t res, const acb_t z,
                const acb_t ez1, const fmpz_t k, int flags, slong prec)
{
    acb_t w, t, oldw, ew;
    mag_t err;
    slong i, wp, accuracy, ebits, kbits, mbits, wp_initial, extraprec;
    int have_ew;

    acb_init(t);
    acb_init(w);
    acb_init(oldw);
    acb_init(ew);
    mag_init(err);

    /* We need higher precision for large k, large exponents, or very close
       to the branch point at -1/e. todo: we should be recomputing
       ez1 to higher precision when close... */
    acb_get_mag(err, z);
    if (fmpz_is_zero(k) && mag_cmp_2exp_si(err, 0) < 0)
        ebits = 0;
    else
        ebits = fmpz_bits(MAG_EXPREF(err));

    if (fmpz_is_zero(k) || (fmpz_is_one(k) && arb_is_negative(acb_imagref(z)))
                        || (fmpz_equal_si(k, -1) && arb_is_nonnegative(acb_imagref(z))))
    {
        acb_get_mag(err, ez1);
        mbits = -MAG_EXP(err);
        mbits = FLINT_MAX(mbits, 0);
        mbits = FLINT_MIN(mbits, prec);
    }
    else
    {
        mbits = 0;
    }

    kbits = fmpz_bits(k);

    extraprec = FLINT_MAX(ebits, kbits);
    extraprec = FLINT_MAX(extraprec, mbits);

    wp = wp_initial = 40 + extraprec;

    accuracy = acb_lambertw_initial(w, z, ez1, k, wp_initial);
    mag_zero(arb_radref(acb_realref(w)));
    mag_zero(arb_radref(acb_imagref(w)));

    /* We should be able to compute e^w for the final certification
       during the Halley iteration. */
    have_ew = 0;

    for (i = 0; i < 5 + FLINT_BIT_COUNT(prec + extraprec); i++)
    {
        /* todo: should we restart? */
        if (!acb_is_finite(w))
            break;

        wp = FLINT_MIN(3 * accuracy, 1.1 * prec + 10);
        wp = FLINT_MAX(wp, 40);
        wp += extraprec;

        acb_set(oldw, w);
        acb_lambertw_halley_step(t, ew, z, w, wp);

        /* estimate the error (conservatively) */
        acb_sub(w, w, t, wp);
        acb_get_mag(err, w);
        acb_set(w, t);
        acb_add_error_mag(t, err);
        accuracy = acb_rel_accuracy_bits(t);

        if (accuracy > 2 * extraprec)
            accuracy *= 2.9;  /* less conservatively */

        accuracy = FLINT_MIN(accuracy, wp);
        accuracy = FLINT_MAX(accuracy, 0);

        if (accuracy > prec + extraprec)
        {
            /* e^w = e^oldw * e^(w-oldw) */
            acb_sub(t, w, oldw, wp);
            acb_exp(t, t, wp);
            acb_mul(ew, ew, t, wp);
            have_ew = 1;
            break;
        }

        mag_zero(arb_radref(acb_realref(w)));
        mag_zero(arb_radref(acb_imagref(w)));
    }

    wp = FLINT_MIN(3 * accuracy, 1.1 * prec + 10);
    wp = FLINT_MAX(wp, 40);
    wp += extraprec;

    if (acb_lambertw_check_branch(w, k, wp))
    {
        acb_t u, r, eu1;
        mag_t err, rad;

        acb_init(u);
        acb_init(r);
        acb_init(eu1);

        mag_init(err);
        mag_init(rad);

        if (have_ew)
            acb_set(t, ew);
        else
            acb_exp(t, w, wp);
        /* t = w e^w */
        acb_mul(t, t, w, wp);

        acb_sub(r, t, z, wp);

        /* Bound W' on the straight line path between t and z */
        acb_union(u, t, z, wp);

        arb_const_e(acb_realref(eu1), wp);
        arb_zero(acb_imagref(eu1));
        acb_mul(eu1, eu1, u, wp);
        acb_add_ui(eu1, eu1, 1, wp);

        if (acb_lambertw_branch_crossing(u, eu1, k))
        {
            mag_inf(err);
        }
        else
        {
            acb_lambertw_bound_deriv(err, u, eu1, k);
            acb_get_mag(rad, r);
            mag_mul(err, err, rad);
        }

        acb_add_error_mag(w, err);

        acb_set(res, w);

        acb_clear(u);
        acb_clear(r);
        acb_clear(eu1);
        mag_clear(err);
        mag_clear(rad);
    }
    else
    {
        acb_indeterminate(res);
    }

    acb_clear(t);
    acb_clear(w);
    acb_clear(oldw);
    acb_clear(ew);
    mag_clear(err);
}

void
acb_lambertw_cleared_cut(acb_t res, const acb_t z, const fmpz_t k, int flags, slong prec)
{
    acb_t ez1;
    acb_init(ez1);

    /* compute e*z + 1 */
    arb_const_e(acb_realref(ez1), prec);
    acb_mul(ez1, ez1, z, prec);
    acb_add_ui(ez1, ez1, 1, prec);

    if (acb_is_exact(z))
    {
        acb_lambertw_main(res, z, ez1, k, flags, prec);
    }
    else
    {
        acb_t zz;
        mag_t err, rad;

        mag_init(err);
        mag_init(rad);
        acb_init(zz);

        acb_lambertw_bound_deriv(err, z, ez1, k);
        mag_hypot(rad, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
        mag_mul(err, err, rad);

        acb_set(zz, z);
        mag_zero(arb_radref(acb_realref(zz)));
        mag_zero(arb_radref(acb_imagref(zz)));  /* todo: recompute ez1? */

        acb_lambertw_main(res, zz, ez1, k, flags, prec);
        acb_add_error_mag(res, err);

        mag_clear(err);
        mag_clear(rad);
        acb_clear(zz);
    }

    acb_clear(ez1);
}

/* Extremely close to the branch point at -1/e, use the series expansion directly. */
int
acb_lambertw_try_near_branch_point(acb_t res, const acb_t z,
    const acb_t ez1, const fmpz_t k, int flags, slong prec)
{
    if (fmpz_is_zero(k) || (fmpz_is_one(k) && arb_is_negative(acb_imagref(z)))
                        || (fmpz_equal_si(k, -1) && arb_is_nonnegative(acb_imagref(z))))
    {
        if (acb_contains_zero(ez1) ||
            (arf_cmpabs_2exp_si(arb_midref(acb_realref(ez1)), -prec / 4.5 - 6) < 0 &&
             arf_cmpabs_2exp_si(arb_midref(acb_imagref(ez1)), -prec / 4.5 - 6) < 0))
        {
            acb_t t;
            acb_init(t);
            acb_mul_2exp_si(t, ez1, 1);
            acb_sqrt(t, t, prec);
            if (!fmpz_is_zero(k))
                acb_neg(t, t);
            acb_lambertw_branchpoint_series(res, t, 1, prec);
            acb_clear(t);
            return 1;
        }
    }

    return 0;
}

void
acb_lambertw_cleared_cut_fix_small(acb_t res, const acb_t z,
    const acb_t ez1, const fmpz_t k, int flags, slong prec)
{
    acb_t zz, zmid, zmide1;
    arf_t eps;

    acb_init(zz);
    acb_init(zmid);
    acb_init(zmide1);
    arf_init(eps);

    arf_mul_2exp_si(eps, arb_midref(acb_realref(z)), -prec);
    acb_set(zz, z);

    if (arf_sgn(arb_midref(acb_realref(zz))) < 0 &&
        (!fmpz_is_zero(k) || arf_sgn(arb_midref(acb_realref(ez1))) < 0) &&
        arf_cmpabs(arb_midref(acb_imagref(zz)), eps) < 0)
    {
        /* now the value must be in [0,2eps] */
        arf_get_mag(arb_radref(acb_imagref(zz)), eps);
        arf_set_mag(arb_midref(acb_imagref(zz)), arb_radref(acb_imagref(zz)));

        if (arf_sgn(arb_midref(acb_imagref(z))) >= 0)
        {
            acb_lambertw_cleared_cut(res, zz, k, flags, prec);
        }
        else
        {
            fmpz_t kk;
            fmpz_init(kk);
            fmpz_neg(kk, k);
            acb_lambertw_cleared_cut(res, zz, kk, flags, prec);
            acb_conj(res, res);
            fmpz_clear(kk);
        }
    }
    else
    {
        acb_lambertw_cleared_cut(res, zz, k, flags, prec);
    }

    acb_clear(zz);
    acb_clear(zmid);
    acb_clear(zmide1);
    arf_clear(eps);
}

void
_acb_lambertw(acb_t res, const acb_t z, const acb_t ez1, const fmpz_t k, int flags, slong prec)
{
    slong goal, ebits, ebits2, ls, lt;
    const fmpz * expo;

    /* Estimated accuracy goal. */
    /* todo: account for exponent bits and bits in k. */
    goal = acb_rel_accuracy_bits(z);
    goal = FLINT_MAX(goal, 10);
    goal = FLINT_MIN(goal, prec);

    /* Handle tiny z directly. For k >= 2, |c_k| <= 4^k / 16. */
    if (fmpz_is_zero(k)
        && arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -goal / 2) < 0
        && arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -goal / 2) < 0)
    {
        mag_t err;
        mag_init(err);
        acb_get_mag(err, z);
        mag_mul_2exp_si(err, err, 2);
        acb_set(res, z);
        acb_submul(res, res, res, prec);
        mag_geom_series(err, err, 3);
        mag_mul_2exp_si(err, err, -4);
        acb_add_error_mag(res, err);
        mag_clear(err);
        return;
    }

    if (arf_cmpabs(arb_midref(acb_realref(z)), arb_midref(acb_imagref(z))) >= 0)
        expo = ARF_EXPREF(arb_midref(acb_realref(z)));
    else
        expo = ARF_EXPREF(arb_midref(acb_imagref(z)));

    ebits = fmpz_bits(expo);

    /* ebits ~= log2(|log(z) + 2 pi i k|) */
    /* ebits2 ~= log2(log(log(z))) */
    ebits = FLINT_MAX(ebits, fmpz_bits(k));
    ebits = FLINT_MAX(ebits, 1) - 1;
    ebits2 = FLINT_BIT_COUNT(ebits);
    ebits2 = FLINT_MAX(ebits2, 1) - 1;

    /* We gain accuracy from the exponent when W ~ log - log log */
    if (fmpz_sgn(expo) > 0 || (fmpz_sgn(expo) < 0 && !fmpz_is_zero(k)))
    {
        goal += ebits - ebits2;
        goal = FLINT_MAX(goal, 10);
        goal = FLINT_MIN(goal, prec);

        /* The asymptotic series with truncation L, M gives us about
           t - max(2+lt+L*(2+ls), M*(2+lt)) bits of accuracy where
           ls = -ebits, lt = ebits2 - ebits. */
        ls = 2 - ebits;
        lt = 2 + ebits2 - ebits;

        if (ebits - FLINT_MAX(lt + 1*ls, 1*lt) > goal)
        {
            acb_lambertw_asymp(res, z, k, 1, 1, goal);
            acb_set_round(res, res, prec);
            return;
        }
        else if (ebits - FLINT_MAX(lt + 3*ls, 5*lt) > goal)
        {
            acb_lambertw_asymp(res, z, k, 3, 5, goal);
            acb_set_round(res, res, prec);
            return;
        }
    }

    /* Extremely close to the branch point at -1/e, use the series expansion directly. */
    if (acb_lambertw_try_near_branch_point(res, z, ez1, k, flags, goal))
    {
        acb_set_round(res, res, prec);
        return;
    }

    /* compute union of both sides */
    if (acb_lambertw_branch_crossing(z, ez1, k))
    {
        acb_t za, zb, eza1, ezb1;
        fmpz_t kk;

        acb_init(za);
        acb_init(zb);
        acb_init(eza1);
        acb_init(ezb1);
        fmpz_init(kk);

        fmpz_neg(kk, k);

        acb_set(za, z);
        acb_conj(zb, z);
        arb_nonnegative_part(acb_imagref(za), acb_imagref(za));
        arb_nonnegative_part(acb_imagref(zb), acb_imagref(zb));

        acb_set(eza1, ez1);
        acb_conj(ezb1, ez1);
        arb_nonnegative_part(acb_imagref(eza1), acb_imagref(eza1));
        arb_nonnegative_part(acb_imagref(ezb1), acb_imagref(ezb1));

        /* Check series expansion again, because now there is no crossing. */
        if (!acb_lambertw_try_near_branch_point(res, za, eza1, k, flags, goal))
            acb_lambertw_cleared_cut_fix_small(za, za, eza1, k, flags, goal);

        if (!acb_lambertw_try_near_branch_point(res, zb, ezb1, kk, flags, goal))
            acb_lambertw_cleared_cut_fix_small(zb, zb, ezb1, kk, flags, goal);

        acb_conj(zb, zb);
        acb_union(res, za, zb, prec);

        acb_clear(za);
        acb_clear(zb);
        acb_clear(eza1);
        acb_clear(ezb1);
        fmpz_clear(kk);
    }
    else
    {
        acb_lambertw_cleared_cut_fix_small(res, z, ez1, k, flags, goal);
        acb_set_round(res, res, prec);
    }
}

void
acb_lambertw_middle(acb_t res, const acb_t z, slong prec)
{
    fmpz_t k;

    if (acb_contains_zero(z))
    {
        acb_indeterminate(res);
        return;
    }

    fmpz_init(k);
    fmpz_set_si(k, -1);

    if (arb_is_positive(acb_imagref(z)))
    {
        acb_lambertw(res, z, k, 0, prec);
    }
    else if (arb_is_negative(acb_imagref(z)))
    {
        acb_conj(res, z);
        acb_lambertw(res, res, k, 0, prec);
        acb_conj(res, res);
    }
    else if (arb_is_negative(acb_realref(z)))
    {
        if (arb_is_nonnegative(acb_imagref(z)))
        {
            acb_lambertw(res, z, k, 0, prec);
        }
        else if (arb_is_negative(acb_imagref(z)))
        {
            acb_conj(res, z);
            acb_lambertw(res, res, k, 0, prec);
            acb_conj(res, res);
        }
        else
        {
            acb_t za, zb;
            acb_init(za);
            acb_init(zb);
            acb_set(za, z);
            acb_conj(zb, z);
            arb_nonnegative_part(acb_imagref(za), acb_imagref(za));
            arb_nonnegative_part(acb_imagref(zb), acb_imagref(zb));
            acb_lambertw(za, za, k, 0, prec);
            acb_lambertw(zb, zb, k, 0, prec);
            acb_conj(zb, zb);
            acb_union(res, za, zb, prec);
            acb_clear(za);
            acb_clear(zb);
        }
    }
    else /* re is positive */
    {
        if (arb_is_positive(acb_imagref(z)))
        {
            acb_lambertw(res, z, k, 0, prec);
        }
        else if (arb_is_nonpositive(acb_imagref(z)))
        {
            acb_conj(res, z);
            acb_lambertw(res, res, k, 0, prec);
            acb_conj(res, res);
        }
        else
        {
            acb_t za, zb;
            acb_init(za);
            acb_init(zb);
            acb_set(za, z);
            acb_conj(zb, z);
            arb_nonnegative_part(acb_imagref(za), acb_imagref(za));
            arb_nonnegative_part(acb_imagref(zb), acb_imagref(zb));
            acb_lambertw(za, za, k, 0, prec);
            acb_lambertw(zb, zb, k, 0, prec);
            acb_conj(zb, zb);
            acb_union(res, za, zb, prec);
            acb_clear(za);
            acb_clear(zb);
        }
    }

    fmpz_clear(k);
}

void
acb_lambertw_left(acb_t res, const acb_t z, const fmpz_t k, slong prec)
{
    if (acb_contains_zero(z) && !(fmpz_equal_si(k, -1) && acb_is_real(z)))
    {
        acb_indeterminate(res);
        return;
    }

    if (arb_is_positive(acb_imagref(z)))
    {
        acb_lambertw(res, z, k, 0, prec);
    }
    else if (arb_is_nonpositive(acb_imagref(z)))
    {
        fmpz_t kk;
        fmpz_init(kk);
        fmpz_add_ui(kk, k, 1);
        fmpz_neg(kk, kk);

        acb_conj(res, z);
        acb_lambertw(res, res, kk, 0, prec);
        acb_conj(res, res);

        fmpz_clear(kk);
    }
    else
    {
        acb_t za, zb;
        fmpz_t kk;

        acb_init(za);
        acb_init(zb);
        fmpz_init(kk);

        acb_set(za, z);
        acb_conj(zb, z);

        arb_nonnegative_part(acb_imagref(za), acb_imagref(za));
        arb_nonnegative_part(acb_imagref(zb), acb_imagref(zb));

        fmpz_add_ui(kk, k, 1);
        fmpz_neg(kk, kk);

        acb_lambertw(za, za, k, 0, prec);
        acb_lambertw(zb, zb, kk, 0, prec);
        acb_conj(zb, zb);

        acb_union(res, za, zb, prec);

        acb_clear(za);
        acb_clear(zb);
        fmpz_clear(kk);
    }
}

void
acb_lambertw(acb_t res, const acb_t z, const fmpz_t k, int flags, slong prec)
{
    acb_t ez1;

    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (flags == ACB_LAMBERTW_LEFT)
    {
        acb_lambertw_left(res, z, k, prec);
        return;
    }

    if (flags == ACB_LAMBERTW_MIDDLE)
    {
        acb_lambertw_middle(res, z, prec);
        return;
    }

    if (acb_contains_zero(z) && !fmpz_is_zero(k))
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(ez1);

    /* precompute z*e + 1 */
    arb_const_e(acb_realref(ez1), prec);
    acb_mul(ez1, ez1, z, prec);
    acb_add_ui(ez1, ez1, 1, prec);

    /* Compute standard branches */

    /* use real code when possible */
    if (acb_is_real(z) && arb_is_positive(acb_realref(ez1)) &&
        (fmpz_is_zero(k) ||
        (fmpz_equal_si(k, -1) && arb_is_negative(acb_realref(z)))))
    {
        arb_lambertw(acb_realref(res), acb_realref(z), !fmpz_is_zero(k), prec);
        arb_zero(acb_imagref(res));
    }
    else
    {
        _acb_lambertw(res, z, ez1, k, flags, prec);
    }

    acb_clear(ez1);
}

