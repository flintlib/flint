/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "fmpq.h"
#include "arb_hypgeom.h"

#define LOG2 0.69314718055994530942
#define EXP1 2.7182818284590452354
#define INV_LOG2 1.4426950408889634074

void
arb_hypgeom_erf_one_eps(arb_t res, const arb_t z)
{
    mag_t t, u;
    mag_init(t);
    mag_init(u);

    arb_get_mag_lower(t, z);
    mag_mul_lower(u, t, t);
    mag_expinv(u, u);
    mag_div(u, u, t);

    /* 1/sqrt(pi) < 289/512 */
    mag_mul_ui(u, u, 289);
    mag_mul_2exp_si(arb_radref(res), u, -9);

    if (mag_cmp_2exp_si(arb_radref(res), 1) > 0)
    {
        mag_one(arb_radref(res));
        mag_mul_2exp_si(arb_radref(res), arb_radref(res), 2);
    }

    if (arf_sgn(arb_midref(z)) > 0)
        arf_one(arb_midref(res));
    else
    {
        arf_one(arb_midref(res));
        arf_neg(arb_midref(res), arb_midref(res));
    }

    mag_clear(t);
    mag_clear(u);
}

void
arb_hypgeom_erf_propagated_error(mag_t err, const arb_t z)
{
    mag_t x;
    mag_init(x);

    /* exp(-z^2) */
    arb_get_mag_lower(x, z);
    mag_mul_lower(x, x, x);
    mag_expinv(err, x);
    mag_mul(err, err, arb_radref(z));

    /* 2/sqrt(pi) < 289/256 */
    mag_mul_ui(err, err, 289);
    mag_mul_2exp_si(err, err, -8);

    /* |erf(a) - erf(b)| <= 2 */
    mag_set_ui(x, 2);
    mag_min(err, err, x);

    mag_clear(x);
}

void
arb_hypgeom_erf_1f1b(arb_t res, const arb_t z, slong prec)
{
    arb_t t, u;
    slong N;
    mag_t err;

    arb_init(t);
    arb_init(u);
    mag_init(err);

    if (arf_cmpabs_2exp_si(arb_midref(z), -32) < 0)
    {
        if (arf_cmpabs_2exp_si(arb_midref(z), -prec) < 0)
            N = 1;
        else
            N = -prec / (2 * ARF_EXP(arb_midref(z))) + 1;
    }
    else
    {
        double u, dz;
        dz = arf_get_d(arb_midref(z), ARF_RND_DOWN);
        dz = fabs(dz);
        u = -dz * dz + prec * LOG2 + log(dz);
        if (dz < 1.0)
          u = FLINT_MAX(u, 1e-6);
        u = u / d_lambertw(u / (EXP1 * dz * dz));
        N = u + 1;
    }

    N = FLINT_MAX(N, 1);

    arb_sqr(t, z, prec);
    _arb_hypgeom_gamma_lower_sum_rs_1(u, 3, 2, t, N, prec);

    /* z^(2k) / rf(3/2,k) <= (z^2)^k / k! */
    arb_get_mag(err, t);
    mag_exp_tail(err, err, N);
    arb_add_error_mag(u, err);

    arb_neg(t, t);
    arb_exp(t, t, prec);
    arb_mul(u, u, t, prec);
    arb_const_sqrt_pi(t, prec);
    arb_div(u, u, t, prec);
    arb_mul(u, u, z, prec);
    arb_mul_2exp_si(res, u, 1);

    arb_clear(t);
    arb_clear(u);
    mag_clear(err);
}

void
arb_hypgeom_erf_asymp(arb_t res, const arb_t z, slong N, int complementary, slong prec, slong prec2)
{
    arb_t t, u;
    int sgn;
    mag_t err, tm;

    if (!arb_is_exact(z) && (arf_cmpabs_ui(arb_midref(z), prec) < 0 ||
        (complementary && arb_rel_accuracy_bits(z) < prec)))
    {
        arb_t zmid;
        mag_t err;

        arb_init(zmid);
        mag_init(err);

        arb_hypgeom_erf_propagated_error(err, z);
        arf_set(arb_midref(zmid), arb_midref(z));
        arb_hypgeom_erf_asymp(res, zmid, N, complementary, prec, prec2);
        arb_add_error_mag(res, err);

        arb_clear(zmid);
        mag_clear(err);
        return;
    }

    arb_init(t);
    arb_init(u);
    mag_init(err);
    mag_init(tm);

    sgn = arf_sgn(arb_midref(z));

    arb_sqr(t, z, prec2);
    arb_neg(t, t);
    _arb_hypgeom_gamma_upper_sum_rs_1(u, 1, 2, t, N, prec2);

    /* Error is bounded by first omitted term, rf(1/2,N) / z^(2N) <= N! / z^(2N) */
    arb_get_mag_lower(err, t);
    mag_inv(err, err);
    mag_pow_ui(err, err, N);
    mag_fac_ui(tm, N);
    mag_mul(err, err, tm);
    arb_add_error_mag(u, err);

    arb_exp(t, t, prec2);
    arb_mul(u, u, t, prec2);
    arb_const_sqrt_pi(t, prec2);
    arb_mul(t, t, z, prec2);
    arb_div(res, u, t, prec2);

    if (!complementary)
    {
        if (sgn == 1)
            arb_sub_ui(res, res, 1, prec);
        else
            arb_add_ui(res, res, 1, prec);

        arb_neg(res, res);
    }

    arb_clear(t);
    arb_clear(u);
    mag_clear(err);
    mag_clear(tm);
}

void
arb_hypgeom_erf_1f1(arb_t res, const arb_t z, slong prec, slong wp)
{
    if (arb_rel_accuracy_bits(z) >= wp)
    {
        arb_hypgeom_erf_1f1b(res, z, wp);
    }
    else
    {
        arb_t zmid;
        mag_t err;

        arb_init(zmid);
        mag_init(err);

        arb_hypgeom_erf_propagated_error(err, z);
        arf_set(arb_midref(zmid), arb_midref(z));

        arb_hypgeom_erf_1f1b(res, zmid, wp);
        arb_add_error_mag(res, err);

        arb_clear(zmid);
        mag_clear(err);
    }

    arb_set_round(res, res, prec);
}

static void
_arf_trunc(arf_t x)
{
    if (arf_sgn(x) < 0)
        arf_ceil(x, x);
    else
        arf_floor(x, x);
}

static void
arb_extract_bits(arb_t t, const arb_t z, slong b)
{
    arb_mul_2exp_si(t, z, b);
    _arf_trunc(arb_midref(t));
    mag_zero(arb_radref(t));
    arb_mul_2exp_si(t, t, -b);
}

/* Compute Gamma(a,z) using the bit-burst algorithm.
   Todo: allow passing precomputed Gamma(a) as input. */
void
_arb_gamma_upper_fmpq_bb(arb_t res, const fmpq_t a, const arb_t z, const mag_t abs_tol, slong prec_lower, slong prec_upper)
{
    slong start_bits, bits, wp, NN;
    arb_t Gz0, Gz1, z0, z1, expmz0, t;
    mag_t AE;

    arb_init(t);
    arb_init(z0);
    arb_init(z1);
    arb_init(Gz0);
    arb_init(Gz1);
    arb_init(expmz0);
    mag_init(AE);

    start_bits = 64;
    wp = prec_upper;

    /* Hack: the error bound for the local Taylor series assumes that the
       step size is much smaller than the expansion point, even when the
       expansion point is close to zero. So when close to zero, we need
       to take more initial bits. */
    while (arf_cmpabs_2exp_si(arb_midref(z), -start_bits / 4) < 0)
    {
        if (start_bits > prec_lower)
        {
            NN = _arb_hypgeom_gamma_lower_fmpq_0_choose_N(AE, a, z, abs_tol);
            _arb_hypgeom_gamma_lower_fmpq_0_bsplit(Gz0, a, z, NN, prec_lower);
            arb_add_error_mag(Gz0, AE);
            arb_gamma_fmpq(t, a, FLINT_MAX(prec_lower, prec_upper));
            arb_sub(res, t, Gz0, prec_upper);
            goto bb_cleanup;
        }

        start_bits *= 2;
    }

    arb_extract_bits(z0, z, start_bits);

    NN = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N(AE, a, z0, abs_tol);

    if (NN != -1)
    {
        _arb_hypgeom_gamma_upper_fmpq_inf_bsplit(Gz0, a, z0, NN, wp);
        arb_add_error_mag(Gz0, AE);
    }
    else
    {
        NN = _arb_hypgeom_gamma_lower_fmpq_0_choose_N(AE, a, z0, abs_tol);
        _arb_hypgeom_gamma_lower_fmpq_0_bsplit(Gz0, a, z0, NN, prec_lower);
        arb_add_error_mag(Gz0, AE);
        arb_gamma_fmpq(t, a, FLINT_MAX(prec_lower, prec_upper));
        arb_sub(Gz0, t, Gz0, prec_upper);
    }

    arb_neg(expmz0, z0);
    arb_exp(expmz0, expmz0, wp);

    for (bits = start_bits * 2; bits < wp / 8; bits *= 2)
    {
        arb_extract_bits(z1, z, bits);
        _arb_gamma_upper_fmpq_step_bsplit(Gz1, a, z0, z1, Gz0, expmz0, abs_tol, wp);
        arb_sub(t, z0, z1, wp);
        arb_exp(t, t, wp);
        arb_mul(expmz0, expmz0, t, wp);
        arb_set(Gz0, Gz1);
        arb_set(z0, z1);
    }

    /* Final step, including error bound */
    _arb_gamma_upper_fmpq_step_bsplit(Gz1, a, z0, z, Gz0, expmz0, abs_tol, wp);

    arb_set(res, Gz1);

bb_cleanup:
    arb_clear(t);
    arb_clear(z0);
    arb_clear(z1);
    arb_clear(Gz0);
    arb_clear(Gz1);
    arb_clear(expmz0);
    mag_clear(AE);
}

int
arb_hypgeom_erf_bb(arb_t res, const arb_t z, int complementary, slong prec)
{
    mag_t tol, tm;
    double x;
    arb_t t;
    fmpq_t a;
    slong wp_lower, wp_upper;
    int sgn;

    /* Avoid bit-burst algorithm for huge input and very close to 0. */
    /* With better error bounds, this exit shouldn't be necessary. */
    if (!arb_is_finite(z) ||
        arf_cmpabs_ui(arb_midref(z), prec) > 0 ||
        arf_cmpabs_2exp_si(arb_midref(z), -prec / 16) < 0)
    {
        return 0;
    }

    sgn = arf_sgn(arb_midref(z));
    x = arf_get_d(arb_midref(z), ARF_RND_DOWN);
    x = fabs(x);

    if (!arb_is_exact(z))
    {
        arb_t zmid;
        mag_t err;
        int success;

        arb_init(zmid);
        mag_init(err);

        arb_hypgeom_erf_propagated_error(err, z);
        arf_set(arb_midref(zmid), arb_midref(z));

        success = arb_hypgeom_erf_bb(res, zmid, complementary, prec);
        if (success)
            arb_add_error_mag(res, err);

        arb_clear(zmid);
        mag_clear(err);
        return success;
    }

    mag_init(tol);
    mag_init(tm);
    arb_init(t);
    fmpq_init(a);

    /* Near 0, need to convert relative to absolute precision for erf */
    if (x < 0.25 && !complementary)
    {
        wp_lower = prec + 20 + 0.001 * prec;

        arb_get_mag(tol, z);
        mag_mul_2exp_si(tol, tol, -wp_lower);

        wp_upper = wp_lower + (-MAG_EXP(tol));
    }
    else if (complementary && sgn == 1 && x > 1.0)
    {
        wp_upper = prec + 20 + 0.001 * prec;

        /* We will have cancellation with the lower series */
        arb_get_mag_lower(tm, z);
        mag_mul(tol, tm, tm);
        mag_expinv(tol, tol);
        mag_div(tol, tol, tm);
        mag_mul_2exp_si(tol, tol, -wp_upper);

        wp_lower = wp_upper + x * x * INV_LOG2;
        if (x >= 1)
            wp_lower = wp_lower - log(x) * INV_LOG2;

        wp_lower = FLINT_MAX(wp_lower, 30);
        wp_upper = FLINT_MAX(wp_upper, 30);
    }
    else
    {
        wp_lower = prec + 20 + 0.001 * prec;
        wp_upper = wp_lower;

        /* Can reduce precision with the upper series */
        mag_set_ui_2exp_si(tol, 1, -wp_lower);

        if (x >= 1)
            wp_upper = wp_upper - x * x * INV_LOG2 - log(x) * INV_LOG2;

        wp_upper = FLINT_MAX(wp_upper, 30);
    }

    fmpq_set_si(a, 1, 2);
    arb_sqr(t, z, FLINT_MAX(wp_lower, wp_upper));
    _arb_gamma_upper_fmpq_bb(res, a, t, tol, wp_lower, wp_upper);

    arb_const_sqrt_pi(t, wp_upper);
    arb_div(res, res, t, wp_upper);

    if (complementary)
    {
        if (sgn < 0)
        {
            arb_sub_ui(res, res, 2, prec);
            arb_neg(res, res);
        }
    }
    else
    {
        arb_sub_ui(res, res, 1, prec);
        if (sgn > 0)
            arb_neg(res, res);
    }

    mag_clear(tol);
    mag_clear(tm);
    arb_clear(t);
    fmpq_clear(a);

    return 1;
}

void
arb_hypgeom_erf(arb_t res, const arb_t z, slong prec)
{
    double abs_z2, log_z;
    double log2_err, err_prev;
    slong acc, prec2, wp, N;
    double x;

    if (!arb_is_finite(z))
    {
        arb_indeterminate(res);
        return;
    }

    if (arb_is_zero(z))
    {
        arb_zero(res);
        return;
    }

    if (arf_cmpabs_2exp_si(arb_midref(z), -prec / 16) < 0)
    {
        wp = prec + 20 + FLINT_BIT_COUNT(prec);
        arb_hypgeom_erf_1f1(res, z, prec, wp);
        return;
    }

    if (arf_cmpabs_2exp_si(arb_midref(z), 60) > 0)
    {
        arb_hypgeom_erf_one_eps(res, z);
        return;
    }

    /* exp(-z^2) / z * N! / z^(2N) < 2^p */
    x = arf_get_d(arb_midref(z), ARF_RND_DOWN);
    x = fabs(x);

    acc = arb_rel_accuracy_bits(z);
    acc = FLINT_MAX(acc, 0);
    acc = FLINT_MIN(acc, prec);
    prec = FLINT_MIN(prec, acc + x * x * INV_LOG2 + 32);

    if (x * x * INV_LOG2 > prec)
    {
        arb_hypgeom_erf_one_eps(res, z);
        return;
    }

    if (prec > 30000 && x > 150.0 / exp(4e-3 * sqrt(prec)) && x < 0.6 * sqrt(prec))
    {
        if (arb_hypgeom_erf_bb(res, z, 0, prec))
            return;
    }

    /* Can we use the asymptotic expansion? */
    if (x > 2.0)
    {
        prec2 = prec + 5 + FLINT_BIT_COUNT(prec);

        abs_z2 = x * x;
        log_z = 0.5 * log(abs_z2);

        if ((x * x + log_z) * INV_LOG2 > prec)
        {
            arb_hypgeom_erf_one_eps(res, z);
            return;
        }

        wp = prec - x * x * INV_LOG2 - log_z * INV_LOG2 + 10;
        wp = FLINT_MAX(wp, 30);

        err_prev = 0.0;
        for (N = 1; ; N++)
        {
            log2_err = -x * x - (2 * N + 1) * log_z + N * (log(N) - 1.0);
            log2_err *= INV_LOG2;

            if (log2_err > err_prev)
                break;

            if (log2_err < -prec2 - 10)
            {
                arb_hypgeom_erf_asymp(res, z, N, 0, prec, wp);
                return;
            }

            err_prev = log2_err;
        }
    }

    wp = prec + 10 + FLINT_BIT_COUNT(prec);
    arb_hypgeom_erf_1f1(res, z, prec, wp);
}

void
arb_hypgeom_erfc(arb_t res, const arb_t z, slong prec)
{
    double x, abs_z2, log_z;
    slong acc, prec2, wp, N;

    if (!arb_is_finite(z))
    {
        arb_indeterminate(res);
        return;
    }

    if (arb_is_zero(z))
    {
        arb_one(res);
        return;
    }

    if (arf_cmp_si(arb_midref(z), 1) <= 0)
    {
        arb_hypgeom_erf(res, z, prec + 5);
        arb_sub_ui(res, res, 1, prec);
        arb_neg(res, res);
        return;
    }

    acc = arb_rel_accuracy_bits(z);
    acc = FLINT_MAX(acc, 0);
    acc = FLINT_MIN(acc, prec);
    prec = FLINT_MIN(prec, acc + 32);

    /* Super-huge -- we only need one term of the asymptotic expansion */
    if (arf_cmpabs_2exp_si(arb_midref(z), prec / 2 + 10) > 0)
    {
        arb_hypgeom_erf_asymp(res, z, 1, 1, prec, prec);
        return;
    }

    x = arf_get_d(arb_midref(z), ARF_RND_DOWN);
    x = fabs(x);

    if (prec > 30000 && x > 150.0 / exp(4e-3 * sqrt(prec)) &&
        x < 0.8 * sqrt(prec) + 0.65e-14 * pow(prec, 3) + 1.5e-33*pow(prec, 6))
    {
        if (arb_hypgeom_erf_bb(res, z, 1, prec))
            return;
    }

    if (arf_cmpabs_2exp_si(arb_midref(z), 30) > 0)
    {
        log_z = ARF_EXP(arb_midref(z)) * LOG2;
    }
    else
    {
        abs_z2 = x * x;
        log_z = 0.5 * log(abs_z2);
    }

    /* Can we use the asymptotic expansion? */
    /* N! / z^(2N) < 2^p */
    if (x > 2.0)
    {
        double log2_err, err_prev;
        prec2 = prec + 5 + FLINT_BIT_COUNT(prec);

        err_prev = 0.0;
        for (N = 1; ; N++)
        {
            log2_err = -(2 * N) * log_z + N * (log(N) - 1.0);
            log2_err *= INV_LOG2;

            if (log2_err > err_prev)
                break;

            if (log2_err < -prec - 5)
            {
                arb_hypgeom_erf_asymp(res, z, N, 1, prec, prec2);
                return;
            }

            err_prev = log2_err;
        }
    }

    /* Compute via 1F1 - with cancellation */
    if (x >= 1.0)
        wp = prec + (x * x + log_z) * INV_LOG2;
    else
        wp = prec - log_z * INV_LOG2;

    wp = wp + 10 + FLINT_BIT_COUNT(prec);

    arb_hypgeom_erf_1f1(res, z, wp, wp);
    arb_sub_ui(res, res, 1, prec);
    arb_neg(res, res);
}

