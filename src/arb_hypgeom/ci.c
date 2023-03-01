/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "arb_hypgeom.h"

#define LOG2 0.69314718055994531
#define INV_LOG2 1.4426950408889634074
#define EXP1 2.7182818284590452

double arf_get_d_log2_abs_approx_clamped(const arf_t x);

void
_arb_hypgeom_ci_asymp(arb_t res, const arb_t z, slong N, slong prec)
{
    arb_t s, c, sz, cz, u;
    fmpq a[1];
    slong wp;
    mag_t err, t;

    N = FLINT_MAX(N, 1);

    arb_init(s);
    arb_init(c);
    arb_init(sz);
    arb_init(cz);
    arb_init(u);
    mag_init(err);
    mag_init(t);

    /* Error is bounded by first omitted term, N! / z^N */
    arb_get_mag_lower(err, z);

    mag_inv(err, err);
    mag_pow_ui(err, err, N);
    mag_fac_ui(t, N);
    mag_mul(err, err, t);

    wp = prec * 1.001 + 5;

    arb_set(u, z);
    *fmpq_numref(&a[0]) = 1;
    *fmpq_denref(&a[0]) = 1;
    arb_hypgeom_sum_fmpq_imag_arb(c, s, a, 1, NULL, 0, u, 1, N, wp);
    arb_add_error_mag(c, err);
    arb_add_error_mag(s, err);

    arb_sin_cos(sz, cz, z, wp);
    arb_mul(c, c, sz, wp);
    arb_submul(c, s, cz, wp);
    arb_div(res, c, z, prec);

    arb_clear(s);
    arb_clear(c);
    arb_clear(sz);
    arb_clear(cz);
    arb_clear(u);
    mag_clear(err);
    mag_clear(t);
}

void
_arb_hypgeom_ci_2f3(arb_t res, const arb_t z, slong N, slong wp, slong prec)
{
    mag_t err, t;
    arb_t s, u;
    fmpq a[1];
    fmpq b[3];

    N = FLINT_MAX(N, 1);

    mag_init(err);
    mag_init(t);
    arb_init(s);
    arb_init(u);

    arb_sqr(u, z, wp);
    arb_mul_2exp_si(u, u, -2);
    arb_neg(u, u);

    *fmpq_numref(&a[0]) = 1;
    *fmpq_denref(&a[0]) = 1;
    *fmpq_numref(&b[0]) = 2;
    *fmpq_denref(&b[0]) = 1;
    *fmpq_numref(&b[1]) = 2;
    *fmpq_denref(&b[1]) = 1;
    *fmpq_numref(&b[2]) = 3;
    *fmpq_denref(&b[2]) = 2;

    /* Terms are bounded by u^N / (4 (N!)^2) */
    arb_get_mag(err, u);

    /* u^N */
    mag_set(t, err);
    mag_pow_ui(t, t, N);

    /* geometric factor for u/N^2 */
    mag_div_ui(err, err, N);
    mag_div_ui(err, err, N);
    mag_geom_series(err, err, 0);
    mag_mul(t, t, err);

    /* 1/(N!)^2 */
    mag_rfac_ui(err, N);
    mag_mul(err, err, err);
    mag_mul(err, err, t);

    /* 1/4 */
    mag_mul_2exp_si(err, err, -2);

    arb_hypgeom_sum_fmpq_arb(s, a, 1, b, 3, u, 0, N, wp);

    arb_add_error_mag(s, err);

    arb_mul(s, s, u, wp);

    arb_log(u, z, wp);
    arb_add(s, s, u, wp);

    arb_const_euler(u, wp);
    arb_add(res, s, u, prec);

    mag_clear(err);
    mag_clear(t);
    arb_clear(u);
    arb_clear(s);
}

void
arb_hypgeom_ci(arb_t res, const arb_t z, slong prec)
{
    slong wp, N, acc;
    double dz, du;

    if (!arb_is_positive(z) || !arb_is_finite(z))
    {
        arb_indeterminate(res);
        return;
    }

    if (ARF_IS_LAGOM(arb_midref(z)))
    {
        acc = arb_rel_accuracy_bits(z);
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        acc += FLINT_MAX(0, -ARF_EXP(arb_midref(z)));
        prec = FLINT_MIN(prec, acc + 32);
    }

    dz = fabs(arf_get_d(arb_midref(z), ARF_RND_DOWN));
    dz = FLINT_MIN(dz, 1e300);

    if (dz > 2.0)
    {
        double log2_err, err_prev, log2dz;

        log2dz = arf_get_d_log2_abs_approx_clamped(arb_midref(z));

        err_prev = 0.0;
        for (N = 1; N < 2 * prec; N++)
        {
            log2_err = ((N + 1.0) * (log(N + 1.0) - 1.0)) * INV_LOG2 - N * log2dz;

            if (log2_err > err_prev)
                break;

            if (log2_err < -prec - 2)
            {
                _arb_hypgeom_ci_asymp(res, z, N, prec);
                return;
            }

            err_prev = log2_err;
        }
    }

    if (arf_cmpabs_2exp_si(arb_midref(z), -30) < 0)
    {
        N = -arf_abs_bound_lt_2exp_si(arb_midref(z));
        wp = prec * 1.001 + 10;
        N = (prec + N - 1) / N;
    }
    else
    {
        du = 0.25 * dz * dz;
        wp = prec * 1.001 + 10;
        if (du > 1.0)
            wp += dz * 1.4426950408889634;

        N = (prec + 5) * LOG2 / (2 * d_lambertw((prec + 5) * LOG2 / (2 * EXP1 * sqrt(du)))) + 1;
    }

    if (arb_is_exact(z))
    {
        _arb_hypgeom_ci_2f3(res, z, N, wp, prec);
    }
    else
    {
        mag_t err;
        mag_init(err);

        /* |ci'(z)| = |cos(z)/z| <= 1/z */
        arb_get_mag_lower(err, z);

        if (mag_cmp_2exp_si(err, 0) >= 0 || 1)
        {
            arb_t zmid;
            arb_init(zmid);
            arb_get_mid_arb(zmid, z);
            mag_inv(err, err);
            mag_mul(err, err, arb_radref(z));
            _arb_hypgeom_ci_2f3(res, zmid, N, wp, prec);
            arb_add_error_mag(res, err);
            arb_clear(zmid);
        }
        else
        {
            _arb_hypgeom_ci_2f3(res, z, N, wp, prec);
        }

        mag_clear(err);
    }
}
