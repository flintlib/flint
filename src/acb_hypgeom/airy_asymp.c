/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb.h"
#include "acb_hypgeom.h"

void acb_hypgeom_mag_chi(mag_t chi, ulong n);

static int
arg_lt_2pi3(const acb_t z, const acb_t zeta)
{
    if (arb_is_nonnegative(acb_realref(z)))
        return 1;

    if (arb_is_positive(acb_imagref(z)) &&
        arb_is_positive(acb_imagref(zeta)))
        return 1;

    if (arb_is_negative(acb_imagref(z)) &&
        arb_is_negative(acb_imagref(zeta)))
        return 1;

    return 0;
}

/* assuming a >= b >= c >= d, e >= f */
void
_acb_mul4div2_ui(acb_t x, ulong a, ulong b, ulong c, ulong d, ulong e, ulong f, slong prec)
{
    if (a < (UWORD(1) << (FLINT_BITS / 4)))
    {
        acb_mul_ui(x, x, a * b * c * d, prec);
    }
    else if (a < (UWORD(1) << (FLINT_BITS / 2)))
    {
        acb_mul_ui(x, x, a * b, prec);
        acb_mul_ui(x, x, c * d, prec);
    }
    else
    {
        acb_mul_ui(x, x, a, prec);
        acb_mul_ui(x, x, b, prec);
        acb_mul_ui(x, x, c, prec);
        acb_mul_ui(x, x, d, prec);
    }

    if (e < (UWORD(1) << (FLINT_BITS / 2)))
    {
        acb_div_ui(x, x, e * f, prec);
    }
    else
    {
        acb_div_ui(x, x, e, prec);
        acb_div_ui(x, x, f, prec);
    }
}

void
acb_hypgeom_airy_asymp_sum(acb_t s0even, acb_t s0odd,
        acb_t s1even, acb_t s1odd,
        mag_t t0n, mag_t t1n,
        const acb_t z, const acb_t z2, slong n, slong prec)
{
    slong m, k, j;
    acb_ptr z2pow;

    acb_get_mag(t0n, z);
    mag_mul_ui(t0n, t0n, 72);
    mag_pow_ui(t0n, t0n, n);
    mag_one(t1n);

    for (k = 1; k <= n; k++)
    {
        mag_mul_ui(t0n, t0n, 6 * k - 1);
        mag_mul_ui(t0n, t0n, 6 * k - 5);
        mag_mul_ui_lower(t1n, t1n, 72 * k);
    }

    mag_div(t0n, t0n, t1n);
    mag_mul_ui(t1n, t0n, 6 * n + 1);
    mag_div_ui(t1n, t1n, 6 * n - 1);

    m = n_sqrt(n / 2);
    m = FLINT_MAX(m, 1);

    z2pow = _acb_vec_init(m + 1);
    _acb_vec_set_powers(z2pow, z2, m + 1, prec);

    if (s0even != NULL)
    {
        acb_zero(s0even);

        for (k = n + (n % 2); k >= 0; k -= 2)
        {
            j = (k / 2) % m;

            if (k < n)
                acb_add(s0even, s0even, z2pow + j, prec);

            if (k > 0)
            {
                _acb_mul4div2_ui(s0even, 6*k-1, 6*k-5, 6*k-7, 6*k-11, k, k-1, prec);
                if (j == 0 && k < n)
                    acb_mul(s0even, s0even, z2pow + m, prec);
            }
        }
    }

    if (s0odd != NULL)
    {
        acb_zero(s0odd);

        for (k = n + 1 + (n % 2); k >= 1; k -= 2)
        {
            j = ((k - 1) / 2) % m;

            if (k < n)
                acb_add(s0odd, s0odd, z2pow + j, prec);

            if (k > 1)
            {
                _acb_mul4div2_ui(s0odd, 6*k-1, 6*k-5, 6*k-7, 6*k-11, k, k-1, prec);
                if (j == 0 && k < n)
                    acb_mul(s0odd, s0odd, z2pow + m, prec);
            }
            else
            {
                acb_mul(s0odd, s0odd, z, prec);
                acb_mul_ui(s0odd, s0odd, 5, prec);
            }
        }
    }

    if (s1even != NULL)
    {
        acb_zero(s1even);

        for (k = n + (n % 2); k >= 0; k -= 2)
        {
            j = (k / 2) % m;

            if (k < n)
                acb_add(s1even, s1even, z2pow + j, prec);

            if (k > 0)
            {
                _acb_mul4div2_ui(s1even, 6*k+1, 6*k-5, 6*k-7, FLINT_ABS(6*k-13), k, k-1, prec);
                if (k == 2)
                    acb_neg(s1even, s1even);
                if (j == 0 && k < n)
                    acb_mul(s1even, s1even, z2pow + m, prec);
            }
        }
    }

    if (s1odd != NULL)
    {
        acb_zero(s1odd);

        for (k = n + 1 + (n % 2); k >= 1; k -= 2)
        {
            j = ((k - 1) / 2) % m;

            if (k < n)
                acb_add(s1odd, s1odd, z2pow + j, prec);

            if (k > 1)
            {
                _acb_mul4div2_ui(s1odd, 6*k+1, 6*k-5, 6*k-7, 6*k-13, k, k-1, prec);
                if (j == 0 && k < n)
                    acb_mul(s1odd, s1odd, z2pow + m, prec);
            }
            else
            {
                acb_mul(s1odd, s1odd, z, prec);
                acb_mul_si(s1odd, s1odd, -7, prec);
            }
        }
    }

    _acb_vec_clear(z2pow, m + 1);
}

void
acb_hypgeom_airy_asymp_bound_factor(mag_t bound, const acb_t z, const acb_t zeta, ulong n)
{
    mag_t t, u, v;

    mag_init(t);
    mag_init(u);
    mag_init(v);

    if (arb_is_nonnegative(acb_realref(z)) && arb_is_nonnegative(acb_realref(zeta)))
    {
        /* 2 exp(7 / (36 |zeta|)) */
        mag_set_ui_2exp_si(t, 50, -8);  /* bound for 7/36 */
        acb_get_mag_lower(u, zeta);
        mag_div(t, t, u);
        mag_exp(t, t);
        mag_mul_2exp_si(bound, t, 1);
    }
    else
    {
        /* 2 exp(7 pi / (72 |zeta|)) */
        mag_set_ui_2exp_si(t, 79, -8);  /* bound for 7 pi/72 */
        acb_get_mag_lower(u, zeta);
        mag_div(t, t, u);
        mag_exp(t, t);
        mag_mul_2exp_si(bound, t, 1);

        /* 4 exp(7 pi / (36 |re(zeta)|)) / |cos(arg(zeta))|^n */
        if (!arg_lt_2pi3(z, zeta))
        {
            arb_get_mag_lower(u, acb_realref(zeta));
            arb_get_mag(v, acb_imagref(zeta));

            /* 4 exp(7 pi / (36 u)) < exp(157/(256 u)) */
            mag_set_ui_2exp_si(t, 157, -8);
            mag_div(t, t, u);
            mag_exp(t, t);
            mag_mul_2exp_si(t, t, 2);

            /* |1/cos(arg(zeta))| = sqrt(1+(v/u)^2) */
            mag_div(v, v, u);
            mag_mul(v, v, v);
            mag_one(u);
            mag_add(v, v, u);
            mag_sqrt(v, v);
            mag_pow_ui(v, v, n);

            mag_mul(t, t, v);
            mag_max(bound, bound, t);
        }

        /* chi(n) */
        acb_hypgeom_mag_chi(t, n);
        mag_mul(bound, bound, t);
    }

    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
}

void
acb_hypgeom_airy_asymp(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, slong n, slong prec)
{
    acb_t t, u, w, z_root, zeta;
    acb_t s0even, s0odd, s1even, s1odd, E1, E2;
    mag_t err1, err2, erru, errv, errtmp;
    int want_d0, want_d1, is_real, upper;

    acb_init(t);
    acb_init(u);
    acb_init(w);
    acb_init(z_root);
    acb_init(zeta);
    acb_init(s0even);
    acb_init(s0odd);
    acb_init(s1even);
    acb_init(s1odd);
    acb_init(E1);
    acb_init(E2);

    mag_init(err1);
    mag_init(err2);
    mag_init(errtmp);
    mag_init(erru);
    mag_init(errv);

    is_real = acb_is_real(z);
    want_d0 = (ai != NULL) || (bi != NULL);
    want_d1 = (aip != NULL) || (bip != NULL);

    /* required for some of the error bounds to be valid */
    n = FLINT_MAX(n, 1);

    /* this is just a heuristic; the sectors are validated later */
    if (arf_sgn(arb_midref(acb_realref(z))) >= 0)
    {
        /* z_root = z^(1/4), zeta = 2 z^(3/2) / 3 */
        acb_sqrt(z_root, z, prec);
        acb_cube(zeta, z_root, prec);
        acb_sqrt(z_root, z_root, prec);
        acb_mul_2exp_si(zeta, zeta, 1);
        acb_div_ui(zeta, zeta, 3, prec);

        /* compute bound factor for z = t = z1, zeta = u = zeta1 */
        acb_hypgeom_airy_asymp_bound_factor(err1, z, zeta, n);

        /* this is just a heuristic; the sectors are validated later */
        upper = arf_sgn(arb_midref(acb_imagref(z))) >= 0;

        if (upper)
        {
            /* check -pi/3 < arg(z) < pi */
            if (arb_is_positive(acb_imagref(z)) ||
                (arb_is_positive(acb_realref(z)) &&
                arb_is_positive(acb_realref(zeta))))
            {
                arb_one(acb_realref(w));
                arb_sqrt_ui(acb_imagref(w), 3, prec);
                acb_mul_2exp_si(w, w, -1);
                acb_neg(w, w);
                acb_mul(t, w, z, prec);
                acb_neg(u, zeta);
                acb_hypgeom_airy_asymp_bound_factor(err2, t, u, n);
            }
            else
            {
                mag_inf(err2);
            }
        }
        else
        {
            /* check -pi < arg(z) < pi/3 */
            if (arb_is_negative(acb_imagref(z)) ||
                (arb_is_positive(acb_realref(z)) &&
                arb_is_positive(acb_realref(zeta))))
            {
                arb_set_si(acb_realref(w), -1);
                arb_sqrt_ui(acb_imagref(w), 3, prec);
                acb_mul_2exp_si(w, w, -1);
                acb_mul(t, w, z, prec);
                acb_neg(u, zeta);
                acb_hypgeom_airy_asymp_bound_factor(err2, t, u, n);
            }
            else
            {
                mag_inf(err2);
            }
        }

        acb_mul_ui(t, zeta, 72, prec);
        acb_inv(t, t, prec);
        acb_mul(u, t, t, prec);

        acb_hypgeom_airy_asymp_sum(want_d0 ? s0even : NULL,
                                   want_d0 ? s0odd : NULL,
                                   want_d1 ? s1even : NULL,
                                   want_d1 ? s1odd : NULL,
            erru, errv, t, u, n, prec);

        /* exponentials */
        acb_exp_invexp(E2, E1, zeta, prec);

        /* common prefactor */
        acb_const_pi(w, prec);
        acb_sqrt(w, w, prec);
        acb_mul_2exp_si(w, w, 1);

        if (ai != NULL || bi != NULL)
        {
            /* t = sum (-1)^k u(k) / (zeta)^k + error */
            acb_sub(t, s0even, s0odd, prec);
            mag_mul(errtmp, err1, erru);
            acb_add_error_mag(t, errtmp);

            /* u = sum (-1)^k u(k) / (-zeta)^k = sum u(k) / (zeta)^k + error */
            acb_add(u, s0even, s0odd, prec);
            mag_mul(errtmp, err2, erru);
            acb_add_error_mag(u, errtmp);

            if (ai != NULL)
            {
                acb_mul(ai, t, E1, prec);
            }

            if (bi != NULL)
            {
                acb_mul(t, t, E1, prec);
                acb_mul(u, u, E2, prec);
                acb_mul_2exp_si(u, u, 1);

                if (upper)
                    acb_mul_onei(t, t);
                else
                    acb_div_onei(t, t);

                acb_add(bi, t, u, prec);
            }

            acb_mul(t, w, z_root, prec);
            acb_inv(t, t, prec);

            if (ai != NULL) acb_mul(ai, ai, t, prec);
            if (bi != NULL) acb_mul(bi, bi, t, prec);
        }

        if (aip != NULL || bip != NULL)
        {
            /* t = sum (-1)^k v(k) / (zeta)^k + error */
            acb_sub(t, s1even, s1odd, prec);
            mag_mul(errtmp, err1, errv);
            acb_add_error_mag(t, errtmp);

            /* u = sum (-1)^k v(k) / (-zeta)^k = sum v(k) / (zeta)^k + error */
            acb_add(u, s1even, s1odd, prec);
            mag_mul(errtmp, err2, errv);
            acb_add_error_mag(u, errtmp);

            if (aip != NULL)
            {
                acb_mul(aip, t, E1, prec);
                acb_neg(aip, aip);
            }

            if (bip != NULL)
            {
                acb_mul(t, t, E1, prec);
                acb_mul(u, u, E2, prec);
                acb_mul_2exp_si(u, u, 1);

                if (upper)
                    acb_div_onei(t, t);
                else
                    acb_mul_onei(t, t);

                acb_add(bip, t, u, prec);
            }

            acb_inv(t, w, prec);
            acb_mul(t, t, z_root, prec);

            if (aip != NULL) acb_mul(aip, aip, t, prec);
            if (bip != NULL) acb_mul(bip, bip, t, prec);
        }
    }
    else
    {
        /* z_root = (-z)^(1/4), zeta = 2 (-z)^(3/2) / 3 */
        acb_neg(t, z);
        acb_sqrt(z_root, t, prec);
        acb_cube(zeta, z_root, prec);
        acb_sqrt(z_root, z_root, prec);
        acb_mul_2exp_si(zeta, zeta, 1);
        acb_div_ui(zeta, zeta, 3, prec);

        if (arg_lt_2pi3(t, zeta))
        {
            /* compute bound factor for i zeta */
            arb_one(acb_realref(w));
            arb_sqrt_ui(acb_imagref(w), 3, prec);
            acb_mul_2exp_si(w, w, -1);
            acb_mul(t, z, w, prec);
            acb_neg(t, t);
            acb_mul_onei(u, zeta);
            acb_hypgeom_airy_asymp_bound_factor(err1, t, u, n);

            /* compute bound factor for -i zeta */
            acb_conj(w, w);
            acb_mul(t, z, w, prec);
            acb_neg(t, t);
            acb_div_onei(u, zeta);
            acb_hypgeom_airy_asymp_bound_factor(err2, t, u, n);
        }
        else
        {
            mag_inf(err1);
            mag_inf(err2);
        }

        /* t = 1/(72 i zeta), u = (1/(72 i zeta))^2 */
        acb_mul_onei(t, zeta);
        acb_mul_ui(t, t, 72, prec);
        acb_inv(t, t, prec);
        acb_mul(u, t, t, prec);

        acb_hypgeom_airy_asymp_sum(want_d0 ? s0even : NULL,
                                   want_d0 ? s0odd : NULL,
                                   want_d1 ? s1even : NULL,
                                   want_d1 ? s1odd : NULL,
            erru, errv, t, u, n, prec);

        /* exponentials */
        acb_const_pi(t, prec);
        acb_mul_2exp_si(t, t, -2);
        acb_sub(t, zeta, t, prec);
        acb_mul_onei(t, t);
        acb_exp_invexp(E2, E1, t, prec);

        /* common prefactor */
        acb_const_pi(w, prec);
        acb_sqrt(w, w, prec);
        acb_mul_2exp_si(w, w, 1);

        if (ai != NULL || bi != NULL)
        {
            /* t = sum (-1)^k u(k) / (i zeta)^k + error */
            acb_sub(t, s0even, s0odd, prec);
            mag_mul(errtmp, err1, erru);
            acb_add_error_mag(t, errtmp);

            /* w = sum (-1)^k u(k) / (-i zeta)^k = sum u(k) / (i zeta)^k + error */
            acb_add(u, s0even, s0odd, prec);
            mag_mul(errtmp, err2, erru);
            acb_add_error_mag(u, errtmp);

            if (ai != NULL)
            {
                acb_mul(ai, t, E1, prec);
                acb_addmul(ai, u, E2, prec);
            }

            if (bi != NULL)
            {
                /* (-i E1) t + (+i E2) u */
                acb_mul(bi, u, E2, prec);
                acb_submul(bi, t, E1, prec);
                acb_mul_onei(bi, bi);
            }

            acb_mul(t, w, z_root, prec);
            acb_inv(t, t, prec);

            if (ai != NULL) acb_mul(ai, ai, t, prec);
            if (bi != NULL) acb_mul(bi, bi, t, prec);
        }

        if (aip != NULL || bip != NULL)
        {
            /* t = sum (-1)^k v(k) / (i zeta)^k + error */
            acb_sub(t, s1even, s1odd, prec);
            mag_mul(errtmp, err1, errv);
            acb_add_error_mag(t, errtmp);

            /* u = sum (-1)^k v(k) / (-i zeta)^k = sum v(k) / (i zeta)^k + error */
            acb_add(u, s1even, s1odd, prec);
            mag_mul(errtmp, err2, errv);
            acb_add_error_mag(u, errtmp);

            if (bip != NULL)
            {
                acb_mul(bip, t, E1, prec);
                acb_addmul(bip, u, E2, prec);
            }

            if (aip != NULL)
            {
                /* i E1 t - i E2 u */
                acb_mul(aip, t, E1, prec);
                acb_submul(aip, u, E2, prec);
                acb_mul_onei(aip, aip);
            }

            acb_inv(t, w, prec);
            acb_mul(t, t, z_root, prec);

            if (aip != NULL) acb_mul(aip, aip, t, prec);
            if (bip != NULL) acb_mul(bip, bip, t, prec);
        }
    }

    if (is_real)
    {
        if (ai != NULL) arb_zero(acb_imagref(ai));
        if (aip != NULL) arb_zero(acb_imagref(aip));
        if (bi != NULL) arb_zero(acb_imagref(bi));
        if (bip != NULL) arb_zero(acb_imagref(bip));
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(w);
    acb_clear(z_root);
    acb_clear(zeta);
    acb_clear(s0even);
    acb_clear(s0odd);
    acb_clear(s1even);
    acb_clear(s1odd);
    acb_clear(E1);
    acb_clear(E2);

    mag_clear(err1);
    mag_clear(err2);
    mag_clear(errtmp);
    mag_clear(erru);
    mag_clear(errv);
}

