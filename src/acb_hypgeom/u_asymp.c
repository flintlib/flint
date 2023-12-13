/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

slong
acb_hypgeom_pfq_choose_n_max(acb_srcptr a, slong p,
                         acb_srcptr b, slong q, const acb_t z,
                         slong prec, slong n_max);

int
acb_hypgeom_u_asymp_determine_region(const mag_t r,
    const mag_t zlo, const acb_t z)
{
    int R;
    mag_t u, zre, zim;

    mag_init(u);
    mag_init(zre);
    mag_init(zim);

    /* lower bounds for |re(z)|, |im(z)|, |z| */
    arb_get_mag_lower(zre, acb_realref(z));
    arb_get_mag_lower(zim, acb_imagref(z));

    R = 0;

    if (mag_cmp(zlo, r) >= 0)
    {
        int znonneg = arb_is_nonnegative(acb_realref(z));

        if (znonneg && mag_cmp(zre, r) >= 0)
        {
            R = 1;
        }
        else if (mag_cmp(zim, r) >= 0 || znonneg)
        {
            R = 2;
        }
        else
        {
            mag_mul_2exp_si(u, r, 1);
            if (mag_cmp(zlo, u) >= 0)
                R = 3;
        }
    }

    mag_clear(u);
    mag_clear(zre);
    mag_clear(zim);

    return R;
}

/* computes the factors that are independent of n (all are upper bounds) */
void
acb_hypgeom_u_asymp_bound_factors(int * R, mag_t alpha,
    mag_t nu, mag_t sigma, mag_t rho, mag_t zinv,
    const acb_t a, const acb_t b, const acb_t z)
{
    mag_t r, u, zlo, sigma_prime;
    acb_t t;

    mag_init(r);
    mag_init(u);
    mag_init(zlo);
    mag_init(sigma_prime);
    acb_init(t);

    /* upper bound for 1/|z| */
    acb_get_mag_lower(zlo, z);
    mag_inv(zinv, zlo);

    /* upper bound for r = |b - 2a| */
    acb_mul_2exp_si(t, a, 1);
    acb_sub(t, b, t, MAG_BITS);
    acb_get_mag(r, t);

    /* determine region */
    *R = acb_hypgeom_u_asymp_determine_region(r, zlo, z);

    if (*R == 0)
    {
        mag_inf(alpha);
        mag_inf(nu);
        mag_inf(sigma);
        mag_inf(rho);
    }
    else
    {
        /* sigma = |(b-2a)/z| */
        mag_mul(sigma, r, zinv);

        /* nu = (1/2 + 1/2 sqrt(1-4 sigma^2))^(-1/2) <= 1 + 2 sigma^2 */
        if (mag_cmp_2exp_si(sigma, -1) <= 0)
        {
            mag_mul(nu, sigma, sigma);
            mag_mul_2exp_si(nu, nu, 1);
            mag_one(u);
            mag_add(nu, nu, u);
        }
        else
        {
            mag_inf(nu);
        }

        /* modified sigma for alpha, beta, rho when in R3 */
        if (*R == 3)
            mag_mul(sigma_prime, sigma, nu);
        else
            mag_set(sigma_prime, sigma);

        /* alpha = 1/(1-sigma') */
        mag_one(alpha);
        mag_sub_lower(alpha, alpha, sigma_prime);
        mag_one(u);
        mag_div(alpha, u, alpha);

        /* rho = |2a^2-2ab+b|/2 + sigma'*(1+sigma'/4)/(1-sigma')^2 */
        mag_mul_2exp_si(rho, sigma_prime, -2);
        mag_one(u);
        mag_add(rho, rho, u);
        mag_mul(rho, rho, sigma_prime);
        mag_mul(rho, rho, alpha);
        mag_mul(rho, rho, alpha);
        acb_sub(t, a, b, MAG_BITS);
        acb_mul(t, t, a, MAG_BITS);
        acb_mul_2exp_si(t, t, 1);
        acb_add(t, t, b, MAG_BITS);
        acb_get_mag(u, t);
        mag_mul_2exp_si(u, u, -1);
        mag_add(rho, rho, u);
    }

    mag_clear(r);
    mag_clear(u);
    mag_clear(zlo);
    mag_clear(sigma_prime);
    acb_clear(t);
}

void
acb_hypgeom_mag_chi(mag_t chi, ulong n)
{
    mag_t p, q;
    ulong k;

    mag_init(p);
    mag_init(q);

    if (n % 2 == 0)
    {
        mag_one(p);
        mag_one(q);
    }
    else
    {
        /* upper bound for pi/2 */
        mag_set_ui_2exp_si(p, 843314857, -28);
        mag_one(q);
    }

    for (k = n; k >= 2; k -= 2)
    {
        mag_mul_ui(p, p, k);
        mag_mul_ui_lower(q, q, k - 1);
    }

    mag_div(chi, p, q);

    mag_clear(p);
    mag_clear(q);
}

static void
acb_hypgeom_mag_Cn(mag_t Cn, int R, const mag_t nu, const mag_t sigma, ulong n)
{
    if (R == 1)
    {
        mag_one(Cn);
    }
    else
    {
        acb_hypgeom_mag_chi(Cn, n);

        if (R == 3)
        {
            mag_t tmp;
            mag_init(tmp);
            mag_mul(tmp, nu, nu);
            mag_mul(tmp, tmp, sigma);
            if (n != 1)
                mag_mul_ui(tmp, tmp, n);
            mag_add(Cn, Cn, tmp);
            mag_pow_ui(tmp, nu, n);
            mag_mul(Cn, Cn, tmp);
            mag_clear(tmp);
        }
    }
}

static int
acb_is_nonpositive_int(const acb_t x)
{
    return acb_is_int(x) && arf_sgn(arb_midref(acb_realref(x))) <= 0;
}

void acb_hypgeom_u_asymp(acb_t res, const acb_t a, const acb_t b,
    const acb_t z, slong n, slong prec)
{
    acb_struct aa[3];
    acb_t s, t, w, winv;
    int R, p, q, is_real, is_terminating;
    slong n_terminating;

    if (!acb_is_finite(a) || !acb_is_finite(b) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(aa);
    acb_init(aa + 1);
    acb_init(aa + 2);
    acb_init(s);
    acb_init(t);
    acb_init(w);
    acb_init(winv);

    is_terminating = 0;
    n_terminating = WORD_MAX;

    /* special case, for incomplete gamma
      [todo: also when they happen to be exact and with difference 1...] */
    if (a == b)
    {
        acb_set(aa, a);
        p = 1;
        q = 0;
    }
    else
    {
        acb_set(aa, a);
        acb_sub(aa + 1, a, b, prec);
        acb_add_ui(aa + 1, aa + 1, 1, prec);
        acb_one(aa + 2);
        p = 2;
        q = 1;
    }

    if (acb_is_nonpositive_int(aa))
    {
        is_terminating = 1;

        if (arf_cmpabs_ui(arb_midref(acb_realref(aa)), prec) < 0)
            n_terminating = 1 - arf_get_si(arb_midref(acb_realref(aa)), ARF_RND_DOWN);
    }

    if (p == 2 && acb_is_nonpositive_int(aa + 1))
    {
        is_terminating = 1;

        if (arf_cmpabs_ui(arb_midref(acb_realref(aa + 1)), n_terminating) < 0)
            n_terminating = 1 - arf_get_si(arb_midref(acb_realref(aa + 1)), ARF_RND_DOWN);
    }

    acb_neg(w, z);
    acb_inv(w, w, prec);
    acb_neg(winv, z);

    /* low degree polynomial -- no need to try to terminate sooner */
    if (is_terminating && n_terminating < 8)
    {
        acb_hypgeom_pfq_sum_invz(s, t, aa, p, aa + p, q, w, winv,
            n_terminating, prec);
        acb_set(res, s);
    }
    else
    {
        mag_t C1, Cn, alpha, nu, sigma, rho, zinv, tmp, err;

        mag_init(C1);
        mag_init(Cn);
        mag_init(alpha);
        mag_init(nu);
        mag_init(sigma);
        mag_init(rho);
        mag_init(zinv);
        mag_init(tmp);
        mag_init(err);

        acb_hypgeom_u_asymp_bound_factors(&R, alpha, nu,
            sigma, rho, zinv, a, b, z);

        is_real = acb_is_real(a) && acb_is_real(b) && acb_is_real(z) &&
            (is_terminating || arb_is_positive(acb_realref(z)));

        if (R == 0)
        {
            /* if R == 0, the error bound is infinite unless terminating */
            if (is_terminating && n_terminating < prec)
            {
                acb_hypgeom_pfq_sum_invz(s, t, aa, p, aa + p, q, w, winv,
                    n_terminating, prec);
                acb_set(res, s);
            }
            else
            {
                acb_indeterminate(res);
            }
        }
        else
        {
            /* C1 */
            acb_hypgeom_mag_Cn(C1, R, nu, sigma, 1);

            /* err = 2 * alpha * exp(...) */
            mag_mul(tmp, C1, rho);
            mag_mul(tmp, tmp, alpha);
            mag_mul(tmp, tmp, zinv);
            mag_mul_2exp_si(tmp, tmp, 1);
            mag_exp(err, tmp);
            mag_mul(err, err, alpha);
            mag_mul_2exp_si(err, err, 1);

            /* choose n automatically */
            if (n < 0)
            {
                slong moreprec;

                /* take err into account when finding truncation point */
                /* we should take Cn into account as well, but this depends
                   on n which is to be determined; it's easier to look
                   only at exp(...) which should be larger anyway */
                if (mag_cmp_2exp_si(err, 10 * prec) > 0)
                    moreprec = 10 * prec;
                else if (mag_cmp_2exp_si(err, 0) < 0)
                    moreprec = 0;
                else
                    moreprec = MAG_EXP(err);

                n = acb_hypgeom_pfq_choose_n_max(aa, p, aa + p, q, w,
                    prec + moreprec, FLINT_MIN(WORD_MAX / 2, 50 + 10 * prec));
            }

            acb_hypgeom_pfq_sum_invz(s, t, aa, p, aa + p, q, w, winv, n, prec);

            /* add error bound, if not terminating */
            if (!(is_terminating && n == n_terminating))
            {
                acb_hypgeom_mag_Cn(Cn, R, nu, sigma, n);
                mag_mul(err, err, Cn);

                /* nth term * factor */
                acb_get_mag(tmp, t);
                mag_mul(err, err, tmp);

                if (is_real)
                    arb_add_error_mag(acb_realref(s), err);
                else
                    acb_add_error_mag(s, err);
            }

            acb_set(res, s);
        }

        mag_clear(C1);
        mag_clear(Cn);
        mag_clear(alpha);
        mag_clear(nu);
        mag_clear(sigma);
        mag_clear(rho);
        mag_clear(zinv);
        mag_clear(tmp);
        mag_clear(err);
    }

    acb_clear(aa);
    acb_clear(aa + 1);
    acb_clear(aa + 2);
    acb_clear(s);
    acb_clear(t);
    acb_clear(w);
    acb_clear(winv);
}

