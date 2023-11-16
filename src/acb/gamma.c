/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bernoulli.h"
#include "acb.h"
#include "acb_hypgeom.h"

void
acb_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t z, int use_reflect, int digamma, slong prec);

void acb_gamma_stirling_bound(mag_ptr err, const acb_t z, slong k0, slong knum, slong n);

void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);


void
acb_gamma_stirling_eval(acb_t s, const acb_t z, slong nterms, int digamma, slong prec)
{
    acb_t t, logz, zinv, zinv2;
    arb_t b;
    mag_t err;

    slong k, term_prec;
    double z_mag, term_mag;

    acb_init(t);
    acb_init(logz);
    acb_init(zinv);
    acb_init(zinv2);
    arb_init(b);

    acb_log(logz, z, prec);
    acb_inv(zinv, z, prec);

    nterms = FLINT_MAX(nterms, 1);

    acb_zero(s);
    if (nterms > 1)
    {
        acb_mul(zinv2, zinv, zinv, prec);

        z_mag = arf_get_d(arb_midref(acb_realref(logz)), ARF_RND_UP) * 1.44269504088896;

        for (k = nterms - 1; k >= 1; k--)
        {
            term_mag = bernoulli_bound_2exp_si(2 * k);
            term_mag -= (2 * k - 1) * z_mag;
            term_prec = prec + term_mag;
            term_prec = FLINT_MIN(term_prec, prec);
            term_prec = FLINT_MAX(term_prec, 10);

            arb_gamma_stirling_coeff(b, k, digamma, term_prec);

            if (prec > 2000)
            {
                acb_set_round(t, zinv2, term_prec);
                acb_mul(s, s, t, term_prec);
            }
            else
                acb_mul(s, s, zinv2, term_prec);

            arb_add(acb_realref(s), acb_realref(s), b, term_prec);
        }

        if (digamma)
            acb_mul(s, s, zinv2, prec);
        else
            acb_mul(s, s, zinv, prec);
    }

    /* remainder bound */
    mag_init(err);
    acb_gamma_stirling_bound(err, z, digamma ? 1 : 0, 1, nterms);
    mag_add(arb_radref(acb_realref(s)), arb_radref(acb_realref(s)), err);
    mag_add(arb_radref(acb_imagref(s)), arb_radref(acb_imagref(s)), err);
    mag_clear(err);

    if (digamma)
    {
        acb_neg(s, s);
        acb_mul_2exp_si(zinv, zinv, -1);
        acb_sub(s, s, zinv, prec);
        acb_add(s, s, logz, prec);
    }
    else
    {
        /* (z-0.5)*log(z) - z + log(2*pi)/2 */
        arb_one(b);
        arb_mul_2exp_si(b, b, -1);
        arb_set(acb_imagref(t), acb_imagref(z));
        arb_sub(acb_realref(t), acb_realref(z), b, prec);
        acb_mul(t, logz, t, prec);
        acb_add(s, s, t, prec);
        acb_sub(s, s, z, prec);
        arb_const_log_sqrt2pi(b, prec);
        arb_add(acb_realref(s), acb_realref(s), b, prec);
    }

    acb_clear(t);
    acb_clear(logz);
    acb_clear(zinv);
    acb_clear(zinv2);
    arb_clear(b);
}

void
acb_gamma(acb_t y, const acb_t x, slong prec)
{
    acb_hypgeom_gamma(y, x, prec);
}

void
acb_rgamma(acb_t y, const acb_t x, slong prec)
{
    acb_hypgeom_rgamma(y, x, prec);
}

void
acb_lgamma(acb_t y, const acb_t x, slong prec)
{
    acb_hypgeom_lgamma(y, x, prec);
}
