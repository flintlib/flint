/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "bernoulli.h"

/* todo: move/cleanup helper functions */

void
acb_gamma_bound_phase(mag_t bound, const acb_t z)
{
    arf_t x, y, t, u;
    int xsign;
    slong prec;

    arf_init(x);
    arf_init(y);
    arf_init(t);
    arf_init(u);

    prec = MAG_BITS;

    /* first compute x, y such that |arg(z)| <= arg(x+yi) */

    /* argument increases with smaller real parts */
    arf_set_mag(x, arb_radref(acb_realref(z)));
    arf_sub(x, arb_midref(acb_realref(z)), x, prec, ARF_RND_FLOOR);

    xsign = arf_sgn(x);

    if (xsign >= 0)  /* argument increases away from the real axis */
        arb_get_abs_ubound_arf(y, acb_imagref(z), prec);
    else  /* argument increases closer to the real axis */
        arb_get_abs_lbound_arf(y, acb_imagref(z), prec);

    if (arf_is_zero(y))
    {
        if (xsign > 0)
            mag_one(bound);
        else
            mag_inf(bound);
    }
    else
    {
        if (xsign >= 0)
        {
            /* compute upper bound for t = y / (sqrt(x^2 + y^2) + x) */
            arf_mul(t, x, x, prec, ARF_RND_DOWN);
            arf_mul(u, y, y, prec, ARF_RND_DOWN);
            arf_add(t, t, u, prec, ARF_RND_DOWN);
            arf_sqrt(t, t, prec, ARF_RND_DOWN);
            arf_add(t, t, x, prec, ARF_RND_DOWN);
            arf_div(t, y, t, prec, ARF_RND_UP);
        }
        else
        {
            /* compute upper bound for t = (sqrt(x^2 + y^2) - x) / y */
            arf_mul(t, x, x, prec, ARF_RND_UP);
            arf_mul(u, y, y, prec, ARF_RND_UP);
            arf_add(t, t, u, prec, ARF_RND_UP);
            arf_sqrt(t, t, prec, ARF_RND_UP);
            arf_sub(t, t, x, prec, ARF_RND_UP);
            arf_div(t, t, y, prec, ARF_RND_UP);
        }

        /* compute upper bound for sqrt(1 + t^2) */
        arf_mul(t, t, t, prec, ARF_RND_UP);
        arf_add_ui(t, t, 1, prec, ARF_RND_UP);
        arf_sqrt(t, t, prec, ARF_RND_UP);

        arf_get_mag(bound, t);
    }

    arf_clear(x);
    arf_clear(y);
    arf_clear(t);
    arf_clear(u);
}

/*
  2 |B_{2n}| G(2n+k-1) / (G(k+1) G(2n+1)) |z| (T |z|^{-1})^(2n+k)
  TODO: CHECK n >= 1 ?
*/
void
acb_gamma_stirling_bound(mag_ptr err, const acb_t z, slong k0, slong knum, slong n)
{
    mag_t c, t, u, v;
    slong i, k;

    if (arb_contains_zero(acb_imagref(z)) &&
        arb_contains_nonpositive(acb_realref(z)))
    {
        for (i = 0; i < knum; i++)
            mag_inf(err + i);
        return;
    }

    mag_init(c);
    mag_init(t);
    mag_init(u);
    mag_init(v);

    /* t = lower bound for |z| */
    acb_get_mag_lower(t, z);
    /* v = upper bound for |z| */
    acb_get_mag(v, z);

    /* c = upper bound for 1/(cos(arg(z)/2) |z|) */
    acb_gamma_bound_phase(c, z);
    mag_div(c, c, t);

    /* numerator: 2 B_{2n} gamma(2n+k-1) |z| */
    mag_bernoulli_div_fac_ui(err, 2 * n);
    mag_mul_2exp_si(err, err, 1);
    mag_fac_ui(u, 2 * n + k0 - 2);
    mag_mul(err, err, u);
    mag_mul(err, err, v);

    /* denominator gamma(k+1) gamma(2n+1) */
    mag_rfac_ui(t, k0);
    mag_mul(err, err, t);

    /* multiply by c^(2n+k) */
    mag_pow_ui(t, c, 2 * n + k0);
    mag_mul(err, err, t);

    for (i = 1; i < knum; i++)
    {
        /* recurrence factor: c * (2n+k-2) / k */
        k = k0 + i;
        mag_mul(err + i, err + i - 1, c);
        mag_mul_ui(err + i, err + i, 2 * n + k - 2);
        mag_div_ui(err + i, err + i, k);
    }

    mag_clear(c);
    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
}

void
arb_gamma_stirling_bound(mag_ptr err, const arb_t x, slong k0, slong knum, slong n)
{
    acb_t z;
    acb_init(z);
    acb_set_arb(z, x);
    acb_gamma_stirling_bound(err, z, k0, knum, n);
    acb_clear(z);
}

void
arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec)
{
    fmpz_t d;
    fmpz_init(d);

    BERNOULLI_ENSURE_CACHED(2 * k);

    arb_set_round_fmpz(b, fmpq_numref(bernoulli_cache + 2 * k), prec);

    if (digamma)
        fmpz_mul_ui(d, fmpq_denref(bernoulli_cache + 2 * k), 2 * k);
    else
        fmpz_mul2_uiui(d, fmpq_denref(bernoulli_cache + 2 * k), 2 * k, 2 * k - 1);

    arb_div_fmpz(b, b, d, prec);
    fmpz_clear(d);
}



void
arb_gamma_stirling_eval(arb_t s, const arb_t z, slong nterms, int digamma, slong prec)
{
    arb_t b, t, logz, zinv, zinv2;
    mag_t err;

    slong k, term_prec;
    double z_mag, term_mag;

    arb_init(b);
    arb_init(t);
    arb_init(logz);
    arb_init(zinv);
    arb_init(zinv2);

    arb_log(logz, z, prec);
    arb_inv(zinv, z, prec);

    nterms = FLINT_MAX(nterms, 1);

    arb_zero(s);

    if (nterms > 1)
    {
        arb_mul(zinv2, zinv, zinv, prec);

        z_mag = arf_get_d(arb_midref(logz), ARF_RND_UP) * 1.44269504088896;

        for (k = nterms - 1; k >= 1; k--)
        {
            term_mag = bernoulli_bound_2exp_si(2 * k);
            term_mag -= (2 * k - 1) * z_mag;
            term_prec = prec + term_mag;
            term_prec = FLINT_MIN(term_prec, prec);
            term_prec = FLINT_MAX(term_prec, 10);

            if (prec > 2000)
            {
                arb_set_round(t, zinv2, term_prec);
                arb_mul(s, s, t, term_prec);
            }
            else
                arb_mul(s, s, zinv2, term_prec);

            arb_gamma_stirling_coeff(b, k, digamma, term_prec);
            arb_add(s, s, b, term_prec);
        }

        if (digamma)
            arb_mul(s, s, zinv2, prec);
        else
            arb_mul(s, s, zinv, prec);
    }

    /* remainder bound */
    mag_init(err);
    arb_gamma_stirling_bound(err, z, digamma ? 1 : 0, 1, nterms);
    mag_add(arb_radref(s), arb_radref(s), err);
    mag_clear(err);

    if (digamma)
    {
        arb_neg(s, s);
        arb_mul_2exp_si(zinv, zinv, -1);
        arb_sub(s, s, zinv, prec);
        arb_add(s, s, logz, prec);
    }
    else
    {
        /* (z-0.5)*log(z) - z + log(2*pi)/2 */
        arb_one(t);
        arb_mul_2exp_si(t, t, -1);
        arb_sub(t, z, t, prec);
        arb_mul(t, logz, t, prec);
        arb_add(s, s, t, prec);
        arb_sub(s, s, z, prec);
        arb_const_log_sqrt2pi(t, prec);
        arb_add(s, s, t, prec);
    }

    arb_clear(t);
    arb_clear(b);
    arb_clear(zinv);
    arb_clear(zinv2);
    arb_clear(logz);
}

void
arb_gamma_fmpq(arb_t y, const fmpq_t x, slong prec)
{
    arb_hypgeom_gamma_fmpq(y, x, prec);
}

void
arb_gamma_fmpz(arb_t y, const fmpz_t x, slong prec)
{
    arb_hypgeom_gamma_fmpz(y, x, prec);
}

void
arb_gamma(arb_t y, const arb_t x, slong prec)
{
    arb_hypgeom_gamma(y, x, prec);
}

void
arb_rgamma(arb_t y, const arb_t x, slong prec)
{
    arb_hypgeom_rgamma(y, x, prec);
}

void
arb_lgamma(arb_t y, const arb_t x, slong prec)
{
    arb_hypgeom_lgamma(y, x, prec);
}

