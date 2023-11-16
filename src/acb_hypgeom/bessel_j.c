/*
    Copyright (C) 2014-2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

/* assumes no aliasing */
/* (+/- iz)^(-1/2-v) * z^v * exp(+/- iz) */
void
acb_hypgeom_bessel_j_asymp_prefactors_fallback(acb_t Ap, acb_t Am, acb_t C,
    const acb_t nu, const acb_t z, slong prec)
{
    acb_t t, u, v;

    acb_init(t);
    acb_init(u);
    acb_init(v);

    /* v = -1/2-nu */
    acb_one(v);
    acb_mul_2exp_si(v, v, -1);
    acb_add(v, v, nu, prec);
    acb_neg(v, v);

    acb_mul_onei(t, z);  /* t = iz */
    acb_neg(u, t);       /* u = -iz */

    /* Ap, Am = (+/- iz)^(-1/2-nu) */
    acb_pow(Ap, t, v, prec);
    acb_pow(Am, u, v, prec);

    /* Ap, Am *= exp(+/- iz) */
    acb_exp_invexp(u, v, t, prec);
    acb_mul(Ap, Ap, u, prec);
    acb_mul(Am, Am, v, prec);

    /* z^nu */
    acb_pow(t, z, nu, prec);
    acb_mul(Ap, Ap, t, prec);
    acb_mul(Am, Am, t, prec);

    /* (2 pi)^(-1/2) */
    acb_const_pi(C, prec);
    acb_mul_2exp_si(C, C, 1);
    acb_rsqrt(C, C, prec);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_hypgeom_bessel_j_asymp_prefactors(acb_t Ap, acb_t Am, acb_t C,
    const acb_t nu, const acb_t z, slong prec)
{
    if (arb_is_positive(acb_realref(z)))
    {
        acb_t t, u;

        acb_init(t);
        acb_init(u);

        /* -(2nu+1)/4 * pi + z */
        acb_mul_2exp_si(t, nu, 1);
        acb_add_ui(t, t, 1, prec);
        acb_mul_2exp_si(t, t, -2);
        acb_neg(t, t);
        acb_const_pi(u, prec);
        acb_mul(t, t, u, prec);
        acb_add(t, t, z, prec);
        acb_mul_onei(t, t);
        acb_exp_invexp(Ap, Am, t, prec);

        /* (2 pi z)^(-1/2) */
        acb_const_pi(C, prec);
        acb_mul_2exp_si(C, C, 1);
        acb_mul(C, C, z, prec);
        acb_rsqrt(C, C, prec);

        acb_clear(t);
        acb_clear(u);
        return;
    }

    acb_hypgeom_bessel_j_asymp_prefactors_fallback(Ap, Am, C, nu, z, prec);
}

void
acb_hypgeom_bessel_j_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    acb_t A1, A2, C, U1, U2, s, t, u;
    int is_real, is_imag;

    /* zero at -inf and +inf when nu is finite */
    if (acb_is_finite(nu) && !acb_is_finite(z) &&
        acb_is_real(z) && !acb_contains_zero(z))
    {
        acb_zero(res);
        return;
    }

    acb_init(A1);
    acb_init(A2);
    acb_init(C);
    acb_init(U1);
    acb_init(U2);
    acb_init(s);
    acb_init(t);
    acb_init(u);

    is_imag = 0;
    is_real = acb_is_real(nu) && acb_is_real(z)
        && (acb_is_int(nu) || arb_is_positive(acb_realref(z)));

    if (!is_real && arb_is_zero(acb_realref(z)) && acb_is_int(nu))
    {
        acb_mul_2exp_si(t, nu, -1);

        if (acb_is_int(t))
            is_real = 1;
        else
            is_imag = 1;
    }

    acb_hypgeom_bessel_j_asymp_prefactors(A1, A2, C, nu, z, prec);

    /* todo: if Ap ~ 2^a and Am = 2^b and U1 ~ U2 ~ 1, change precision? */

    if (!acb_is_finite(A1) || !acb_is_finite(A2) || !acb_is_finite(C))
    {
        acb_indeterminate(res);
    }
    else
    {
        /* s = 1/2 + nu */
        acb_one(s);
        acb_mul_2exp_si(s, s, -1);
        acb_add(s, s, nu, prec);

        /* t = 1 + 2 nu */
        acb_mul_2exp_si(t, nu, 1);
        acb_add_ui(t, t, 1, prec);

        acb_mul_onei(u, z);
        acb_mul_2exp_si(u, u, 1);
        acb_hypgeom_u_asymp(U2, s, t, u, -1, prec);
        acb_neg(u, u);
        acb_hypgeom_u_asymp(U1, s, t, u, -1, prec);

        acb_mul(res, A1, U1, prec);
        acb_addmul(res, A2, U2, prec);
        acb_mul(res, res, C, prec);

        if (is_real)
            arb_zero(acb_imagref(res));
        if (is_imag)
            arb_zero(acb_realref(res));
    }

    acb_clear(A1);
    acb_clear(A2);
    acb_clear(C);
    acb_clear(U1);
    acb_clear(U2);
    acb_clear(s);
    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    acb_struct b[2];
    acb_t w, c, t;

    if (acb_is_int(nu) && arb_is_negative(acb_realref(nu)))
    {
        acb_init(t);
        acb_neg(t, nu);

        acb_hypgeom_bessel_j_0f1(res, t, z, prec);

        acb_mul_2exp_si(t, t, -1);
        if (!acb_is_int(t))
            acb_neg(res, res);

        acb_clear(t);
        return;
    }

    acb_init(b + 0);
    acb_init(b + 1);
    acb_init(w);
    acb_init(c);
    acb_init(t);

    acb_add_ui(b + 0, nu, 1, prec);
    acb_one(b + 1);

    /* (z/2)^nu / gamma(nu+1) */
    acb_mul_2exp_si(c, z, -1);
    acb_pow(c, c, nu, prec);
    acb_rgamma(t, b + 0, prec);
    acb_mul(c, t, c, prec);

    /* -z^2/4 */
    acb_mul(w, z, z, prec);
    acb_mul_2exp_si(w, w, -2);
    acb_neg(w, w);

    acb_hypgeom_pfq_direct(t, NULL, 0, b, 2, w, -1, prec);

    acb_mul(res, t, c, prec);

    acb_clear(b + 0);
    acb_clear(b + 1);
    acb_clear(w);
    acb_clear(c);
    acb_clear(t);
}

/*

The asymptotic series can be used roughly when

[(1+log(2))/log(2) = 2.44269504088896] * z > p

We are a bit more conservative and use the factor 2.
*/

void
acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    mag_t zmag;

    mag_init(zmag);
    acb_get_mag(zmag, z);

    if (mag_cmp_2exp_si(zmag, 4) < 0 ||
        (mag_cmp_2exp_si(zmag, 64) < 0 && 2 * mag_get_d(zmag) < prec))
        acb_hypgeom_bessel_j_0f1(res, nu, z, prec);
    else
        acb_hypgeom_bessel_j_asymp(res, nu, z, prec);

    mag_clear(zmag);
}

