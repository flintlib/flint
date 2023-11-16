/*
    Copyright (C) 2014-2015, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_bessel_i_asymp_prefactors(acb_t A, acb_t B, acb_t C,
    const acb_t nu, const acb_t z, int scaled, slong prec)
{
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    /* C = (2 pi z)^(-1/2) */
    acb_const_pi(C, prec);
    acb_mul_2exp_si(C, C, 1);
    acb_mul(C, C, z, prec);
    acb_rsqrt(C, C, prec);

    if (arb_is_positive(acb_imagref(z)) ||
        (arb_is_zero(acb_imagref(z)) && arb_is_negative(acb_realref(z))))
    {
        acb_exp_pi_i(t, nu, prec);
        acb_mul_onei(t, t);
    }
    else if (arb_is_negative(acb_imagref(z)) ||
        (arb_is_zero(acb_imagref(z)) && arb_is_positive(acb_realref(z))))
    {
        acb_neg(t, nu);
        acb_exp_pi_i(t, t, prec);
        acb_mul_onei(t, t);
        acb_neg(t, t);
    }
    else
    {
        acb_exp_pi_i(t, nu, prec);
        acb_mul_onei(t, t);
        acb_neg(u, nu);
        acb_exp_pi_i(u, u, prec);
        acb_mul_onei(u, u);
        acb_neg(u, u);

        arb_union(acb_realref(t), acb_realref(t), acb_realref(u), prec);
        arb_union(acb_imagref(t), acb_imagref(t), acb_imagref(u), prec);
    }

    if (scaled)
    {
        acb_neg(u, z);
        acb_mul_2exp_si(u, u, 1);
        acb_exp(u, u, prec);
        acb_mul(A, t, u, prec);
        acb_one(B);
    }
    else
    {
        acb_exp_invexp(B, A, z, prec);
        acb_mul(A, A, t, prec);
    }

    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_bessel_i_asymp(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    acb_t A1, A2, C, U1, U2, s, t, u;
    int is_real, is_imag;

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

    if (!is_real && !scaled && arb_is_zero(acb_realref(z)) && acb_is_int(nu))
    {
        acb_mul_2exp_si(t, nu, -1);

        if (acb_is_int(t))
            is_real = 1;
        else
            is_imag = 1;
    }

    if (scaled)
        is_imag = 0;

    acb_hypgeom_bessel_i_asymp_prefactors(A1, A2, C, nu, z, scaled, prec);

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

        acb_mul_2exp_si(u, z, 1);
        acb_hypgeom_u_asymp(U1, s, t, u, -1, prec);
        acb_neg(u, u);
        acb_hypgeom_u_asymp(U2, s, t, u, -1, prec);

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
acb_hypgeom_bessel_i_0f1(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    acb_struct b[2];
    acb_t w, c, t;

    if (acb_is_int(nu) && arb_is_negative(acb_realref(nu)))
    {
        acb_init(t);
        acb_neg(t, nu);
        acb_hypgeom_bessel_i_0f1(res, t, z, scaled, prec);
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

    /* z^2/4 */
    acb_mul(w, z, z, prec);
    acb_mul_2exp_si(w, w, -2);

    acb_hypgeom_pfq_direct(t, NULL, 0, b, 2, w, -1, prec);

    if (scaled)
    {
        acb_neg(w, z);
        acb_exp(w, w, prec);
        acb_mul(t, t, w, prec);
    }

    acb_mul(res, t, c, prec);

    acb_clear(b + 0);
    acb_clear(b + 1);
    acb_clear(w);
    acb_clear(c);
    acb_clear(t);
}

void
acb_hypgeom_bessel_i_nointegration(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    mag_t zmag;

    mag_init(zmag);
    acb_get_mag(zmag, z);

    if (mag_cmp_2exp_si(zmag, 4) < 0 ||
        (mag_cmp_2exp_si(zmag, 64) < 0 && 2 * mag_get_d(zmag) < prec))
        acb_hypgeom_bessel_i_0f1(res, nu, z, scaled, prec);
    else
        acb_hypgeom_bessel_i_asymp(res, nu, z, scaled, prec);

    mag_clear(zmag);
}

void
_acb_hypgeom_bessel_i(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    acb_t res2;
    slong acc, max, t;

    acb_init(res2);

    acb_hypgeom_bessel_i_nointegration(res2, nu, z, scaled, prec);

    acc = acb_rel_accuracy_bits(res2);

    if (acc < 0.5 * prec)
    {
        max = prec;
        t = acb_rel_accuracy_bits(z);
        max = FLINT_MIN(max, t);
        t = acb_rel_accuracy_bits(nu);
        max = FLINT_MIN(max, t);

        if (max > 2 && acc < 0.5 * max)
        {
            if (acb_is_real(nu) && acb_is_real(z) && arf_cmp_2exp_si(arb_midref(acb_realref(nu)), -1) > 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 60) < 0)
            {

                arb_hypgeom_bessel_i_integration(acb_realref(res),
                        acb_realref(nu), acb_realref(z), scaled, prec);
                arb_zero(acb_imagref(res));

                if (acb_rel_accuracy_bits(res) > acb_rel_accuracy_bits(res2) ||
                    (acb_is_finite(res) && !acb_is_finite(res2)))
                {
                    acb_swap(res, res2);
                }
            }
        }
    }

    acb_swap(res, res2);
    acb_clear(res2);
}

void
acb_hypgeom_bessel_i(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    _acb_hypgeom_bessel_i(res, nu, z, 0, prec);
}

void
acb_hypgeom_bessel_i_scaled(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    _acb_hypgeom_bessel_i(res, nu, z, 1, prec);
}
