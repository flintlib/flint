/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

static void
acb_log_sin_pi_half(acb_t res, const acb_t z, slong prec, int upper)
{
    acb_t t, u, zmid;
    arf_t n;
    arb_t pi;

    acb_init(t);
    acb_init(u);
    acb_init(zmid);
    arf_init(n);
    arb_init(pi);

    arf_set(arb_midref(acb_realref(zmid)), arb_midref(acb_realref(z)));
    arf_set(arb_midref(acb_imagref(zmid)), arb_midref(acb_imagref(z)));

    arf_floor(n, arb_midref(acb_realref(zmid)));
    arb_sub_arf(acb_realref(zmid), acb_realref(zmid), n, prec);

    arb_const_pi(pi, prec);

    if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(zmid)), 2) < 1)
    {
        acb_sin_pi(t, zmid, prec);
        acb_log(t, t, prec);
    }
    else  /* i*pi*(z-0.5) + log((1-exp(-2i*pi*z))/2) */
    {
        acb_mul_2exp_si(t, zmid, 1);
        acb_neg(t, t);

        if (upper)
            acb_conj(t, t);

        acb_exp_pi_i(t, t, prec);
        acb_sub_ui(t, t, 1, prec);
        acb_neg(t, t);

        acb_mul_2exp_si(t, t, -1);

        acb_log(t, t, prec);
        acb_one(u);
        acb_mul_2exp_si(u, u, -1);
        acb_sub(u, zmid, u, prec);
        if (upper)
            acb_conj(u, u);
        acb_mul_onei(u, u);
        acb_addmul_arb(t, u, pi, prec);
        if (upper)
            acb_conj(t, t);
    }

    if (upper)
        arb_submul_arf(acb_imagref(t), pi, n, prec);
    else
        arb_addmul_arf(acb_imagref(t), pi, n, prec);

    /* propagated error bound from the derivative pi cot(pi z) */
    if (!acb_is_exact(z))
    {
        mag_t zm, um;

        mag_init(zm);
        mag_init(um);

        acb_cot_pi(u, z, prec);
        acb_mul_arb(u, u, pi, prec);

        mag_hypot(zm, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
        acb_get_mag(um, u);
        mag_mul(um, um, zm);

        acb_add_error_mag(t, um);

        mag_clear(zm);
        mag_clear(um);
    }

    acb_set(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(zmid);
    arf_clear(n);
    arb_clear(pi);
}

void
acb_log_sin_pi(acb_t res, const acb_t z, slong prec)
{
    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (arb_is_positive(acb_imagref(z)) ||
        (arb_is_zero(acb_imagref(z)) && arb_is_negative(acb_realref(z))))
    {
        acb_log_sin_pi_half(res, z, prec, 1);
    }
    else if (arb_is_negative(acb_imagref(z)) ||
        (arb_is_zero(acb_imagref(z)) && arb_is_positive(acb_realref(z))))
    {
        acb_log_sin_pi_half(res, z, prec, 0);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_log_sin_pi_half(t, z, prec, 1);
        acb_log_sin_pi_half(res, z, prec, 0);
        arb_union(acb_realref(res), acb_realref(res), acb_realref(t), prec);
        arb_union(acb_imagref(res), acb_imagref(res), acb_imagref(t), prec);
        acb_clear(t);
    }
}

