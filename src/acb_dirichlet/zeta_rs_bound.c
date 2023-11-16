/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* Arias de Reyna, Theorem 4.2 */
void
acb_dirichlet_zeta_rs_bound(mag_t err, const acb_t s, slong K)
{
    arb_t a;
    mag_t c1, c2, c3;
    slong e;

    if (!arb_is_positive(acb_imagref(s)) || K < 1 || !acb_is_finite(s))
    {
        mag_inf(err);
        return;
    }

    arb_init(a);

    arb_add_ui(a, acb_realref(s), K, MAG_BITS);
    arb_sub_ui(a, a, 2, MAG_BITS);

    if (!arb_is_nonnegative(acb_realref(s)) && !arb_is_nonnegative(a))
    {
        mag_inf(err);
        arb_clear(a);
        return;
    }

    mag_init(c1);
    mag_init(c2);
    mag_init(c3);

    /* c1 = 1/7 2^(3 sigma / 2)  re(sigma) >= 0 */
    /* c1 < 1/2                  re(sigma) < 0  */
    arf_set_mag(arb_midref(a), arb_radref(acb_realref(s)));
    arf_add(arb_midref(a), arb_midref(a), arb_midref(acb_realref(s)), MAG_BITS, ARF_RND_CEIL);

    if (arf_sgn(arb_midref(a)) <= 0)
    {
        mag_set_ui_2exp_si(c1, 1, -1);
    }
    else if (arf_cmp_2exp_si(arb_midref(a), FLINT_BITS - 4) < 0)
    {
        mag_one(c1);
        mag_div_ui(c1, c1, 7);
        e = arf_get_si(arb_midref(a), ARF_RND_CEIL);
        mag_mul_2exp_si(c1, c1, (3 * e + 1) / 2);
        if (mag_cmp_2exp_si(c1, -1) < 0)
            mag_set_ui_2exp_si(c1, 1, -1);
    }
    else
    {
        mag_inf(c1);
    }

    /* c2 = 1 / ((10/11) sqrt(t/(2pi))) = (11/10) sqrt((2pi)/t) */
    arb_get_mag_lower(c3, acb_imagref(s));
    mag_const_pi(c2);
    mag_mul_2exp_si(c2, c2, 1);
    mag_div(c2, c2, c3);
    mag_sqrt(c2, c2);
    mag_mul_ui(c2, c2, 11);
    mag_div_ui(c2, c2, 10);

    /* c2 = c2^(K+1) */
    mag_pow_ui(c2, c2, K + 1);

    /* c3 = gamma((K+1)/2) */
    mag_fac_ui(c3, K / 2);

    mag_mul(err, c1, c2);
    mag_mul(err, err, c3);

    mag_clear(c1);
    mag_clear(c2);
    mag_clear(c3);
    arb_clear(a);
}

