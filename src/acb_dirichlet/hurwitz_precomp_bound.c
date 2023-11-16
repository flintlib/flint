/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_hurwitz_precomp_bound(mag_t res, const acb_t s,
    slong A, slong K, slong N)
{
    acb_t s1;
    mag_t x, t, TK, C;
    slong sigmaK;
    arf_t u;

    if (A < 1 || K < 1 || N < 1)
    {
        mag_inf(res);
        return;
    }

    /* sigmaK = re(s) + K, floor bound */
    arf_init(u);
    arf_set_mag(u, arb_radref(acb_realref(s)));
    arf_sub(u, arb_midref(acb_realref(s)), u, MAG_BITS, ARF_RND_FLOOR);
    arf_add_ui(u, u, K, MAG_BITS, ARF_RND_FLOOR);
    if (arf_cmp_ui(u, 2) < 0 || arf_cmp_2exp_si(u, FLINT_BITS - 2) > 0)
    {
        mag_inf(res);
        arf_clear(u);
        return;
    }
    sigmaK = arf_get_si(u, ARF_RND_FLOOR);
    arf_clear(u);

    acb_init(s1);
    mag_init(x);
    mag_init(t);
    mag_init(TK);
    mag_init(C);

    /* With N grid points, we will have |x| <= 1 / (2N). */
    mag_one(x);
    mag_div_ui(x, x, 2 * N);

    /* T(K) = |x|^K |(s)_K| / K! * [1/A^(sigma+K) + ...] */
    mag_pow_ui(TK, x, K);
    acb_rising_ui_get_mag(t, s, K);
    mag_mul(TK, TK, t);
    mag_rfac_ui(t, K);
    mag_mul(TK, TK, t);
    /* Note: here we assume that mag_hurwitz_zeta_uiui uses an error bound
       that is at least as large as the one used in the proof. */
    mag_hurwitz_zeta_uiui(t, sigmaK, A);
    mag_mul(TK, TK, t);

    /* C = |x|/A (1 + 1/(K+sigma+A-1)) (1 + |s-1|/(K+1)) */

    mag_div_ui(C, x, A);

    mag_one(t);
    mag_div_ui(t, t, sigmaK + A - 1);
    mag_add_ui(t, t, 1);
    mag_mul(C, C, t);

    acb_sub_ui(s1, s, 1, MAG_BITS);
    acb_get_mag(t, s1);
    mag_div_ui(t, t, K + 1);
    mag_add_ui(t, t, 1);
    mag_mul(C, C, t);

    mag_geom_series(t, C, 0);
    mag_mul(res, TK, t);

    acb_clear(s1);
    mag_clear(x);
    mag_clear(t);
    mag_clear(TK);
    mag_clear(C);
}

