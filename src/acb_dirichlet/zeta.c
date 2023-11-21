/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_dirichlet.h"

#ifdef __GNUC__
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

void acb_zeta_si(acb_t z, slong s, slong prec);

static void
_acb_dirichlet_zeta(acb_t res, const acb_t s, slong prec)
{
    acb_t a;
    acb_init(a);
    acb_one(a);
    _acb_poly_zeta_cpx_series(res, s, a, 0, 1, prec);
    acb_clear(a);
}

/* todo: use euler product for complex s, and check efficiency
   for large negative integers */
void
acb_dirichlet_zeta(acb_t res, const acb_t s, slong prec)
{
    double cutoff;

    if (acb_is_int(s) &&
        arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), FLINT_BITS - 1) < 0)
    {
        acb_zeta_si(res, arf_get_si(arb_midref(acb_realref(s)), ARF_RND_DOWN), prec);
        return;
    }

    if (arb_contains_zero(acb_imagref(s)) && arb_contains_si(acb_realref(s), 1))
    {
        acb_indeterminate(res);
        return;
    }

    cutoff = 24.0 * prec * sqrt(prec);

    if (arf_cmpabs_d(arb_midref(acb_imagref(s)), cutoff) >= 0 &&
        arf_cmpabs_d(arb_midref(acb_realref(s)), 10 + prec * 0.1) <= 0)
    {
        acb_dirichlet_zeta_rs(res, s, 0, prec);
        return;
    }

    if ((arf_sgn(arb_midref(acb_realref(s))) < 0) &&
        !acb_contains_zero(s))
    {
        acb_t t, u, v;
        slong wp = prec + 6;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_sub_ui(t, s, 1, wp);

        /* 2 * (2pi)^(s-1) */
        arb_const_pi(acb_realref(u), wp);
        acb_mul_2exp_si(u, u, 1);
        acb_pow(u, u, t, wp);
        acb_mul_2exp_si(u, u, 1);

        /* sin(pi*s/2) */
        acb_mul_2exp_si(v, s, -1);
        acb_sin_pi(v, v, wp);
        acb_mul(u, u, v, wp);

        /* gamma(1-s) zeta(1-s) */
        acb_neg(t, t);
        acb_gamma(v, t, wp);
        acb_mul(u, u, v, wp);
        _acb_dirichlet_zeta(v, t, prec);
        acb_mul(res, u, v, prec);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
    }
    else
    {
        _acb_dirichlet_zeta(res, s, prec);
    }
}

