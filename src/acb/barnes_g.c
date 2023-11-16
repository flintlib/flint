/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_poly.h"

void
_arb_const_zeta_minus_one_eval(arb_t y, slong prec)
{
    acb_struct z[2];
    acb_t s, a;
    acb_init(z + 0);
    acb_init(z + 1);
    acb_init(s);
    acb_init(a);
    acb_set_si(s, -1);
    acb_one(a);
    _acb_poly_zeta_cpx_series(z, s, a, 0, 2, prec + 20);
    arb_set(y, acb_realref(z + 1));
    acb_clear(z + 0);
    acb_clear(z + 1);
    acb_clear(s);
    acb_clear(a);
}

ARB_DEF_CACHED_CONSTANT(_arb_const_zeta_minus_one, _arb_const_zeta_minus_one_eval)

/* LogG(z) = (z-1) LogGamma(z) - zeta'(-1,z) + zeta'(-1)
   LogG'(z) = (1/2)(-2z + 1 + log(2pi)) + (z-1) digamma(z) */

void
_acb_log_barnes_g_zeta(acb_t res, const acb_t z, slong prec)
{
    acb_struct t[3];
    acb_init(t + 0);
    acb_init(t + 1);
    acb_init(t + 2);

    acb_set_si(t + 2, -1);
    _acb_poly_zeta_cpx_series(t, t + 2, z, 0, 2, prec);

    _arb_const_zeta_minus_one(acb_realref(t), prec);
    arb_zero(acb_imagref(t));
    acb_sub(t, t, t + 1, prec);

    acb_lgamma(t + 1, z, prec);
    acb_sub_ui(t + 2, z, 1, prec);
    acb_addmul(t, t + 1, t + 2, prec);

    acb_set(res, t);

    acb_clear(t + 0);
    acb_clear(t + 1);
    acb_clear(t + 2);
}

void
_acb_barnes_g_ui_rec(acb_t res, ulong n, slong prec)
{
    acb_t t;
    ulong k;

    acb_init(t);

    acb_one(res);
    acb_one(t);

    for (k = 3; k < n; k++)
    {
        acb_mul_ui(t, t, k - 1, prec);
        acb_mul(res, res, t, prec);
    }

    acb_clear(t);
}

void
acb_log_barnes_g(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_int(z))
    {
        if (arb_is_nonpositive(acb_realref(z)))
        {
            acb_indeterminate(res);
            return;
        }

        if (arf_cmpabs_ui(arb_midref(acb_realref(z)), prec) < 0)
        {
            _acb_barnes_g_ui_rec(res,
                arf_get_si(arb_midref(acb_realref(z)), ARF_RND_DOWN), prec);
            acb_log(res, res, prec);
            return;
        }
    }

    _acb_log_barnes_g_zeta(res, z, prec);
}

void
acb_barnes_g(acb_t res, const acb_t z, slong prec)
{
    int real;

    real = acb_is_real(z);

    if (acb_is_int(z))
    {
        if (arb_is_nonpositive(acb_realref(z)))
        {
            acb_zero(res);
            return;
        }

        if (arf_cmpabs_ui(arb_midref(acb_realref(z)), prec) < 0)
        {
            _acb_barnes_g_ui_rec(res,
                arf_get_si(arb_midref(acb_realref(z)), ARF_RND_DOWN), prec);
            return;
        }
    }

    _acb_log_barnes_g_zeta(res, z, prec);
    acb_exp(res, res, prec);

    if (real)
        arb_zero(acb_imagref(res));
}

