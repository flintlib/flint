/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"
#include "arb_hypgeom.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

static void
_acb_log_rising_correct_branch(acb_t res,
        const acb_t t_wrong, const acb_t z, ulong r, slong prec)
{
    acb_t f;
    arb_t pi, u, v;
    fmpz_t pi_mult;
    slong i, argprec;

    acb_init(f);

    arb_init(u);
    arb_init(pi);
    arb_init(v);

    fmpz_init(pi_mult);

    argprec = FLINT_MIN(prec, 40);

    arb_zero(u);
    for (i = 0; i < r; i++)
    {
        acb_add_ui(f, z, i, argprec);
        acb_arg(v, f, argprec);
        arb_add(u, u, v, argprec);
    }

    if (argprec == prec)
    {
        arb_set(acb_imagref(res), u);
    }
    else
    {
        arb_sub(v, u, acb_imagref(t_wrong), argprec);
        arb_const_pi(pi, argprec);
        arb_div(v, v, pi, argprec);

        if (arb_get_unique_fmpz(pi_mult, v))
        {
            arb_const_pi(v, prec);
            arb_mul_fmpz(v, v, pi_mult, prec);
            arb_add(acb_imagref(res), acb_imagref(t_wrong), v, prec);
        }
        else
        {
            arb_zero(u);
            for (i = 0; i < r; i++)
            {
                acb_add_ui(f, z, i, prec);
                acb_arg(v, f, prec);
                arb_add(u, u, v, prec);
            }
            arb_set(acb_imagref(res), u);
        }
    }

    arb_set(acb_realref(res), acb_realref(t_wrong));

    acb_clear(f);

    arb_clear(u);
    arb_clear(v);
    arb_clear(pi);

    fmpz_clear(pi_mult);
}

void
acb_hypgeom_log_rising_ui_jet_fallback(acb_ptr res, const acb_t z, slong r, slong len, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_set(t, z);

    if (len == 1)
    {
        acb_hypgeom_rising_ui_rec(res, t, r, prec);
        acb_log(res, res, prec);
    }
    else
    {
        acb_hypgeom_rising_ui_jet(res, t, r, len, prec);
        _acb_poly_log_series(res, res, FLINT_MIN(len, r + 1), len, prec);
    }

    _acb_log_rising_correct_branch(res, res, t, r, prec);

    acb_clear(t);
}

void
acb_hypgeom_log_rising_ui_jet(acb_ptr res, const acb_t z, ulong r, slong len, slong prec)
{
    double za, zb, sa, sb, ta, tb, ma, mb, zak;
    slong k, correction;
    int neg;

    if (r == 0 || len == 0)
    {
        _acb_vec_zero(res, len);
        return;
    }

    if (r == 1)
    {
        if (len == 1)
        {
            acb_log(res, z, prec);
        }
        else
        {
            acb_set(res, z);
            acb_one(res + 1);
            _acb_poly_log_series(res, res, 2, len, prec);
        }
        return;
    }

    if (arb_is_zero(acb_imagref(z)))
    {
        if (arb_is_positive(acb_realref(z)))
        {
            acb_hypgeom_rising_ui_jet(res, z, r, len, prec);
            _acb_poly_log_series(res, res, FLINT_MIN(len, r + 1), len, prec);
        }
        else if (arb_contains_int(acb_realref(z)))
        {
            _acb_vec_indeterminate(res, len);
        }
        else
        {
            arb_t t, u;

            arb_init(t);
            arb_init(u);

            arb_floor(u, acb_realref(z), prec);
            arb_neg(u, u);

            arb_set_ui(t, r);
            arb_min(u, u, t, prec);
            arb_const_pi(t, prec);
            arb_mul(t, u, t, prec);

            acb_hypgeom_rising_ui_jet(res, z, r, len, prec);
            _acb_vec_neg(res, res, FLINT_MIN(len, r + 1));
            _acb_poly_log_series(res, res, FLINT_MIN(len, r + 1), len, prec);

            arb_swap(acb_imagref(res), t);

            arb_clear(t);
            arb_clear(u);
        }

        return;
    }

    /* We use doubles if it is safe.
       - No overflow/underflow possible.
       - Input is accurate enough (and not too close to the real line).
         Note: the relative error for a complex floating-point
         multiplication is bounded by sqrt(5) * eps, and we basically only
         need to determine the result to within one quadrant. */

    /* todo: wide */
    if (prec <= 20 || acb_rel_accuracy_bits(z) < 30 || arb_rel_accuracy_bits(acb_imagref(z)) < 30)
    {
        acb_hypgeom_log_rising_ui_jet_fallback(res, z, r, len, prec);
        return;
    }

    za = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    zb = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);

    if (!(r <= 1e6 && za <= 1e6 && za >= -1e6 && zb <= 1e6 && zb >= -1e6 && (zb > 1e-6 || zb < -1e-6)))
    {
        acb_hypgeom_log_rising_ui_jet_fallback(res, z, r, len, prec);
        return;
    }

    sa = za;
    sb = zb;
    correction = 0;
    neg = 0;

    for (k = 1; k < r; k++)
    {
        zak = za + k;

        ta = sa * zak - sb * zb;
        tb = sb * zak + sa * zb;

        if (zb > 0.0)
        {
            if (sb >= 0.0 && tb < 0.0)
                correction += 2;
        }
        else
        {
            if (sb < 0.0 && tb >= 0.0)
                correction += 2;
        }

        sa = ta;
        sb = tb;

        if (k % 4 == 0)
        {
            ma = fabs(sa);
            mb = fabs(sb);

            /* Rescale to protect against overflow. */
            if (ma > mb)
                ma = 1.0 / ma;
            else
                ma = 1.0 / mb;

            sa *= ma;
            sb *= ma;
        }
    }

    if (sa < 0.0)
    {
        neg = 1;

        if ((zb > 0.0 && sb >= 0.0) || (zb < 0.0 && sb < 0.0))
            correction += 1;
        else
            correction -= 1;
    }

    if (len == 1)
    {
        acb_hypgeom_rising_ui_rec(res, z, r, prec);
        if (neg)
            acb_neg(res, res);
        acb_log(res, res, prec);
    }
    else
    {
        acb_hypgeom_rising_ui_jet(res, z, r, len, prec);
        if (neg)
            _acb_vec_neg(res, res, FLINT_MIN(len, r + 1));
        _acb_poly_log_series(res, res, FLINT_MIN(len, r + 1), len, prec);
    }

    if (zb < 0.0)
        correction = -correction;

    if (correction != 0)
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, prec);
        arb_addmul_si(acb_imagref(res), t, correction, prec);
        arb_clear(t);
    }
}

void
acb_hypgeom_log_rising_ui(acb_ptr res, const acb_t z, ulong r, slong prec)
{
    acb_hypgeom_log_rising_ui_jet(res, z, r, 1, prec);
}

