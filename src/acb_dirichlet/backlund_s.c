/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* S(t) = N(t) - theta(t)/pi - 1 */
static void
_backlund_s(arb_t res, const arb_t t, slong prec)
{
    arb_t N;
    acb_t z;

    arb_init(N);
    acb_init(z);

    acb_dirichlet_zeta_nzeros(N, t, prec);
    acb_set_arb(z, t);
    acb_dirichlet_hardy_theta(z, z, NULL, NULL, 1, prec);
    arb_const_pi(acb_imagref(z), prec);
    arb_div(acb_realref(z), acb_realref(z), acb_imagref(z), prec);
    arb_sub(res, N, acb_realref(z), prec);
    arb_sub_ui(res, res, 1, prec);

    arb_clear(N);
    acb_clear(z);
}

void
acb_dirichlet_backlund_s(arb_t res, const arb_t t, slong prec)
{
    mag_t m, b;

    mag_init(m);
    mag_init(b);

    arb_get_mag(m, t);

    if (!arb_is_nonnegative(t))
    {
        arb_indeterminate(res);
    }
    else if (mag_cmp_2exp_si(m, 6) < 0) /* If t is small then N(t) is fast. */
    {
        _backlund_s(res, t, prec);
    }
    else
    {
        /*
         * If the error radius of t is greater than 8/log(t) then use a fast
         * bound. The threshold doesn't matter for correctness, but note that
         * the average spacing of zeros is 2*pi/log(t) asymptotically.
         * Otherwise if the radius is smaller use S(t) = N(t) - theta(t)/pi - 1.
         */
        mag_log(b, m);
        mag_mul_2exp_si(b, b, -3);
        mag_inv(b, b);
        if (mag_cmp(arb_radref(t), b) > 0)
        {
            arb_zero(res);
            acb_dirichlet_backlund_s_bound(arb_radref(res), t);
        }
        else
        {
            acb_t z;
            acb_init(z);
            /*
             * The difference (N(t) - theta(t)/pi) loses precision from
             * cancellation. Add extra bits to mitigate the loss.
             */
            acb_set_arb(z, t);
            acb_dirichlet_hardy_theta(z, z, NULL, NULL, 1, prec);
            _backlund_s(res, t, prec + mag_get_d_log2_approx(m));
            acb_clear(z);
        }
    }

    mag_clear(m);
    mag_clear(b);
}
