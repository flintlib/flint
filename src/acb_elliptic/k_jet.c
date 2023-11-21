/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"

void
acb_elliptic_k_jet(acb_ptr w, const acb_t m, slong len, slong prec)
{
    acb_t t, u, msub1m, m2sub1;
    slong k, n;

    if (len < 1)
        return;

    if (len == 1)
    {
        acb_elliptic_k(w, m, prec);
        return;
    }

    if (acb_is_zero(m))
    {
        acb_const_pi(w, prec);
        acb_mul_2exp_si(w, w, -1);

        for (k = 1; k < len; k++)
        {
            acb_mul_ui(w + k, w + k - 1, (2 * k - 1) * (2 * k - 1), prec);
            acb_div_ui(w + k, w + k, 4 * k * k, prec);
        }

        return;
    }

    acb_init(t);
    acb_init(u);
    acb_init(msub1m);
    acb_init(m2sub1);

    acb_sub_ui(msub1m, m, 1, prec);
    acb_neg(t, msub1m);
    acb_sqrt(t, t, prec);
    acb_mul(msub1m, msub1m, m, prec);

    acb_mul_2exp_si(m2sub1, m, 1);
    acb_sub_ui(m2sub1, m2sub1, 1, prec);

    acb_agm1_cpx(w, t, 2, prec);

    /* pi M'(t) / (4 t M(t)^2) */
    acb_mul(u, w, w, prec);
    acb_mul(t, t, u, prec);
    acb_div(w + 1, w + 1, t, prec);

    acb_const_pi(u, prec);
    acb_mul(w + 1, w + 1, u, prec);
    acb_mul_2exp_si(w + 1, w + 1, -2);

    /* pi / (2 M(t)) */
    acb_const_pi(u, prec);
    acb_div(w, u, w, prec);
    acb_mul_2exp_si(w, w, -1);

    acb_inv(t, msub1m, prec);

    for (k = 2; k < len; k++)
    {
        n = k - 2;

        acb_mul_ui(w + k, w + n, (2 * n + 1) * (2 * n + 1), prec);

        acb_mul(u, w + n + 1, m2sub1, prec);
        acb_addmul_ui(w + k, u, (n + 1) * (n + 1) * 4, prec);

        acb_mul(w + k, w + k, t, prec);
        acb_div_ui(w + k, w + k, 4 * (n + 1) * (n + 2), prec);
        acb_neg(w + k, w + k);
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(msub1m);
    acb_clear(m2sub1);
}

