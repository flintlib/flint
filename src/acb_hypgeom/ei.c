/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_ei_asymp(acb_t res, const acb_t z, slong prec)
{
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    acb_one(t);
    acb_neg(u, z);

    acb_hypgeom_u_asymp(u, t, t, u, -1, prec);
    acb_div(u, u, z, prec);
    acb_exp(t, z, prec);
    acb_mul(u, u, t, prec);

    if (arb_is_zero(acb_imagref(z)))
    {
        arb_zero(acb_imagref(u));
    }
    else if (arb_is_positive(acb_imagref(z)))
    {
        acb_const_pi(t, prec);
        arb_add(acb_imagref(u), acb_imagref(u), acb_realref(t), prec);
    }
    else if (arb_is_negative(acb_imagref(z)))
    {
        acb_const_pi(t, prec);
        arb_sub(acb_imagref(u), acb_imagref(u), acb_realref(t), prec);
    }
    else
    {
        /* add [-pi,pi] i */
        acb_const_pi(t, prec);
        arb_add_error(acb_imagref(u), acb_realref(t));
    }

    acb_swap(res, u);

    acb_clear(t);
    acb_clear(u);
}

/*
Ei(z) = z 2F2(1,1,2,2,z) + 0.5[log(z)-log(1/z)] + gamma
*/
void
acb_hypgeom_ei_2f2(acb_t res, const acb_t z, slong prec)
{
    acb_t a, t;
    acb_struct b[2];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(t);

    acb_one(a);
    acb_set_ui(b, 2);
    acb_set_ui(b + 1, 2);

    acb_hypgeom_pfq_direct(t, a, 1, b, 2, z, -1, prec);
    acb_mul(t, t, z, prec);

    arb_const_euler(acb_realref(a), prec);
    arb_add(acb_realref(t), acb_realref(t), acb_realref(a), prec);

    if (arb_is_zero(acb_imagref(z)))
    {
        if (arb_is_positive(acb_realref(z)))
        {
            acb_log(a, z, prec);
        }
        else
        {
            /* ok if overlapping zero -- will be indeterminate either way */
            acb_neg(a, z);
            acb_log(a, a, prec);
            arb_zero(acb_imagref(a));
        }

        acb_add(t, t, a, prec);
    }
    else if (arb_is_nonzero(acb_imagref(z))
        || arb_is_nonnegative(acb_realref(z))) /* not overlapping (-inf,0] */
    {
        acb_log(a, z, prec);
        acb_add(t, t, a, prec);
    }
    else
    {
        acb_log(a, z, prec);
        arb_zero(acb_imagref(a));

        acb_const_pi(b, prec);
        arb_add_error(acb_imagref(a), acb_realref(b));

        acb_add(t, t, a, prec);
    }

    acb_swap(res, t);

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(t);
}

void
acb_hypgeom_ei(acb_t res, const acb_t z, slong prec)
{
    if (acb_hypgeom_u_use_asymp(z, prec))
        acb_hypgeom_ei_asymp(res, z, prec);
    else
        acb_hypgeom_ei_2f2(res, z, prec);
}

