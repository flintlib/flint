/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_agm_sqrt_lowprec(acb_t r, const acb_t a, const acb_t root, slong prec)
{
    acb_t res;
    acb_t neg;

    acb_init(res);
    acb_init(neg);

    /* Take any square root, avoiding potentially massive precision loss
       if a intersects the default branch cut */

    if (arb_contains_zero(acb_imagref(a)) && arb_is_negative(acb_realref(a)))
    {
        acb_neg(res, a);
        acb_sqrt(res, res, prec);
        acb_mul_onei(res, res);
    }
    else
    {
        acb_sqrt(res, a, prec);
    }
    acb_neg(neg, res);

    if (acb_overlaps(root, res))
    {
        if (acb_overlaps(root, neg))
        {
            acb_indeterminate(r);
        }
        else
        {
            acb_set(r, res);
        }
    }
    else /* (!acb_overlaps(root, res)) */
    {
        if (!acb_ovelaps(root, neg))
        {
            acb_indeterminate(r);
        }
        else
        {
            acb_set(r, neg);
        }
    }
    
    acb_clear(res);
    acb_clear(neg);
}
