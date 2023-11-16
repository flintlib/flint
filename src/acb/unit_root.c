/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "acb.h"

static void
_acb_unit_root(acb_t res, ulong order, slong prec)
{
    fmpq_t t;
    fmpq_init(t);
    fmpq_set_si(t, 2, order);
    arb_sin_cos_pi_fmpq(acb_imagref(res), acb_realref(res), t, prec);
    fmpq_clear(t);
}

void
acb_unit_root(acb_t res, ulong order, slong prec)
{
    switch (order)
    {
       case 1:
           acb_one(res);
           break;
       case 2:
           acb_set_si(res, -1);
           break;
       case 4:
           acb_onei(res);
           break;
       default:
           _acb_unit_root(res, order, prec);
           break;
    }
}
