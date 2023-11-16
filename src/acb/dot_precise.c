/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_dot_precise(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep,
    acb_srcptr y, slong ystep, slong len, slong prec)
{
    arb_ptr tmp;
    slong i;

    tmp = flint_malloc(sizeof(arb_struct) * (4 * len));

    for (i = 0; i < len; i++)
    {
        tmp[0 * len + i] = *acb_realref(x + i * xstep);
        tmp[1 * len + i] = *acb_imagref(x + i * xstep);
        tmp[2 * len + i] = *acb_realref(y + i * ystep);
        arb_init(tmp + 3 * len + i);
        arb_neg(tmp + 3 * len + i, acb_imagref(y + i * ystep));
    }

    arb_dot_precise(acb_realref(res), initial == NULL ? NULL : acb_realref(initial), subtract,
        tmp, 1, tmp + 2 * len, 1, 2 * len, prec);

    for (i = 0; i < len; i++)
        arb_clear(tmp + 3 * len + i);

    for (i = 0; i < len; i++)
    {
        tmp[0 * len + i] = *acb_realref(x + i * xstep);
        tmp[1 * len + i] = *acb_imagref(x + i * xstep);
        tmp[2 * len + i] = *acb_imagref(y + i * ystep);
        tmp[3 * len + i] = *acb_realref(y + i * ystep);
    }

    arb_dot_precise(acb_imagref(res), initial == NULL ? NULL : acb_imagref(initial), subtract,
        tmp, 1, tmp + 2 * len, 1, 2 * len, prec);

    flint_free(tmp);
}
