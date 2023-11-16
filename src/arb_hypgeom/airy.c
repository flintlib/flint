/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
arb_hypgeom_airy(arb_t ai, arb_t aip, arb_t bi, arb_t bip, const arb_t z, slong prec)
{
    acb_struct tmp[5];

    if (ai != NULL) acb_init(tmp);
    if (aip != NULL) acb_init(tmp + 1);
    if (bi != NULL) acb_init(tmp + 2);
    if (bip != NULL) acb_init(tmp + 3);

    acb_init(tmp + 4);
    acb_set_arb(tmp + 4, z);

    acb_hypgeom_airy(ai ? tmp : NULL,
                     aip ? tmp + 1 : NULL,
                     bi ? tmp + 2 : NULL,
                     bip ? tmp + 3 : NULL,
                     tmp + 4, prec);

    if (ai != NULL) arb_set(ai, acb_realref(tmp));
    if (aip != NULL) arb_set(aip, acb_realref(tmp + 1));
    if (bi != NULL) arb_set(bi, acb_realref(tmp + 2));
    if (bip != NULL) arb_set(bip, acb_realref(tmp + 3));

    if (ai != NULL) acb_clear(tmp);
    if (aip != NULL) acb_clear(tmp + 1);
    if (bi != NULL) acb_clear(tmp + 2);
    if (bip != NULL) acb_clear(tmp + 3);

    acb_clear(tmp + 4);
}

