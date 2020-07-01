/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "ca.h"

void
ca_get_acb(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)
{
    slong wp, initial, maxprec;

    acb_indeterminate(res);

    initial = prec * 1.05 + 30;
    maxprec = FLINT_MAX(2 * initial, ctx->options[CA_OPT_PREC_LIMIT]);

    for (wp = initial; wp <= maxprec; wp *= 2)
    {
        ca_get_acb_raw(res, x, wp, ctx);

        if (acb_rel_accuracy_bits(res) >= prec)
            break;
    }
}

void
ca_get_acb_accurate_parts(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx)
{
    slong wp, initial, maxprec;

    acb_indeterminate(res);

    initial = prec * 1.05 + 30;
    maxprec = FLINT_MAX(2 * initial, ctx->options[CA_OPT_PREC_LIMIT]);

    for (wp = initial; wp <= maxprec; wp *= 2)
    {
        ca_get_acb_raw(res, x, wp, ctx);

        if (arb_rel_accuracy_bits(acb_realref(res)) >= prec &&
            arb_rel_accuracy_bits(acb_imagref(res)) >= prec)
            break;
    }
}
