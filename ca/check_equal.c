/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

truth_t
ca_check_equal(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    acb_t u, v;
    ca_t t;
    truth_t res;
    truth_t x_alg, y_alg;
    slong prec;

    if (x->field == CA_FIELD_ID_QQ && y->field == CA_FIELD_ID_QQ)
    {
        return fmpq_equal(CA_FMPQ(x), CA_FMPQ(y)) ? T_TRUE : T_FALSE;
    }

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        if ((x->field & CA_UNKNOWN) || (y->field & CA_UNKNOWN))
            return T_UNKNOWN;

        if ((x->field & CA_SIGNED_INF) && (y->field & CA_SIGNED_INF))
        {
            ca_t xsign, ysign;

            *xsign = *x;
            *ysign = *y;

            xsign->field &= ~CA_SPECIAL;
            ysign->field &= ~CA_SPECIAL;

            return ca_check_equal(xsign, ysign, ctx);
        }

        if (x->field == y->field)
            return T_TRUE;
        else
            return T_FALSE;
    }

    if (ca_equal_repr(x, y, ctx))
        return T_TRUE;

    /* same algebraic number field ==> sufficient to compare representation */
    if (x->field == y->field && ctx->fields[x->field].type == CA_FIELD_TYPE_NF)
        return T_FALSE;

    res = T_UNKNOWN;

    acb_init(u);
    acb_init(v);

    /* for (prec = 64; (prec <= ctx->options[CA_OPT_PREC_LIMIT]) && (res == T_UNKNOWN); prec *= 2) */
    prec = 64;

    {
        ca_get_acb_raw(u, x, prec, ctx);
        ca_get_acb_raw(v, y, prec, ctx);

        if (!acb_overlaps(u, v))
        {
            res = T_FALSE;
        }
    }

    acb_clear(u);
    acb_clear(v);

    x_alg = ca_check_is_algebraic(x, ctx);
    y_alg = ca_check_is_algebraic(y, ctx);

    if ((x_alg == T_TRUE && y_alg == T_FALSE) ||
        (x_alg == T_FALSE && y_alg == T_TRUE))
        return T_FALSE;

    /* todo: try qqbar computation */
    /* we may want to do this selectively; in some cases cancellation in
       computing x-y will be helpful; in other cases, subtracting the
       terms will make life more difficult */
    if (x_alg == T_TRUE && y_alg == T_TRUE)
    {
        /* ... */
    }

    if (res == T_UNKNOWN)
    {
        /* check_is_zero may have additional heuristics */
        ca_init(t, ctx);
        ca_sub(t, x, y, ctx);
        res = ca_check_is_zero(t, ctx);
        ca_clear(t, ctx);
    }

    return res;
}

