/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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

    if (CA_IS_QQ(x, ctx) && CA_IS_QQ(y, ctx))
    {
        return fmpq_equal(CA_FMPQ(x), CA_FMPQ(y)) ? T_TRUE : T_FALSE;
    }

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        if (CA_IS_UNKNOWN(x) || CA_IS_UNKNOWN(y))
            return T_UNKNOWN;

        if (CA_IS_SIGNED_INF(x) && CA_IS_SIGNED_INF(y))
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
    if (x->field == y->field && CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
        return T_FALSE;

    /* Rational number field elements *should* have been demoted to QQ
       automatically, but let's do a comparison as a precaution. */
    if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)) && CA_IS_QQ(y, ctx))
        return nf_elem_equal_fmpq(CA_NF_ELEM(x), CA_FMPQ(y), CA_FIELD_NF(CA_FIELD(x, ctx))) ? T_TRUE : T_FALSE;

    if (CA_FIELD_IS_NF(CA_FIELD(y, ctx)) && CA_IS_QQ(x, ctx))
        return nf_elem_equal_fmpq(CA_NF_ELEM(y), CA_FMPQ(x), CA_FIELD_NF(CA_FIELD(y, ctx))) ? T_TRUE : T_FALSE;

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
    if (0 && x_alg == T_TRUE && y_alg == T_TRUE)
    {
        /* ...
        qqbar_t a, b;
        qqbar_init(a);
        qqbar_init(b);

        if (ca_get_qqbar(a, x, ctx))
        {
            if (ca_get_qqbar(b, y, ctx))
            {
                int eq = qqbar_equal(a, b);
                qqbar_clear(a);
                qqbar_clear(b);
                return eq ? T_TRUE : T_FALSE;
            }
        }

        qqbar_clear(a);
        qqbar_clear(b);
        */
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
