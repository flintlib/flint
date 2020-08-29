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
ca_check_le(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    acb_t v, w;
    truth_t res;
    slong prec, prec_limit;

    if (CA_IS_QQ(x, ctx) && CA_IS_QQ(y, ctx))
    {
        return fmpq_cmp(CA_FMPQ(x), CA_FMPQ(y)) <= 0 ? T_TRUE : T_FALSE;
    }

    if (CA_IS_SPECIAL(x) || CA_IS_SPECIAL(y))
    {
        /* Not implemented */
        return T_UNKNOWN;
    }

    res = T_UNKNOWN;

    acb_init(v);
    acb_init(w);

    prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
    prec_limit = FLINT_MAX(prec_limit, 64);

    for (prec = 64; (prec <= prec_limit) && (res == T_UNKNOWN); prec *= 2)
    {
        ca_get_acb_raw(v, x, prec, ctx);
        ca_get_acb_raw(w, y, prec, ctx);

        if (!arb_contains_zero(acb_imagref(v)) ||
            !arb_contains_zero(acb_imagref(w)))
        {
            res = T_FALSE;
            break;
        }

        if (arb_is_zero(acb_imagref(v)) && arb_is_zero(acb_imagref(w)))
        {
            if (arb_le(acb_realref(v), acb_realref(w)))
            {
                res = T_TRUE;
                break;
            }
            else if (arb_gt(acb_realref(v), acb_realref(w)))
            {
                res = T_FALSE;
                break;
            }
        }

        /* try qqbar computation */
        /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
        if (prec == 64)
        {
            if (ca_can_evaluate_qqbar(x, ctx) && ca_can_evaluate_qqbar(y, ctx))
            {
                qqbar_t t, u;

                qqbar_init(t);
                qqbar_init(u);

                if (ca_get_qqbar(t, x, ctx))
                {
                    if (!qqbar_is_real(t))
                    {
                        res = T_FALSE;
                    }
                    else if (ca_get_qqbar(u, y, ctx))
                    {
                        if (!qqbar_is_real(u))
                        {
                            res = T_FALSE;
                        }
                        else
                        {
                            res = qqbar_cmp_re(t, u) <= 0 ? T_TRUE : T_FALSE;
                        }
                    }
                }

                qqbar_clear(t);
                qqbar_clear(u);
            }
        }
    }

    acb_clear(v);
    acb_clear(w);

    return res;
}
