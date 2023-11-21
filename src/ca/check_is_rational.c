/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

/* todo: fast check in number field */
truth_t
ca_check_is_rational(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        return T_TRUE;
    }
    else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        return nf_elem_is_rational(CA_NF_ELEM(x), CA_FIELD_NF(CA_FIELD(x, ctx))) ? T_TRUE : T_FALSE;
    }
    else
    {
        acb_t t;
        truth_t res;
        slong prec, prec_limit;

        res = T_UNKNOWN;

        acb_init(t);

        prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
        prec_limit = FLINT_MAX(prec_limit, 64);

        for (prec = 64; (prec <= prec_limit) && (res == T_UNKNOWN); prec *= 2)
        {
            ca_get_acb_raw(t, x, prec, ctx);

            if (!arb_contains_zero(acb_imagref(t)))
            {
                res = T_FALSE;
                break;
            }

            /* try qqbar computation */
            /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
            if (prec == 64)
            {
                qqbar_t a;
                qqbar_init(a);

                if (ca_get_qqbar(a, x, ctx))
                    res = qqbar_is_rational(a) ? T_TRUE : T_FALSE;

                qqbar_clear(a);
            }
        }

        acb_clear(t);

        return res;
    }
}

