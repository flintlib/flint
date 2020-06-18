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
ca_check_is_zero(const ca_t x, ca_ctx_t ctx)
{
    acb_t v;
    truth_t res;
    slong prec;

    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }

    if (x->field == CA_FIELD_ID_QQ)
    {
        if (fmpq_is_zero(CA_FMPQ(x)))
            return T_TRUE;
        else
            return T_FALSE;
    }

    if (x->field == CA_FIELD_ID_QQ_I)
    {
        const fmpz *n;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));

        if (fmpz_is_zero(n) && fmpz_is_zero(n + 1))
            return T_TRUE;

        return T_FALSE;
    }

    if ((ctx->fields + x->field)->type == CA_FIELD_TYPE_NF)
    {
        if (nf_elem_is_zero(CA_NF_ELEM(x), CA_FIELD_NF(ctx->fields + x->field)))
            return T_TRUE;
        else
            return T_FALSE;
    }

    res = T_UNKNOWN;

    acb_init(v);

    for (prec = 64; (prec <= ctx->options[CA_OPT_PREC_LIMIT]) && (res == T_UNKNOWN); prec *= 2)
    {
        ca_get_acb_raw(v, x, prec, ctx);

        if (!acb_contains_zero(v))
        {
            res = T_FALSE;
        }
    }

    /* todo: try exact simplifications (e.g. in qqbar) */

    acb_clear(v);

    return res;
}

