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
_ca_check_is_zero_qqbar(const ca_t x, ca_ctx_t ctx)
{
    qqbar_t t;
    truth_t res;
    qqbar_init(t);

    if (ca_get_qqbar(t, x, ctx))
        res = qqbar_is_zero(t) ? T_TRUE : T_FALSE;
    else
        res = T_UNKNOWN;

    qqbar_clear(t);
    return res;
}

truth_t
ca_check_is_zero(const ca_t x, ca_ctx_t ctx)
{
    acb_t v;
    truth_t res;
    slong prec, prec_limit;

    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }

    if (CA_IS_QQ(x, ctx))
    {
        if (fmpq_is_zero(CA_FMPQ(x)))
            return T_TRUE;
        else
            return T_FALSE;
    }

    if (CA_IS_QQ_I(x, ctx))
    {
        const fmpz *n;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));

        if (fmpz_is_zero(n) && fmpz_is_zero(n + 1))
            return T_TRUE;

        return T_FALSE;
    }

    if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        if (nf_elem_is_zero(CA_NF_ELEM(x), CA_FIELD_NF(CA_FIELD(x, ctx))))
            return T_TRUE;
        else
            return T_FALSE;
    }

    res = T_UNKNOWN;

    /* todo: in the following, we should extract the numerator; the denominator is irrelevant */

    acb_init(v);

    prec_limit = ctx->options[CA_OPT_PREC_LIMIT];
    prec_limit = FLINT_MAX(prec_limit, 64);

    for (prec = 64; (prec <= prec_limit) && (res == T_UNKNOWN); prec *= 2)
    {
        ca_get_acb_raw(v, x, prec, ctx);

        if (!acb_contains_zero(v))
        {
            res = T_FALSE;
            break;
        }

        /* try qqbar computation */
        /* todo: precision to do this should depend on complexity of the polynomials, degree of the elements... */
        if (prec == 64)
        {
            res = _ca_check_is_zero_qqbar(x, ctx);
        }
    }

    acb_clear(v);

    return res;
}

