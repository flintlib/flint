/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_inv(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    truth_t is_zero;
    slong field_index;
    ulong xfield;
    ca_field_srcptr res_field;

    xfield = x->field;

    if (xfield == CA_FIELD_ID_QQ)
    {
        if (fmpq_is_zero(CA_FMPQ(x)))
        {
            ca_uinf(res, ctx);
        }
        else
        {
            _ca_make_fmpq(res, ctx);
            fmpq_inv(CA_FMPQ(res), CA_FMPQ(x));
        }
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        if ((xfield & CA_SIGNED_INF) || (xfield & CA_UNSIGNED_INF))
            ca_zero(res, ctx);
        else
            ca_set(res, x, ctx);
        return;
    }

    is_zero = ca_check_is_zero(x, ctx);

    if (is_zero == T_TRUE)
    {
        ca_uinf(res, ctx);
        return;
    }
    else if (is_zero == T_UNKNOWN)
    {
        ca_unknown(res, ctx);
        return;
    }

    field_index = xfield;

    _ca_make_field_element(res, field_index, ctx);
    res_field = CA_FIELD(res, ctx);

    if (CA_FIELD_IS_QQ(res_field))
    {
        fmpq_inv(CA_FMPQ(res), CA_FMPQ(x));
    }
    else if (CA_FIELD_IS_NF(res_field))
    {
        nf_elem_inv(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_FIELD_NF(res_field));
    }
    else
    {
        fmpz_mpoly_q_inv(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_FIELD_MCTX(res_field, ctx));
    }
}

