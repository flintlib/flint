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
ca_neg(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    slong field_index;
    ca_field_type_t type;
    ulong xfield;

    xfield = x->field;

    if (xfield == CA_FIELD_ID_QQ)
    {
        _ca_make_fmpq(res, ctx);
        fmpq_neg(CA_FMPQ(res), CA_FMPQ(x));
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        if ((xfield & CA_UNKNOWN) || (xfield & CA_UNDEFINED) || (xfield & CA_UNSIGNED_INF))
        {
            ca_set(res, x, ctx);
            return;
        }
    }

    field_index = xfield & ~CA_SPECIAL;
    type = ctx->fields[field_index].type;

    _ca_make_field_element(res, field_index, ctx);
    res->field = xfield;  /* set special flags */

    if (type == CA_FIELD_TYPE_QQ)
    {
        fmpq_neg(CA_FMPQ(res), CA_FMPQ(x));
    }
    else if (type == CA_FIELD_TYPE_NF)
    {
        nf_elem_neg(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_FIELD_NF(ctx->fields + field_index));
    }
    else if (type == CA_FIELD_TYPE_MPOLY_Q)
    {
        fmpz_mpoly_q_neg(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_FIELD_MCTX(ctx->fields + field_index));
    }
}

