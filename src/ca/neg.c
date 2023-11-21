/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_neg(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_field_srcptr field;
    ulong field_flags;

    if (CA_IS_QQ(x, ctx))
    {
        _ca_make_fmpq(res, ctx);
        fmpq_neg(CA_FMPQ(res), CA_FMPQ(x));
        return;
    }

    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_UNKNOWN(x) || CA_IS_UNDEFINED(x) || CA_IS_UNSIGNED_INF(x))
        {
            ca_set(res, x, ctx);
            return;
        }
    }

    field_flags = x->field;
    _ca_make_field_element(res, (ca_field_ptr) (field_flags & ~CA_SPECIAL), ctx);
    field = CA_FIELD(res, ctx);
    res->field = field_flags;  /* set special flags */

    if (field == ctx->field_qq)
    {
        fmpq_neg(CA_FMPQ(res), CA_FMPQ(x));
    }
    else if (CA_FIELD_IS_NF(field))
    {
        nf_elem_neg(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_FIELD_NF(field));
    }
    else
    {
        fmpz_mpoly_q_neg(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx));
    }
}

