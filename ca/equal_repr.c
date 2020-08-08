/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

int
ca_equal_repr(const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    slong field_index;
    ca_field_srcptr field;

    /* by assumption: cached field objects are unique */
    if (x->field != y->field)
        return 0;

    field_index = x->field & ~CA_SPECIAL;
    field = ctx->fields + field_index;

    if (field_index == CA_FIELD_ID_QQ)
    {
        return fmpq_equal(CA_FMPQ(x), CA_FMPQ(y));
    }
    else if (CA_FIELD_IS_NF(field))
    {
        return nf_elem_equal(CA_NF_ELEM(x), CA_NF_ELEM(y), CA_FIELD_NF(field));
    }
    else
    {
        return fmpz_mpoly_q_equal(CA_MPOLY_Q(x), CA_MPOLY_Q(y), CA_FIELD_MCTX(field, ctx));
    }
}
