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
ca_clear(ca_t x, ca_ctx_t ctx)
{
    slong index;

    index = x->field & ~CA_SPECIAL;

    if (index == CA_FIELD_ID_QQ)
    {
        fmpz_clear(CA_FMPQ_NUMREF(x));
        fmpz_clear(CA_FMPQ_DENREF(x));
    }
    else if (ctx->fields[index].type == CA_FIELD_TYPE_NF)
    {
        nf_elem_clear(CA_NF_ELEM(x), CA_FIELD_NF(ctx->fields + index));
    }
    else if (ctx->fields[index].type == CA_FIELD_TYPE_MULTI ||
             ctx->fields[index].type == CA_FIELD_TYPE_FUNC)
    {
        fmpz_mpoly_q_clear(CA_MPOLY_Q(x), CA_FIELD_MCTX(ctx->fields + index, ctx));
    }
}

