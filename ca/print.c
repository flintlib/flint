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
ca_print(ca_t x, ca_ctx_t ctx)
{
    slong field;

    field = x->field;

    if (field == CA_FIELD_ID_QQ)
    {
        fmpq_print(CA_FMPQ(x));
    }
    else if (ctx->fields[field].type == CA_FIELD_TYPE_NF)
    {
        nf_elem_print_pretty(CA_NF_ELEM(x), &(ctx->fields[field].nf_ext->data.qqbar.nf), "x1");
    }
    else
    {
        fmpz_mpoly_q_print_pretty(CA_MPOLY_Q(x), NULL, CA_FIELD_MCTX(ctx->fields + field));
    }

    flint_printf("  in  ");
    ca_field_print(ctx->fields + field);
}

