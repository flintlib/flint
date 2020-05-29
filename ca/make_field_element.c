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
_ca_make_field_element(ca_t x, slong new_index, ca_ctx_t ctx)
{
    slong old_index;
    ca_field_type_t old_type, new_type;

    old_index = x->field;

    if (old_index == new_index)
        return;

    if (new_index < 0 || new_index >= ctx->fields_len)
    {
        flint_printf("_ca_make_field_element: field index out of range\n");
        flint_abort();
    }

    old_type = ctx->fields[old_index].type;
    new_type = ctx->fields[new_index].type;

    /* todo: recycle storage when compatible */
    old_type = old_type;
    ca_clear(x, ctx);

    if (new_type == CA_FIELD_TYPE_QQ)
    {
        *CA_FMPQ_NUMREF(x) = 0;
        *CA_FMPQ_DENREF(x) = 1;
    }
    else if (new_type == CA_FIELD_TYPE_NF)
    {
        nf_elem_init(CA_NF_ELEM(x), &(ctx->fields[new_index].nf_ext->data.qqbar.nf));
    }
    else if (new_type == CA_FIELD_TYPE_MPOLY_Q)
    {
        fmpz_mpoly_q_init(CA_MPOLY_Q(x), CA_FIELD_MCTX(ctx->fields + new_index));
    }

    x->field = new_index;
}

