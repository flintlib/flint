/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_field.h"

void ca_clear_unchecked(ca_t x, ca_ctx_t ctx);

/* todo: recycle storage when compatible */
void
_ca_make_field_element(ca_t x, ca_field_srcptr field, ca_ctx_t ctx)
{
    ca_field_srcptr old_field;

    if (field == NULL)
    {
        flint_throw(FLINT_ERROR, "NULL in _ca_make_field_element\n");
    }

    old_field = (ca_field_srcptr) (x->field & ~CA_SPECIAL);

    if (old_field == field)
    {
        x->field = (ulong) field;  /* unset special status */
        return;
    }

    ca_clear_unchecked(x, ctx);

    if (field == ctx->field_qq)
    {
        *CA_FMPQ_NUMREF(x) = 0;
        *CA_FMPQ_DENREF(x) = 1;
    }
    else if (CA_FIELD_IS_NF(field))
    {
        nf_elem_init(CA_NF_ELEM(x), CA_FIELD_NF(field));
    }
    else
    {
        x->elem.mpoly_q = (fmpz_mpoly_q_struct *) flint_malloc(sizeof(fmpz_mpoly_q_struct));
        fmpz_mpoly_q_init(CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx));
    }

    x->field = (ulong) field;
}
