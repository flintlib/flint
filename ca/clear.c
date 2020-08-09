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
    ca_field_srcptr field;

    field = (ca_field_srcptr) (x->field & ~CA_SPECIAL);

    if (field != NULL)
    {
        if (field == ctx->field_qq)
        {
            fmpz_clear(CA_FMPQ_NUMREF(x));
            fmpz_clear(CA_FMPQ_DENREF(x));
        }
        else if (CA_FIELD_IS_NF(field))
        {
            nf_elem_clear(CA_NF_ELEM(x), CA_FIELD_NF(field));
        }
        else
        {
            fmpz_mpoly_q_clear(CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx));
            flint_free(x->elem.mpoly_q);
        }
    }
}

