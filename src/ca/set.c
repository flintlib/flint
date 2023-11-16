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
ca_set(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (res != x)
    {
        ulong field_flags;
        ca_field_srcptr field;

        field_flags = x->field;
        field = (ca_field_srcptr) (x->field & ~CA_SPECIAL);

        /* for Undefined, Unknown, UnsignedInfinity */
        if (field == NULL)
        {
            ca_clear(res, ctx);
            res->field = field_flags;
            return;
        }

        _ca_make_field_element(res, field, ctx);
        res->field = field_flags;  /* set special flags */

        if (field != NULL)
        {
            if (CA_FIELD_IS_QQ(field))
            {
                fmpq_set(CA_FMPQ(res), CA_FMPQ(x));
            }
            else if (CA_FIELD_IS_NF(field))
            {
                nf_elem_set(CA_NF_ELEM(res), CA_NF_ELEM(x), CA_FIELD_NF(field));
            }
            else
            {
                fmpz_mpoly_q_set(CA_MPOLY_Q(res), CA_MPOLY_Q(x), CA_FIELD_MCTX(field, ctx));
            }
        }
    }
}

