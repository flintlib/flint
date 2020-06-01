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
ca_field_clear(ca_field_t K)
{
    slong i;

    if (K->type == CA_FIELD_TYPE_QQ)
        return;

    if (K->type == CA_FIELD_TYPE_NF)
    {
        qqbar_clear(&K->data.nf.x);
        nf_clear(&K->data.nf.nf);
    }

    if (K->type == CA_FIELD_TYPE_MULTI)
    {
        flint_free(K->data.multi.ext);

        if (K->data.multi.ideal_len != 0)
        {
            /* todo: retrieve cached mctx from ctx! */
            fmpz_mpoly_ctx_t mctx;
            fmpz_mpoly_ctx_init(mctx, K->data.multi.len, CA_MPOLY_ORD);

            for (i = 0; i < K->data.multi.ideal_len; i++)
                fmpz_mpoly_clear(K->data.multi.ideal + i, mctx);
        }

        flint_free(K->data.multi.ideal);
    }
}

