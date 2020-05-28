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

    if (K->type == CA_FIELD_QQ || K->type == CA_FIELD_NF)
        return;

    if (K->type == CA_FIELD_MPOLY_Q)
    {
        if (K->len != 0)
        {
            flint_free(K->ext);

            for (i = 0; i < K->ideal_len; i++)
                fmpz_mpoly_clear(K->ideal + i, CA_FIELD_MCTX(K));

            flint_free(K->ideal);
        }
    }
}

