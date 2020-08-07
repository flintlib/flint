/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_field.h"

void
ca_field_clear(ca_field_t K, ca_ctx_t ctx)
{
    slong i, length, ideal_length;

    length = CA_FIELD_LENGTH(K);

    if (length == 0)
        return;

    flint_free(CA_FIELD_EXT(K));

    if (CA_FIELD_IS_NF(K))
        return;

    ideal_length = CA_FIELD_IDEAL_LENGTH(K);

    if (ideal_length > 0)
    {
        for (i = 0; i < ideal_length; i++)
            fmpz_mpoly_clear(CA_FIELD_IDEAL_ELEM(K, i), CA_FIELD_MCTX(K, ctx));

        flint_free(CA_FIELD_IDEAL(K));
    }
}

