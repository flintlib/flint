/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_field.h"

void
ca_field_clear(ca_field_t K, ca_ctx_t ctx)
{
    slong length;

    length = CA_FIELD_LENGTH(K);

    if (length == 0)
        return;

    flint_free(CA_FIELD_EXT(K));

    if (CA_FIELD_IS_NF(K))
        return;

    fmpz_mpoly_vec_clear(CA_FIELD_IDEAL(K), CA_FIELD_MCTX(K, ctx));
}

