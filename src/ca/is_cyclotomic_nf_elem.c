/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

int
ca_is_cyclotomic_nf_elem(slong * p, ulong * q, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
        return 0;

    if (CA_IS_QQ(x, ctx))
        return 0;

    if (CA_IS_QQ_I(x, ctx))
    {
        if (p != NULL) p[0] = 1;
        if (q != NULL) q[0] = 4;
        return 1;
    }

    return CA_FIELD_IS_NF(CA_FIELD(x, ctx)) &&
        qqbar_is_root_of_unity(p, q, CA_FIELD_NF_QQBAR(CA_FIELD(x, ctx)));
}
