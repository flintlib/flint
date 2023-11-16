/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

ca_ext_ptr
ca_is_gen_as_ext(const ca_t x, ca_ctx_t ctx)
{
    ca_field_ptr K;

    if (CA_IS_SPECIAL(x))
        return NULL;

    K = CA_FIELD(x, ctx);

    if (CA_FIELD_IS_QQ(K))
        return NULL;

    if (CA_FIELD_IS_NF(K))
    {
        if (!nf_elem_is_gen(CA_NF_ELEM(x), CA_FIELD_NF(K)))
            return NULL;

        return CA_FIELD_EXT_ELEM(K, 0);
    }

    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(CA_MPOLY_Q(x)), CA_FIELD_MCTX(K, ctx)))
    {
        if (fmpz_mpoly_is_gen(fmpz_mpoly_q_numref(CA_MPOLY_Q(x)), -1, CA_FIELD_MCTX(K, ctx)))
        {
            slong i;

            for (i = 0; ; i++)
                if (fmpz_mpoly_is_gen(fmpz_mpoly_q_numref(CA_MPOLY_Q(x)), i, CA_FIELD_MCTX(K, ctx)))
                    return CA_FIELD_EXT_ELEM(K, i);
        }
    }

    return NULL;
}
