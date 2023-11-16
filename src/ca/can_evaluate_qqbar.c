/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

int
ca_ext_can_evaluate_qqbar(const ca_ext_t x, ca_ctx_t ctx)
{
    if (CA_EXT_IS_QQBAR(x))
        return 1;

    if (CA_EXT_HEAD(x) == CA_Sqrt ||
        CA_EXT_HEAD(x) == CA_Sign ||
        CA_EXT_HEAD(x) == CA_Abs ||
        CA_EXT_HEAD(x) == CA_Re ||
        CA_EXT_HEAD(x) == CA_Im ||
        CA_EXT_HEAD(x) == CA_Conjugate ||
        CA_EXT_HEAD(x) == CA_Floor ||
        CA_EXT_HEAD(x) == CA_Ceil)
        return ca_can_evaluate_qqbar(CA_EXT_FUNC_ARGS(x), ctx);

    if (CA_EXT_HEAD(x) == CA_Pow)
        return ca_can_evaluate_qqbar(CA_EXT_FUNC_ARGS(x), ctx) &&
               CA_IS_QQ(CA_EXT_FUNC_ARGS(x) + 1, ctx);

    return 0;
}

int
ca_can_evaluate_qqbar(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        return 0;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        return 1;
    }
    else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        return 1;
    }
    else
    {
        slong nvars, i;
        int res;
        int * used;

        nvars = CA_FIELD_LENGTH(CA_FIELD(x, ctx));
        used = flint_calloc(nvars, sizeof(int));

        fmpz_mpoly_q_used_vars(used, CA_MPOLY_Q(x), CA_FIELD_MCTX(CA_FIELD(x, ctx), ctx));
        res = 1;

        /* todo: exclude extension numbers that are not actually used */
        for (i = 0; i < nvars; i++)
        {
            if (used[i])
            {
                if (!ca_ext_can_evaluate_qqbar(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), i), ctx))
                {
                    res = 0;
                    break;
                }
            }
        }

        flint_free(used);

        return res;
    }
}
