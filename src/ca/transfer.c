/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_transfer(ca_t res, ca_ctx_t res_ctx, const ca_t src, ca_ctx_t src_ctx)
{
    if (res_ctx == src_ctx)
    {
        ca_set(res, src, res_ctx);
    }
    else if (CA_IS_QQ(src, src_ctx))
    {
        _ca_make_fmpq(res, res_ctx);
        fmpq_set(CA_FMPQ(res), CA_FMPQ(src));
    }
    else
    {
        fexpr_t expr;
        fexpr_init(expr);

        /* todo: optimizations, e.g. direct transfer of number field
           elements where permissible */
        ca_get_fexpr(expr, src, CA_FEXPR_SERIALIZATION, src_ctx);

        if (!ca_set_fexpr(res, expr, res_ctx))
        {
            flint_throw(FLINT_ERROR, "ca_transfer: failed to recreate from expression!\n");
        }

        fexpr_clear(expr);
    }
}
