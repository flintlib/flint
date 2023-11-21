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
ca_sgn(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_SIGNED_INF(x))
        {
            ca_set(res, x, ctx);
            res->field &= ~CA_INF;
        }
        else if (CA_IS_UNKNOWN(x))
        {
            ca_unknown(res, ctx);
        }
        else
        {
            ca_undefined(res, ctx);
        }
        return;
    }
    else if (CA_IS_QQ(x, ctx))
    {
        ca_set_si(res, fmpz_sgn(fmpq_numref(CA_FMPQ(x))), ctx);
    }
    else
    {
        qqbar_t t;
        qqbar_init(t);

        if (ca_get_qqbar(t, x, ctx))
        {
            qqbar_sgn(t, t);

            if (qqbar_within_limits(t, ctx->options[CA_OPT_QQBAR_DEG_LIMIT], 0))
                ca_set_qqbar(res, t, ctx);
            else
                _ca_function_fx(res, CA_Sign, x, ctx);
        }
        else
        {
            _ca_function_fx(res, CA_Sign, x, ctx);
        }

        qqbar_clear(t);
    }
}

