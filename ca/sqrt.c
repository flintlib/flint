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
ca_sqrt_inert(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        /* todo: could compute inert sign + sqrt for signed inf */
        ca_sqrt(res, x, ctx);
    }
    else
    {
        _ca_function_fx(res, CA_Sqrt, x, ctx);
    }
}

void
_ca_sqrt_nofactor(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_t y, tmp;
    ca_init(y, ctx);
    ca_init(tmp, ctx);

    _ca_function_fx(y, CA_Sqrt, x, ctx);
    ca_merge_fields(tmp, res, x, y, ctx);

    ca_clear(y, ctx);
    ca_clear(tmp, ctx);
}

void
ca_sqrt_nofactor(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_SIGNED_INF(x))
        {
            ca_sgn(res, x, ctx);
            ca_sqrt(res, res, ctx);
            if (!ca_is_unknown(res, ctx))
                res->field |= CA_INF;
        }
        else
        {
            ca_set(res, x, ctx);
        }
    }
    else
    {
        qqbar_t t;
        slong deg;

        qqbar_init(t);

        if (ca_get_qqbar(t, x, ctx))
        {
            deg = qqbar_degree(t);

            qqbar_sqrt(t, t);

            /* use? qqbar_within_limits(t, ctx->options[CA_OPT_QQBAR_DEG_LIMIT], 0) */

            if (qqbar_degree(t) <= FLINT_MAX(2, deg))
                ca_set_qqbar(res, t, ctx);
            else
                _ca_sqrt_nofactor(res, x, ctx);
        }
        else
        {
            _ca_sqrt_nofactor(res, x, ctx);
        }

        qqbar_clear(t);
    }
}

void
ca_sqrt(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_sqrt_nofactor(res, x, ctx);
}

