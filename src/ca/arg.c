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
ca_arg(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_SIGNED_INF(x))
        {
            ca_sgn(res, x, ctx);
            ca_arg(res, res, ctx);
        }
        else if (CA_IS_UNKNOWN(x))
        {
            ca_unknown(res, ctx);
        }
        else
        {
            ca_undefined(res, ctx);
        }
    }
    else if (CA_IS_QQ(x, ctx))
    {
        if (fmpz_sgn(CA_FMPQ_NUMREF(x)) >= 0)
        {
            ca_zero(res, ctx);
        }
        else
        {
            ca_pi(res, ctx);
            ca_neg(res, res, ctx);
        }
    }
    else
    {
        ca_t s;
        qqbar_t t;
        slong p;
        ulong q;

        ca_init(s, ctx);
        qqbar_init(t);

        ca_sgn(s, x, ctx);

        if (ca_get_qqbar(t, s, ctx) && qqbar_log_pi_i(&p, &q, t))
        {
            ca_pi(res, ctx);
            ca_mul_si(res, res, p, ctx);
            ca_div_ui(res, res, q, ctx);
        }
        else
        {
            _ca_function_fx(res, CA_Arg, x, ctx);
        }

        ca_clear(s, ctx);
        qqbar_clear(t);
    }
}
