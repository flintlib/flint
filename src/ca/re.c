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
ca_re(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (CA_IS_UNKNOWN(x))
            ca_unknown(res, ctx);
        else
            ca_undefined(res, ctx);
    }
    else if (CA_IS_QQ(x, ctx))
    {
        ca_set(res, x, ctx);
    }
    else if (CA_IS_QQ_I(x, ctx))
    {
        const fmpz *n, *d;
        fmpq_t t;

        n = QNF_ELEM_NUMREF(CA_NF_ELEM(x));
        d = QNF_ELEM_DENREF(CA_NF_ELEM(x));

        fmpq_init(t);
        fmpq_set_fmpz_frac(t, n, d);
        ca_set_fmpq(res, t, ctx);
        fmpq_clear(t);
    }
    else if (ca_check_is_real(x, ctx) == T_TRUE)   /* todo: avoid duplicate computations with is_real/is_imaginary */
    {
        ca_set(res, x, ctx);
    }
    else if (ca_check_is_imaginary(x, ctx) == T_TRUE)
    {
        ca_zero(res, ctx);
    }
    else
    {
        _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_Re, x), ctx);
        fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
    }
}
