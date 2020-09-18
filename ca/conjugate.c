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
ca_conjugate(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        /* todo */
        ca_unknown(res, ctx);
    }
    else if (CA_IS_QQ(x, ctx))
    {
        ca_set(res, x, ctx);
    }
    else if (CA_IS_QQ_I(x, ctx))
    {
        ca_set(res, x, ctx);
        fmpz_neg(QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1, QNF_ELEM_NUMREF(CA_NF_ELEM(res)) + 1);
    }
    else if (ca_check_is_real(x, ctx) == T_TRUE)   /* todo: avoid duplicate computations with is_real/is_imaginary */
    {
        ca_set(res, x, ctx);
    }
    else if (ca_check_is_imaginary(x, ctx) == T_TRUE)
    {
        ca_neg(res, x, ctx);
    }
    else
    {
        _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_Conjugate, x), ctx);
        fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
    }
}
