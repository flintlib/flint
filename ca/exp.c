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
ca_exp(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_ext_ptr ext;

    if (CA_IS_SPECIAL(x))
    {
        /* todo: complex signed infinity -> 0, undefined, uinf */

        if (ca_check_is_pos_inf(x, ctx) == T_TRUE)
            ca_pos_inf(res, ctx);
        else if (ca_check_is_neg_inf(x, ctx) == T_TRUE)
            ca_zero(res, ctx);
        else if (ca_check_is_undefined(x, ctx) == T_TRUE || ca_check_is_uinf(x, ctx) == T_TRUE)
            ca_undefined(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }

    ext = ca_is_gen_as_ext(x, ctx);

    /* exp(log(z)) = z */
    if (ext != NULL && CA_EXT_HEAD(ext) == CA_Log)
    {
        ca_set(res, CA_EXT_FUNC_ARGS(ext), ctx);
        return;
    }

    /* todo: evaluate at simple rational multiples of pi i, etc. */
    if (ca_check_is_zero(x, ctx) == T_TRUE)
    {
        ca_one(res, ctx);
        return;
    }

    _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_Exp, x), ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
    /* todo: detect simple values instead of creating extension element in the first place */
    _ca_mpoly_q_reduce_ideal(CA_MPOLY_Q(res), CA_FIELD(res, ctx), ctx);
    ca_condense_field(res, ctx);
}

