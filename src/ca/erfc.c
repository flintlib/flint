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
ca_erfc(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        /* todo: other signed infinities */

        if (ca_check_is_pos_inf(x, ctx) == T_TRUE)
            ca_zero(res, ctx);
        else if (ca_check_is_neg_inf(x, ctx) == T_TRUE)
            ca_set_ui(res, 2, ctx);
        else if (ca_check_is_pos_i_inf(x, ctx) == T_TRUE)
            ca_neg_i_inf(res, ctx);
        else if (ca_check_is_neg_i_inf(x, ctx) == T_TRUE)
            ca_pos_i_inf(res, ctx);
        else if (ca_check_is_undefined(x, ctx) == T_TRUE || ca_check_is_uinf(x, ctx) == T_TRUE)
            ca_undefined(res, ctx);
        else
            ca_unknown(res, ctx);
        return;
    }

    if (ca_check_is_zero(x, ctx) == T_TRUE)
    {
        ca_one(res, ctx);
        return;
    }

    _ca_make_field_element(res, _ca_ctx_get_field_fx(ctx, CA_Erfc, x), ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_MCTX_1(ctx));
}
