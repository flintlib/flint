/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"
#include "ca.h"

void _ca_function_fx(ca_t res, calcium_func_code func, const ca_t x, ca_ctx_t ctx)
{
    ca_field_srcptr field = _ca_ctx_get_field_fx(ctx, func, x);
    _ca_make_field_element(res, field, ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_FIELD_MCTX(field, ctx));
}

void _ca_function_fxy(ca_t res, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_field_srcptr field = _ca_ctx_get_field_fxy(ctx, func, x, y);
    _ca_make_field_element(res, field, ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_FIELD_MCTX(field, ctx));
}
