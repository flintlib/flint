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
ca_i(ca_t x, ca_ctx_t ctx)
{
    _ca_make_field_element(x, ctx->field_qq_i, ctx);
    nf_elem_gen(CA_NF_ELEM(x), CA_FIELD_NF(ctx->field_qq_i));
}
