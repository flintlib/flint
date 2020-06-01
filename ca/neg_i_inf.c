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
ca_neg_i_inf(ca_t x, ca_ctx_t ctx)
{
    ca_i(x, ctx);
    nf_elem_neg(CA_NF_ELEM(x), CA_NF_ELEM(x), CA_FIELD_NF(ctx->fields + CA_FIELD_ID_QQ_I));
    x->field |= CA_SIGNED_INF;
}

