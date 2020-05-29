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
ca_i(ca_t x, ca_ctx_t ctx)
{
    _ca_make_field_element(x, CA_FIELD_ID_QQ_I, ctx);
    nf_elem_gen(CA_NF_ELEM(x), &(ctx->fields[CA_FIELD_ID_QQ_I].nf_ext->data.qqbar.nf));
}

