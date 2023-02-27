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
ca_printn(const ca_t x, slong n, ca_ctx_t ctx)
{
    ulong save_flags;
    save_flags = ctx->options[CA_OPT_PRINT_FLAGS];
    ctx->options[CA_OPT_PRINT_FLAGS] = CA_PRINT_N | (CA_PRINT_DIGITS * n);
    ca_print(x, ctx);
    ctx->options[CA_OPT_PRINT_FLAGS] = save_flags;
}
