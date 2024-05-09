/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"

void ca_ext_get_arg(ca_t res, const ca_ext_t x, slong i, ca_ctx_t ctx)
{
    if (CA_EXT_HEAD(x) == CA_QQBar || i < 0 || i >= CA_EXT_FUNC_NARGS(x))
    {
        flint_throw(FLINT_ERROR, "ca_ext_get_arg: index out of range\n");
    }
    else
    {
        ca_set(res, CA_EXT_FUNC_ARGS(x) + i, ctx);
    }
}
