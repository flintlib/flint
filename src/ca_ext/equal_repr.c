/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"

int
ca_ext_equal_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx)
{
    slong i, nargs;

    if (CA_EXT_HASH(x) != CA_EXT_HASH(y))
        return 0;

    if (CA_EXT_HEAD(x) != CA_EXT_HEAD(y))
        return 0;

    if (CA_EXT_HEAD(x) == CA_QQBar)
        return qqbar_equal(CA_EXT_QQBAR(x), CA_EXT_QQBAR(y));

    nargs = CA_EXT_FUNC_NARGS(x);
    if (nargs != CA_EXT_FUNC_NARGS(y))
        return 0;

    for (i = 0; i < nargs; i++)
        if (!ca_equal_repr(CA_EXT_FUNC_ARGS(x) + i, CA_EXT_FUNC_ARGS(y) + i, ctx))
            return 0;

    return 1;
}
