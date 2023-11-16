/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

int
fexpr_is_builtin_call(const fexpr_t expr, slong i)
{
    fexpr_t func;

    if (fexpr_is_atom(expr))
        return 0;

    fexpr_view_func(func, expr);

    return fexpr_is_builtin_symbol(func, i);
}
