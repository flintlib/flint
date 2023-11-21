/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"

#define FEXPR_IS_ARITHMETIC_OP(h) \
    ((h) == FEXPR_SYMBOL_Add || (h) == FEXPR_SYMBOL_Sub || \
     (h) == FEXPR_SYMBOL_Mul || (h) == FEXPR_SYMBOL_Div || \
     (h) == FEXPR_SYMBOL_Neg || (h) == FEXPR_SYMBOL_Pos)

int
fexpr_is_arithmetic_operation(const fexpr_t expr)
{
    ulong head;
    fexpr_t func;
    head = expr->data[0];

    switch (FEXPR_TYPE(head))
    {
        case FEXPR_TYPE_CALL0:
        case FEXPR_TYPE_CALL1:
        case FEXPR_TYPE_CALL2:
        case FEXPR_TYPE_CALL3:
        case FEXPR_TYPE_CALL4:
            return FEXPR_IS_ARITHMETIC_OP(expr->data[1]);
        case FEXPR_TYPE_CALLN:
            fexpr_view_func(func, expr);
            return FEXPR_IS_ARITHMETIC_OP(func->data[0]);
        default:
            return 0;
    }
}
