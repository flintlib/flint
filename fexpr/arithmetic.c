/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"

void
fexpr_neg(fexpr_t res, const fexpr_t a)
{
    /* todo: handle aliasing in call1 */
    if (res == a)
    {
        fexpr_t tmp;
        fexpr_init(tmp);
        fexpr_set(tmp, a);
        fexpr_neg(res, tmp);
        fexpr_clear(tmp);
    }
    else
    {
        fexpr_t tmp;
        ulong tmp_head = FEXPR_SYMBOL_Neg;
        tmp->data = &tmp_head;
        tmp->alloc = 1;
        fexpr_call1(res, tmp, a);
    }
}

void
fexpr_add(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Add;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

void
fexpr_sub(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Sub;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

void
fexpr_mul(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Mul;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

void
fexpr_div(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Div;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

void
fexpr_pow(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Pow;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

int
fexpr_is_Pow(const fexpr_t expr)
{
    fexpr_t op;
    ulong op_head;

    if (fexpr_is_atom(expr))
        return 0;

    fexpr_view_func(op, expr);
    op_head = op->data[0];

    return op_head == FEXPR_SYMBOL_Pow;
}
