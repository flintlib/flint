/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

ulong
dlog_precomp(const dlog_precomp_t pre, ulong b)
{
    if (b == 1)
        return 0;

    switch (pre->type)
    {
        case DLOG_MODPE:
            return dlog_modpe(pre->t.modpe, b);
        case DLOG_CRT:
            return dlog_crt(pre->t.crt, b);
        case DLOG_POWER:
            return dlog_power(pre->t.power, b);
        case DLOG_TABLE:
            return dlog_table(pre->t.table, b);
        case DLOG_BSGS:
            return dlog_bsgs(pre->t.bsgs, b);
        case DLOG_23:
            return dlog_order23(pre->t.order23, b);
        default:
            flint_throw(FLINT_ERROR, "(dlog_precomp): Unknown option");
    }
}
