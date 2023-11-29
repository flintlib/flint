/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

void
dlog_precomp_clear(dlog_precomp_t pre)
{
    if (pre == NULL)
        return;
    switch (pre->type)
    {
        case DLOG_MODPE:
            dlog_modpe_clear(pre->t.modpe);
            break;
        case DLOG_CRT:
            dlog_crt_clear(pre->t.crt);
            break;
        case DLOG_POWER:
            dlog_power_clear(pre->t.power);
            break;
        case DLOG_TABLE:
            dlog_table_clear(pre->t.table);
            break;
        case DLOG_BSGS:
            dlog_bsgs_clear(pre->t.bsgs);
            break;
        case DLOG_23:
            dlog_order23_clear(pre->t.order23);
            break;
        default:
            flint_throw(FLINT_ERROR, "dlog_precomp_clear: unknown type %d\n", pre->type);
    }
}
