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
ca_extension_print(const ca_extension_t ext)
{
    if (ext->type == CA_EXT_QQBAR)
    {
        flint_printf("Algebraic number ");
        qqbar_printnd(&ext->data.qqbar.x, 10);
    }
    else if (ext->type == CA_EXT_FUNCTION)
    {
        switch (ext->data.function.func)
        {
            case CA_Pi:
                flint_printf("Pi");
                break;
            default:
                flint_printf("<unknown extension>");
        }
    }
    else
    {
        flint_printf("<unknown extension>");
    }
}

