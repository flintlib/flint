/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"

void
ca_ext_print(const ca_ext_t x, ca_ctx_t ctx)
{
    if (x->head == CA_QQBar)
    {
        flint_printf("Algebraic ");

        if (qqbar_is_i(CA_EXT_QQBAR(x)))
            flint_printf("I");
        else
        {
/*
            flint_printf("Algebraic [deg %wd] ", qqbar_degree(CA_EXT_QQBAR(x)));
            qqbar_printn(CA_EXT_QQBAR(x), 10);
*/

            qqbar_printn(CA_EXT_QQBAR(x), 8);
/*
            flint_printf(" (");
            fmpz_poly_print_pretty(QQBAR_POLY(CA_EXT_QQBAR(x)), "a");
            flint_printf("=0)");
*/
        }
    }
    else
    {
        flint_printf("%s", calcium_func_name(CA_EXT_HEAD(x)));

        if (CA_EXT_FUNC_NARGS(x) != 0)
        {
            slong i;
            flint_printf("(");
            for (i = 0; i < CA_EXT_FUNC_NARGS(x); i++)
            {
                ca_print(CA_EXT_FUNC_ARGS(x) + i, ctx);

                if (i < CA_EXT_FUNC_NARGS(x) - 1)
                    flint_printf(", ");
            }
            flint_printf(")");
        }
    }
}

