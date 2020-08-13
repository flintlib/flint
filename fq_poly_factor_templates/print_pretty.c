/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
void
TEMPLATE(T, poly_factor_print_pretty) (const TEMPLATE(T, poly_factor_t) fac,
                                       const char *var,
                                       const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = 0; i < fac->num; i++)
    {
        TEMPLATE(T, poly_print_pretty) (fac->poly + i, var, ctx);
        flint_printf(" ^ %wd\n", fac->exp[i]);
    }
}


#endif
