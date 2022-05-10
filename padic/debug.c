/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "flint-impl.h"
#include "fmpz.h"
#include "padic.h"

void
padic_debug(const padic_t op)
{
    printf("(");
    fmpz_print(padic_unit(op)); 
    printf(" " WORD_FMT "d " WORD_FMT "d)", padic_val(op), padic_prec(op));
}
