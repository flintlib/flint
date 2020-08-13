/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_init(padic_t rop)
{
    fmpz_init(padic_unit(rop));
    padic_val(rop)  = 0;
    padic_prec(rop) = PADIC_DEFAULT_PREC;
}

void padic_init2(padic_t rop, slong N)
{
    fmpz_init(padic_unit(rop));
    padic_val(rop)  = 0;
    padic_prec(rop) = N;
}

