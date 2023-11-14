/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nf.h"

void nf_print(const nf_t nf)
{
    flint_printf("Number field with defining polynomial ");
    fmpq_poly_print_pretty(nf->pol, "x");
}

