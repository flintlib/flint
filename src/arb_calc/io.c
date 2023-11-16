/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "arb_calc.h"

/* printing *******************************************************************/

void
arf_interval_fprintd(FILE * file, const arf_interval_t v, slong n)
{
    flint_fprintf(file, "[");
    arf_fprintd(file, &v->a, n);
    flint_fprintf(file, ", ");
    arf_fprintd(file, &v->b, n);
    flint_fprintf(file, "]");
}

void arf_interval_printd(const arf_interval_t v, slong n) { arf_interval_fprintd(stdout, v, n); }
