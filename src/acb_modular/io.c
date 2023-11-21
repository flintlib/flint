/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "acb_modular.h"

/* printing *******************************************************************/

void psl2z_fprint(FILE * file, const psl2z_t g)
{
    flint_fprintf(file, "[");
    fmpz_fprint(file, &g->a); flint_fprintf(file, " ");
    fmpz_fprint(file, &g->b); flint_fprintf(file, "; ");
    fmpz_fprint(file, &g->c); flint_fprintf(file, " ");
    fmpz_fprint(file, &g->d); flint_fprintf(file, "]");
}

void psl2z_print(const psl2z_t g) { psl2z_fprint(stdout, g); }
