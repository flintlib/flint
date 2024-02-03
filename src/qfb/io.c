/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "qfb.h"

/* printing *******************************************************************/

void qfb_print(qfb_t q)
{
    printf("(");
    fmpz_print(q->a); printf(", ");
    fmpz_print(q->b); printf(", ");
    fmpz_print(q->c); printf(")");
}
