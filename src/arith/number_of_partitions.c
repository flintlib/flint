/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"
#include "partitions.h"

/* compatibility wrapper */
void
arith_number_of_partitions(fmpz_t x, ulong n)
{
    partitions_fmpz_ui(x, n);
}
