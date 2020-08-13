/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"

double
arith_bell_number_size(ulong n)
{
    if (n == 0)
        return 2;

    return n * log(0.792 * n/log(n+1)) * 1.44269504088896 + 2;
}
