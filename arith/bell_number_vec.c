/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void
arith_bell_number_vec(fmpz * res, slong n)
{
    if (n < 5000)
        arith_bell_number_vec_recursive(res, n);
    else
        arith_bell_number_vec_multi_mod(res, n);
}
