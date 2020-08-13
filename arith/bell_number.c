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
arith_bell_number(fmpz_t b, ulong n)
{
    if (n < BELL_NUMBER_TAB_SIZE)
        fmpz_set_ui(b, bell_number_tab[n]);
    else if (n < 5000)
        arith_bell_number_bsplit(b, n);
    else
        arith_bell_number_multi_mod(b, n);
}
