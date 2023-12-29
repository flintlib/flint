/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

void fmpzi_print(const fmpzi_t x)
{
    fmpz_print(fmpzi_realref(x));
    if (fmpz_sgn(fmpzi_imagref(x)) >= 0)
        flint_printf("+");
    fmpz_print(fmpzi_imagref(x));
    flint_printf("*I");
}
