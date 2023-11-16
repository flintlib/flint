/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_frexp(arf_t man, fmpz_t exp, const arf_t x)
{
    arf_set(man, x);
    fmpz_zero(exp);

    if (!arf_is_special(man))
        fmpz_swap(exp, ARF_EXPREF(man));
}

