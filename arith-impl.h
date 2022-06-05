/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARITH_IMPL_H
#define ARITH_IMPL_H

#include <math.h>
#include "fmpz_factor.h"
#include "arith.h"

/* defined in bell_number_dobinski.c */
void _fmpz_ui_pow_ui(fmpz_t x, ulong b, ulong e);

/* defined in bell_number_multi_mod.c */
void divisor_table(unsigned int * tab, slong len);

#endif
